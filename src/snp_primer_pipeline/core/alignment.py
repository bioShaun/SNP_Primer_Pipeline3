#!/usr/bin/env python3
"""
Multiple sequence alignment module for SNP Primer Pipeline.

This module handles multiple sequence alignment and variant site identification.
"""

import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

from ..exceptions import AlignmentError


# V2-compatible helper functions
def _tm_simple(seq: str) -> int:
    """Simple Tm calculator (V2 style)."""
    t = 0
    for a in seq.upper():
        if a in 'AT':
            t += 2
        if a in 'CG':
            t += 4
    return t


def _find_longest_substring(s1: str, s2: str) -> tuple:
    """
    Find the segment with largest Tm (V2 style).
    Returns: (longestStart, longestEnd, nL, nR)
    """
    longest_start = 0
    longest_end = 0
    largest_tm = 0
    start = 0
    
    gaps = [i for i, c in enumerate(s1) if c == '-' or s2[i] == '-']
    gaps.append(len(s1))
    
    for gap in gaps:
        end = gap
        tm = _tm_simple(s1[start:end])
        if tm > largest_tm:
            longest_start = start
            longest_end = end
            largest_tm = tm
        start = gap + 1
    
    n_left = len(s1[:longest_start].replace("-", ""))
    n_right = len(s1[longest_end:].replace("-", ""))
    return longest_start, longest_end, n_left, n_right


def _score_pairwise(seq1: str, seq2: str, gapopen: float = -4.0, 
                    gapext: float = -1.0, match: float = 1.0, 
                    mismatch: float = -1.0) -> float:
    """Alignment score (V2 style)."""
    score = 0.0
    gap = False
    
    for i in range(len(seq1)):
        pair = (seq1[i], seq2[i])
        if not gap:
            if '-' in pair:
                gap = True
                score += gapopen
            elif seq1[i] == seq2[i]:
                score += match
            else:
                score += mismatch
        else:
            if '-' not in pair:
                gap = False
                if seq1[i] == seq2[i]:
                    score += match
                else:
                    score += mismatch
            else:
                score += gapext
    return score


def _gap_diff(seq1: str, seq2: str) -> int:
    """Count gap differences for an alignment (V2 style)."""
    ngap = 0
    for i in range(len(seq1)):
        pair = seq1[i] + seq2[i]
        if pair.count("-") == 1:
            ngap += 1
    return ngap


def _get_homeo_seq(fasta: dict, target: str, ids: list, 
                   align_left: int, align_right: int) -> list:
    """
    Get the list of sequences in the homeolog groups for comparison (V2 style).
    
    Args:
        fasta: Dictionary of sequence name -> aligned sequence
        target: Target sequence name
        ids: List of other sequence names
        align_left: Left alignment position
        align_right: Right alignment position
        
    Returns:
        List of sequences for comparison (one per homeolog)
    """
    s1 = fasta[target]
    seq2comp = []
    
    for k in ids:
        s2 = fasta[k]
        target_seq = s1[align_left:(align_right + 1)]
        homeo_seq = s2[align_left:(align_right + 1)]
        score1 = _score_pairwise(target_seq, homeo_seq)
        
        # Get the sequences for comparison
        index_l, index_r, n_left, n_right = _find_longest_substring(target_seq, homeo_seq)
        index_l += align_left
        index_r += align_left
        
        seq_l = s2[:index_l].replace("-", "")
        seq_r = s2[index_r:].replace("-", "")
        
        # Handle case where not enough bases available
        if len(seq_l) < n_left:
            seq_l = "-" * (n_left - len(seq_l)) + seq_l
        if len(seq_r) < n_right:
            seq_r = seq_r + "-" * (n_right - len(seq_r))
        
        seqk = seq_l[::-1][:n_left][::-1] + s2[index_l:index_r] + seq_r[:n_right]
        
        score2 = _score_pairwise(target_seq.replace("-", ""), seqk)
        
        # If score in alignment is better and gaps < 4, use gap-removed version
        if score1 > score2 and _gap_diff(target_seq, homeo_seq) < 4:
            seqk = "".join([homeo_seq[i] for i, c in enumerate(target_seq) if c != '-'])
        
        seq2comp.append(seqk)
    
    return seq2comp


@dataclass

class AlignedSequence:
    """Aligned sequence with gap information."""
    
    name: str
    sequence: str  # Sequence with gaps
    
    @property
    def clean_sequence(self) -> str:
        """Return sequence without gaps."""
        return self.sequence.replace("-", "")
    
    @property
    def length(self) -> int:
        """Return aligned sequence length."""
        return len(self.sequence)


class MultipleSequenceAlignment:
    """Multiple sequence alignment result."""
    
    def __init__(self, sequences: List[AlignedSequence]):
        """
        Initialize alignment.
        
        Args:
            sequences: List of aligned sequences
        """
        self.sequences = sequences
        self._t2a: Dict[int, int] = {}  # template to alignment mapping
        self._a2t: Dict[int, int] = {}  # alignment to template mapping
        self._target_name: Optional[str] = None
        
        if sequences:
            self.alignment_length = sequences[0].length
            # Verify all sequences have same length
            for seq in sequences:
                if seq.length != self.alignment_length:
                    raise AlignmentError(
                        f"Sequence {seq.name} has length {seq.length}, "
                        f"expected {self.alignment_length}"
                    )
    
    def set_target_sequence(self, target_name: str) -> None:
        """
        Set target sequence and build coordinate mappings.
        
        Args:
            target_name: Name of target sequence
            
        Raises:
            AlignmentError: If target sequence not found
        """
        target_seq = None
        for seq in self.sequences:
            if seq.name == target_name:
                target_seq = seq
                break
        
        if target_seq is None:
            raise AlignmentError(f"Target sequence '{target_name}' not found")
        
        self._target_name = target_name
        self._build_coordinate_mappings(target_seq)
    
    def _build_coordinate_mappings(self, target_seq: AlignedSequence) -> None:
        """Build coordinate mappings between template and alignment positions."""
        template_pos = 0
        
        for align_pos, char in enumerate(target_seq.sequence):
            if char != "-":
                # 0-based positions
                self._t2a[template_pos] = align_pos
                self._a2t[align_pos] = template_pos
                template_pos += 1
    
    def template_to_alignment(self, template_pos: int) -> Optional[int]:
        """
        Convert template position to alignment position.
        
        Args:
            template_pos: Position in template sequence (0-based)
            
        Returns:
            Alignment position (0-based) or None if invalid
        """
        return self._t2a.get(template_pos)
    
    def alignment_to_template(self, align_pos: int) -> Optional[int]:
        """
        Convert alignment position to template position.
        
        Args:
            align_pos: Position in alignment (0-based)
            
        Returns:
            Template position (0-based) or None if invalid
        """
        return self._a2t.get(align_pos)
    
    def find_variant_sites(
        self,
        target_name: str,
        min_gap_left: int = 20,
        min_gap_right: int = 20,
    ) -> Tuple[List[int], List[int]]:
        """
        Find variant sites in the alignment using V2 logic.
        
        Args:
            target_name: Name of target sequence
            min_gap_left: Minimum gap-free positions to the left (not used in V2 logic)
            min_gap_right: Minimum gap-free positions to the right (not used in V2 logic)
            
        Returns:
            Tuple of (sites_diff_all, sites_diff_any):
            - sites_diff_all: Template positions where target differs from ALL other sequences
            - sites_diff_any: Template positions where target differs from ANY other sequence
            
            Note: Returns TEMPLATE coordinates (0-based), not alignment coordinates.
            This matches V2 behavior.
            
        Raises:
            AlignmentError: If target sequence not found
        """
        # Set target sequence if not already set
        if self._target_name != target_name:
            self.set_target_sequence(target_name)
        
        target_seq = None
        other_seqs = []
        
        for seq in self.sequences:
            if seq.name == target_name:
                target_seq = seq
            else:
                other_seqs.append(seq)
        
        if target_seq is None:
            raise AlignmentError(f"Target sequence '{target_name}' not found")
        
        if not other_seqs:
            # No other sequences to compare
            return [], []
        
        # Get clean template sequence length for bounds checking
        template_length = len(target_seq.clean_sequence)
        
        # Find gap boundaries similar to V2
        target_gaps = target_seq.sequence
        gap_left = max([len(seq.sequence) - len(seq.sequence.lstrip('-')) for seq in self.sequences])
        gap_right = min([len(seq.sequence.rstrip('-')) for seq in self.sequences])
        
        # Build t2a and a2t mappings for the target sequence
        t2a = {}  # template to alignment mapping
        a2t = {}  # alignment to template mapping
        template_pos = 0
        
        for align_pos, char in enumerate(target_seq.sequence):
            if char != "-":
                t2a[template_pos] = align_pos
                a2t[align_pos] = template_pos
                template_pos += 1
        
        sites_diff_all = []  # Target differs from ALL others (template coords)
        sites_diff_any = []  # Target differs from ANY others (template coords)
        
        # V2-style loop: only check within gap boundaries
        for i in range(gap_left, gap_right):
            b1 = target_seq.sequence[i]  # target base at alignment position i
            
            if b1 == "-":  # Skip gaps in target
                continue
                
            pos_template = a2t.get(i)  # position in the target template (no gaps)
            if pos_template is None:
                continue
                
            # V2-style boundary check
            if pos_template < 20 or pos_template > template_length - 20:
                continue
            
            # Compare with other sequences at the same alignment position
            diff_from_all = True
            diff_from_any = False
            
            for other_seq in other_seqs:
                b2 = other_seq.sequence[i]
                
                # Skip gaps in other sequences
                if b2 == "-":
                    continue
                
                if b1 != b2:
                    diff_from_any = True
                else:
                    diff_from_all = False
            
            if diff_from_any:
                if pos_template not in sites_diff_any:  # avoid duplicates from gaps
                    sites_diff_any.append(pos_template)
            
            if diff_from_all:
                if pos_template not in sites_diff_all:  # avoid duplicates from gaps
                    sites_diff_all.append(pos_template)
        
        return sorted(sites_diff_all), sorted(sites_diff_any)
    
    def find_variant_sites_v2(
        self,
        target_name: str,
        snp_position: int,
    ) -> Tuple[List[int], List[int], Dict[int, List[int]]]:
        """
        Find variant sites using exact V2 logic with get_homeo_seq.
        
        This method implements V2's exact algorithm that uses 20bp window
        comparison around each position relative to the SNP position.
        
        Args:
            target_name: Name of target sequence
            snp_position: SNP position in template (0-based)
            
        Returns:
            Tuple of (sites_diff_all, sites_diff_any, diffarray):
            - sites_diff_all: Template positions differing from ALL homeologs
            - sites_diff_any: Template positions differing from ANY homeolog
            - diffarray: Dict mapping template pos -> list of 0/1 for each homeolog
        """
        # Build fasta dictionary for get_homeo_seq
        fasta = {}
        target = None
        ids = []
        
        for seq in self.sequences:
            fasta[seq.name] = seq.sequence
            if seq.name == target_name:
                target = seq.name
            else:
                ids.append(seq.name)
        
        if target is None:
            raise AlignmentError(f"Target sequence '{target_name}' not found")
        
        if not ids:
            return [], [], {}
        
        target_seq_obj = self.get_sequence_by_name(target_name)
        seq_template = target_seq_obj.clean_sequence
        template_length = len(seq_template)
        
        # Build t2a and a2t mappings
        t2a = {}
        a2t = {}
        template_pos = 0
        for align_pos, char in enumerate(target_seq_obj.sequence):
            if char != "-":
                t2a[template_pos] = align_pos
                a2t[align_pos] = template_pos
                template_pos += 1
        
        # Calculate gap boundaries (V2 style)
        gap_left = max([len(seq.sequence) - len(seq.sequence.lstrip('-')) 
                       for seq in self.sequences])
        gap_right = min([len(seq.sequence.rstrip('-')) for seq in self.sequences])
        
        variation = []  # Sites that differ from ALL
        variation2 = []  # Sites that differ from at least 1
        diffarray = {}
        
        # V2-style loop through alignment positions
        for i in range(gap_left, gap_right):
            b1 = fasta[target][i]
            
            if b1 == "-":
                continue
            
            pos_template = a2t.get(i)
            if pos_template is None:
                continue
            
            # V2 boundary check
            if pos_template < 20 or pos_template > template_length - 20:
                continue
            
            nd = 0  # number of differences
            da = [0] * len(ids)  # differ array
            
            # Calculate alignment range for get_homeo_seq (V2 style)
            if pos_template < snp_position:
                # 20 bp left of current position
                left_pos = pos_template - 19
                if left_pos < 0:
                    left_pos = 0
                if left_pos not in t2a:
                    continue
                align_left = t2a[left_pos]
                align_right = i
            else:
                # 20 bp right of current position
                right_pos = pos_template + 19
                if right_pos >= template_length:
                    right_pos = template_length - 1
                if right_pos not in t2a:
                    continue
                align_left = i
                align_right = t2a[right_pos]
            
            # Get homeolog sequences for comparison
            seq2comp = _get_homeo_seq(fasta, target, ids, align_left, align_right)
            
            # Compare bases (V2 style)
            for m, k in enumerate(seq2comp):
                if pos_template < snp_position:
                    b2 = k[-1] if k else '-'  # Last base of extracted sequence
                else:
                    b2 = k[0] if k else '-'   # First base of extracted sequence
                
                if b1 != b2:
                    nd += 1
                    da[m] = 1
            
            # Store in diffarray
            diffarray[pos_template] = da
            
            # Check if differs from all
            if nd == len(ids):
                if pos_template not in variation:
                    variation.append(pos_template)
            
            # Check if differs from any
            if nd > 0:
                if pos_template not in variation2:
                    variation2.append(pos_template)
        
        return sorted(variation), sorted(variation2), diffarray

    
    def _is_gap_free_region(
        self,
        sequence: AlignedSequence,
        center_pos: int,
        left_margin: int,
        right_margin: int
    ) -> bool:
        """Check if region around position is gap-free."""
        start = max(0, center_pos - left_margin)
        end = min(len(sequence.sequence), center_pos + right_margin + 1)
        
        region = sequence.sequence[start:end]
        return "-" not in region
    
    def get_sequence_by_name(self, name: str) -> Optional[AlignedSequence]:
        """Get sequence by name."""
        for seq in self.sequences:
            if seq.name == name:
                return seq
        return None
    
    def to_fasta(self, output_file: Path) -> Path:
        """
        Write alignment to FASTA file.
        
        Args:
            output_file: Output file path
            
        Returns:
            Path to output file
        """
        try:
            with open(output_file, 'w') as f:
                for seq in self.sequences:
                    f.write(f">{seq.name}\n{seq.sequence}\n")
            return Path(output_file)
        except IOError as e:
            raise AlignmentError(f"Failed to write alignment file: {e}") from e


class MultipleSequenceAligner:
    """Multiple sequence aligner using external tools."""
    
    def __init__(self, algorithm: str = "muscle", software_paths=None):
        """
        Initialize aligner.
        
        Args:
            algorithm: Alignment algorithm ("muscle" or "mafft")
            software_paths: SoftwarePaths object with tool paths
        """
        self.algorithm = algorithm.lower()
        if self.algorithm not in ["muscle", "mafft"]:
            raise AlignmentError(f"Unsupported algorithm: {algorithm}")
        
        # Get software paths
        if software_paths is None:
            from ..config import SoftwarePaths
            self.software_paths = SoftwarePaths.auto_detect()
        else:
            self.software_paths = software_paths
    
    def align(self, sequences: Dict[str, str]) -> MultipleSequenceAlignment:
        """
        Align sequences.
        
        Args:
            sequences: Dictionary mapping sequence names to sequences
            
        Returns:
            MultipleSequenceAlignment object
            
        Raises:
            AlignmentError: If alignment fails
        """
        if len(sequences) < 2:
            raise AlignmentError("At least 2 sequences required for alignment")
        
        # Create temporary input file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_input:
            for name, seq in sequences.items():
                tmp_input.write(f">{name}\n{seq}\n")
            tmp_input_path = Path(tmp_input.name)
        
        try:
            # Create temporary output file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_output:
                tmp_output_path = Path(tmp_output.name)
            
            # Run alignment
            if self.algorithm == "muscle":
                self._run_muscle(tmp_input_path, tmp_output_path)
            elif self.algorithm == "mafft":
                self._run_mafft(tmp_input_path, tmp_output_path)
            
            # Parse result
            return self._parse_alignment_file(tmp_output_path)
            
        finally:
            # Clean up temporary files
            try:
                tmp_input_path.unlink()
                tmp_output_path.unlink()
            except:
                pass
    
    def align_file(self, fasta_file: Path, output_file: Path) -> MultipleSequenceAlignment:
        """
        Align sequences from FASTA file.
        
        Args:
            fasta_file: Input FASTA file
            output_file: Output alignment file
            
        Returns:
            MultipleSequenceAlignment object
        """
        if self.algorithm == "muscle":
            self._run_muscle(fasta_file, output_file)
        elif self.algorithm == "mafft":
            self._run_mafft(fasta_file, output_file)
        
        return self._parse_alignment_file(output_file)
    
    def _run_muscle(self, input_file: Path, output_file: Path) -> None:
        """Run MUSCLE alignment."""
        cmd = [str(self.software_paths.muscle_path), "-in", str(input_file), "-out", str(output_file)]
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
        except subprocess.CalledProcessError as e:
            raise AlignmentError(f"MUSCLE alignment failed: {e.stderr}") from e
        except FileNotFoundError as e:
            raise AlignmentError(f"MUSCLE not found at {self.software_paths.muscle_path}: {e}") from e
    
    def _run_mafft(self, input_file: Path, output_file: Path) -> None:
        """Run MAFFT alignment."""
        cmd = ["mafft", "--auto", str(input_file)]
        
        try:
            with open(output_file, 'w') as f:
                result = subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True
                )
        except subprocess.CalledProcessError as e:
            raise AlignmentError(f"MAFFT alignment failed: {e.stderr}") from e
        except FileNotFoundError as e:
            raise AlignmentError(f"MAFFT not found in PATH: {e}") from e
    
    def _parse_alignment_file(self, alignment_file: Path) -> MultipleSequenceAlignment:
        """Parse alignment file into MultipleSequenceAlignment object."""
        sequences = []
        current_name = None
        current_seq = []
        
        try:
            with open(alignment_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    if line.startswith('>'):
                        # Save previous sequence
                        if current_name is not None:
                            sequences.append(AlignedSequence(
                                name=current_name,
                                sequence=''.join(current_seq)
                            ))
                        
                        # Start new sequence
                        current_name = line[1:]
                        current_seq = []
                    else:
                        current_seq.append(line)
                
                # Save last sequence
                if current_name is not None:
                    sequences.append(AlignedSequence(
                        name=current_name,
                        sequence=''.join(current_seq)
                    ))
            
            if not sequences:
                raise AlignmentError("No sequences found in alignment file")
            
            return MultipleSequenceAlignment(sequences)
            
        except IOError as e:
            raise AlignmentError(f"Failed to read alignment file: {e}") from e