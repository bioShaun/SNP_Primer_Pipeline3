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
        Find variant sites in the alignment.
        
        Args:
            target_name: Name of target sequence
            min_gap_left: Minimum gap-free positions to the left
            min_gap_right: Minimum gap-free positions to the right
            
        Returns:
            Tuple of (sites_diff_all, sites_diff_any):
            - sites_diff_all: Positions where target differs from ALL other sequences
            - sites_diff_any: Positions where target differs from ANY other sequence
            
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
        
        sites_diff_all = []  # Target differs from ALL others
        sites_diff_any = []  # Target differs from ANY others
        
        for pos in range(self.alignment_length):
            # Skip positions with insufficient flanking regions
            if pos < min_gap_left or pos >= self.alignment_length - min_gap_right:
                continue
            
            target_char = target_seq.sequence[pos]
            
            # Skip gaps in target
            if target_char == "-":
                continue
            
            # Check if flanking regions are gap-free
            if not self._is_gap_free_region(target_seq, pos, min_gap_left, min_gap_right):
                continue
            
            # Compare with other sequences
            diff_from_all = True
            diff_from_any = False
            
            for other_seq in other_seqs:
                other_char = other_seq.sequence[pos]
                
                # Skip gaps in other sequences
                if other_char == "-":
                    continue
                
                # Check if flanking regions are gap-free
                if not self._is_gap_free_region(other_seq, pos, min_gap_left, min_gap_right):
                    continue
                
                if target_char != other_char:
                    diff_from_any = True
                else:
                    diff_from_all = False
            
            if diff_from_any:
                sites_diff_any.append(pos)
            
            if diff_from_all:
                sites_diff_all.append(pos)
        
        return sites_diff_all, sites_diff_any
    
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