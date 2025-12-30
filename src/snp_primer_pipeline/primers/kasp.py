#!/usr/bin/env python3
"""
KASP primer design module for SNP Primer Pipeline.

This module handles KASP primer design using Primer3 and multiple sequence alignment.
"""

import copy
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import re

from ..models import Primer, PrimerPair
from ..core.primer3_parser import Primer3Input, Primer3Runner, Primer3OutputParser
from ..core.alignment import MultipleSequenceAlignment
from ..exceptions import PrimerDesignError


class KASPDesigner:
    """KASP primer designer."""
    
    # IUPAC code mapping
    IUPAC_MAP = {
        "R": "AG", "Y": "TC", "S": "GC", 
        "W": "AT", "K": "TG", "M": "AC"
    }
    
    def __init__(
        self,
        primer3_path: Optional[Path] = None,
        config_path: Optional[Path] = None,
        max_tm: float = 63.0,
        max_size: int = 25,
        pick_anyway: bool = False,
    ):
        """
        Initialize KASP designer.
        
        Args:
            primer3_path: Path to primer3_core executable
            config_path: Path to Primer3 configuration directory
            max_tm: Maximum primer Tm
            max_size: Maximum primer size
            pick_anyway: Pick primers anyway even if constraints violated
        """
        # Auto-detect primer3 path if not provided
        if primer3_path is None:
            from ..config import SoftwarePaths
            software_paths = SoftwarePaths.auto_detect()
            primer3_path = software_paths.primer3_path
            
        # Auto-detect config path if not provided
        if config_path is None:
            config_path = Path(__file__).parent.parent.parent / "bin" / "primer3_config"
            
        # Don't pass config_path as settings_file - we handle thermodynamic parameters in Primer3Input
        self.primer3_runner = Primer3Runner(primer3_path, settings_file=None)
        self.primer3_parser = Primer3OutputParser()
        self.max_tm = max_tm
        self.max_size = max_size
        self.pick_anyway = pick_anyway
    
    def design_primers(
        self,
        template_sequence: str,
        snp_position: int,
        snp_alleles: Tuple[str, str],
        product_size_range: Tuple[int, int] = (50, 250),
        variant_sites: Optional[List[int]] = None,
        alignment: Optional[MultipleSequenceAlignment] = None,
        target_name: Optional[str] = None,
        output_dir: Optional[Path] = None,
    ) -> List[PrimerPair]:
        """
        Design KASP primers for a SNP.
        
        Args:
            template_sequence: Template DNA sequence
            snp_position: SNP position in template (0-based)
            snp_alleles: Tuple of (allele_a, allele_b)
            product_size_range: Min and max product sizes
            variant_sites: List of variant positions for primer specificity
            alignment: Multiple sequence alignment for specificity checking
            target_name: Name of target sequence in alignment
            output_dir: Output directory for intermediate files
            
        Returns:
            List of designed PrimerPair objects
        """
        if variant_sites is None:
            variant_sites = []
        
        # Determine which allele to use as template (prefer A/T for lower Tm)
        allele_a, allele_b = snp_alleles
        template_allele = allele_a if allele_a in "AT" else allele_b
        
        # Modify template sequence with preferred allele
        modified_template = (
            template_sequence[:snp_position] + 
            template_allele + 
            template_sequence[snp_position + 1:]
        )
        
        # Generate Primer3 inputs for each variant site
        primer3_inputs = self._generate_primer3_inputs(
            modified_template,
            snp_position,
            variant_sites,
            product_size_range
        )
        
        if not primer3_inputs:
            # No variant sites found, design primers with forced ends at SNP
            primer3_inputs = self._generate_snp_forced_inputs(
                modified_template,
                snp_position,
                product_size_range
            )
        
        # Run Primer3 for all inputs
        primer_pairs = []
        for input_data in primer3_inputs:
            try:
                output_string = self.primer3_runner.run_string(input_data)
                results = self.primer3_parser.parse_string(output_string, num_pairs=1)
                
                for seq_id, pairs in results.items():
                    primer_pairs.extend(pairs)
            except PrimerDesignError:
                # Skip failed designs
                continue
        
        # Filter and score primer pairs
        filtered_pairs = self._filter_and_score_primers(
            primer_pairs,
            snp_position,
            variant_sites,
            alignment,
            target_name
        )
        
        # Convert to KASP format (add allele-specific primers)
        kasp_primers = self._convert_to_kasp_format(
            filtered_pairs,
            snp_position,
            snp_alleles
        )
        
        return kasp_primers
    
    def _generate_primer3_inputs(
        self,
        template: str,
        snp_position: int,
        variant_sites: List[int],
        product_size_range: Tuple[int, int]
    ) -> List[str]:
        """Generate Primer3 input strings for variant sites."""
        inputs = []
        product_min, product_max = product_size_range
        
        for var_site in variant_sites:
            if var_site == snp_position:
                continue
            
            # Determine primer positions
            if var_site < snp_position:
                left_end = var_site
                right_end = snp_position
            else:
                left_end = snp_position
                right_end = var_site
            
            # Skip if product would be too large
            if right_end - left_end > product_max - 35:  # Account for primer lengths
                continue
            
            # Create Primer3 input
            p3_input = Primer3Input()
            p3_input.set_template(template)
            p3_input.set_product_size_range([(product_min, product_max)])
            p3_input.set_force_left_end(left_end + 1)  # Convert to 1-based
            p3_input.set_force_right_end(right_end + 1)
            p3_input.set_setting("PRIMER_MAX_SIZE", self.max_size)
            p3_input.set_setting("PRIMER_MAX_TM", self.max_tm)
            p3_input.set_setting("PRIMER_PICK_ANYWAY", 1 if self.pick_anyway else 0)
            p3_input.set_setting("PRIMER_NUM_RETURN", 1)
            
            sequence_id = f"var_{var_site + 1}"  # 1-based for output
            inputs.append(p3_input.generate(sequence_id))
        
        return inputs
    
    def _generate_snp_forced_inputs(
        self,
        template: str,
        snp_position: int,
        product_size_range: Tuple[int, int]
    ) -> List[str]:
        """Generate Primer3 inputs with forced ends at SNP position."""
        inputs = []
        
        # Left primer ending at SNP
        p3_input_left = Primer3Input()
        p3_input_left.set_template(template)
        p3_input_left.set_product_size_range([product_size_range])
        p3_input_left.set_force_left_end(snp_position + 1)  # 1-based
        p3_input_left.set_setting("PRIMER_MAX_SIZE", self.max_size)
        p3_input_left.set_setting("PRIMER_MAX_TM", self.max_tm)
        p3_input_left.set_setting("PRIMER_PICK_ANYWAY", 1 if self.pick_anyway else 0)
        p3_input_left.set_setting("PRIMER_NUM_RETURN", 5)
        
        inputs.append(p3_input_left.generate("left"))
        
        # Right primer ending at SNP
        p3_input_right = Primer3Input()
        p3_input_right.set_template(template)
        p3_input_right.set_product_size_range([product_size_range])
        p3_input_right.set_force_right_end(snp_position + 1)  # 1-based
        p3_input_right.set_setting("PRIMER_MAX_SIZE", self.max_size)
        p3_input_right.set_setting("PRIMER_MAX_TM", self.max_tm)
        p3_input_right.set_setting("PRIMER_PICK_ANYWAY", 1 if self.pick_anyway else 0)
        p3_input_right.set_setting("PRIMER_NUM_RETURN", 5)
        
        inputs.append(p3_input_right.generate("right"))
        
        return inputs
    
    def _filter_and_score_primers(
        self,
        primer_pairs: List[PrimerPair],
        snp_position: int,
        variant_sites: List[int],
        alignment: Optional[MultipleSequenceAlignment] = None,
        target_name: Optional[str] = None
    ) -> List[PrimerPair]:
        """Filter and score primer pairs."""
        filtered_pairs = []
        
        for pair in primer_pairs:
            if pair.product_size == 0:
                continue
            
            # Check if primers contain variant sites
            self._annotate_primer_variants(pair.left, variant_sites)
            self._annotate_primer_variants(pair.right, variant_sites)
            
            # Calculate score
            pair.score = self.calculate_score(pair, snp_position, variant_sites)
            
            # Check specificity if alignment provided
            if alignment and target_name:
                if self._check_primer_specificity(pair, alignment, target_name):
                    pair.left.diff_three_all = True
                    pair.right.diff_three_all = True
            
            filtered_pairs.append(pair)
        
        # Sort by score (higher is better)
        filtered_pairs.sort(key=lambda x: x.score, reverse=True)
        
        return filtered_pairs
    
    def _annotate_primer_variants(self, primer: Primer, variant_sites: List[int]) -> None:
        """Annotate primer with variant site information."""
        if primer.direction == "LEFT":
            primer_range = list(range(primer.start, primer.end + 1))
        else:
            primer_range = list(range(primer.end, primer.start + 1))
        
        var_sites_in_primer = set(variant_sites).intersection(primer_range)
        primer.diff_num = len(var_sites_in_primer)
        
        # Check if 3' end can differentiate all homeologs
        if primer.direction == "LEFT":
            three_prime_pos = primer.end
        else:
            three_prime_pos = primer.end
        
        primer.diff_three_all = three_prime_pos in variant_sites
    
    def _check_primer_specificity(
        self,
        pair: PrimerPair,
        alignment: MultipleSequenceAlignment,
        target_name: str
    ) -> bool:
        """Check if primers can differentiate target from homeologs."""
        # This is a simplified check - in practice, you'd want to implement
        # the full homeolog sequence comparison logic from the original code
        return True
    
    def calculate_score(
        self,
        primer_pair: PrimerPair,
        snp_position: int,
        variant_sites: List[int]
    ) -> float:
        """
        Calculate primer pair score for ranking.
        
        Args:
            primer_pair: PrimerPair to score
            snp_position: SNP position
            variant_sites: List of variant sites
            
        Returns:
            Score (higher is better)
        """
        score = 0.0
        
        # Product size score (smaller is better)
        score += 150.0 / primer_pair.product_size
        
        # 3' differentiation capability
        if primer_pair.left.diff_three_all or primer_pair.right.diff_three_all:
            score += 5.0
        
        # Tm difference penalty (smaller difference is better)
        tm_diff = abs(primer_pair.left.tm - primer_pair.right.tm)
        score -= tm_diff / 10.0
        
        # Variant sites in primers (more is better)
        total_variants = primer_pair.left.diff_num + primer_pair.right.diff_num
        score += total_variants / 10.0
        
        return score
    
    def _convert_to_kasp_format(
        self,
        primer_pairs: List[PrimerPair],
        snp_position: int,
        snp_alleles: Tuple[str, str]
    ) -> List[PrimerPair]:
        """Convert primer pairs to KASP format with allele-specific primers.
        
        V2 algorithm:
        - If left primer ends at SNP position: left becomes allele-specific, right is common
        - If right primer ends at SNP position: right becomes allele-specific, left is common
        - The allele-specific primer's last base is REPLACED with the SNP allele
        """
        kasp_primers = []
        allele_a, allele_b = snp_alleles
        
        for pair in primer_pairs:
            # Determine which primer ends at the SNP position
            # V2 uses 1-based positions, so snp_position + 1 in 1-based = snp_position in 0-based
            # For LEFT primers: end is the 3' position (0-based)
            # For RIGHT primers: end is the 5' position (0-based), start is 3' position
            
            left_ends_at_snp = (pair.left.end == snp_position)
            right_ends_at_snp = (pair.right.end == snp_position)
            
            if left_ends_at_snp:
                # Left primer is allele-specific, right is common
                primer_a = copy.deepcopy(pair.left)
                primer_b = copy.deepcopy(pair.left)
                common_primer = copy.deepcopy(pair.right)
                
                # Replace last base with allele (V2 style)
                primer_a.sequence = primer_a.sequence[:-1] + allele_a
                primer_b.sequence = primer_b.sequence[:-1] + allele_b
                
            elif right_ends_at_snp:
                # Right primer is allele-specific, left is common
                primer_a = copy.deepcopy(pair.right)
                primer_b = copy.deepcopy(pair.right)
                common_primer = copy.deepcopy(pair.left)
                
                # For right primers, use reverse complement of allele
                primer_a.sequence = primer_a.sequence[:-1] + self._reverse_complement(allele_a)
                primer_b.sequence = primer_b.sequence[:-1] + self._reverse_complement(allele_b)
                
            else:
                # Neither primer ends at SNP - use left as allele-specific by default
                primer_a = copy.deepcopy(pair.left)
                primer_b = copy.deepcopy(pair.left)
                common_primer = copy.deepcopy(pair.right)
                
                # Append allele to end
                primer_a.sequence = primer_a.sequence + allele_a
                primer_b.sequence = primer_b.sequence + allele_b
                primer_a.length = len(primer_a.sequence)
                primer_b.length = len(primer_b.sequence)
                primer_a.end = primer_a.start + primer_a.length - 1
                primer_b.end = primer_b.start + primer_b.length - 1
            
            # Update primer names
            primer_a.name = f"Allele_{allele_a}"
            primer_b.name = f"Allele_{allele_b}"
            common_primer.name = "Common"
            
            # Create KASP primer pairs
            kasp_pair_a = PrimerPair(
                left=primer_a,
                right=common_primer,
                product_size=pair.product_size,
                penalty=pair.penalty,
                compl_any=pair.compl_any,
                compl_end=pair.compl_end,
                score=pair.score
            )
            kasp_pair_b = PrimerPair(
                left=primer_b,
                right=common_primer,
                product_size=pair.product_size,
                penalty=pair.penalty,
                compl_any=pair.compl_any,
                compl_end=pair.compl_end,
                score=pair.score
            )
            
            kasp_primers.extend([kasp_pair_a, kasp_pair_b])
        
        return kasp_primers
    
    def _reverse_complement(self, seq: str) -> str:
        """Get reverse complement of DNA sequence."""
        complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
        return "".join(complement.get(base, base) for base in reversed(seq.upper()))
    
    def format_output(
        self,
        primer_pairs: List[PrimerPair],
        snp_name: str,
        output_file: Path,
        variant_sites: Optional[List[int]] = None
    ) -> None:
        """
        Format and write KASP primer results to file.
        
        Args:
            primer_pairs: List of designed primer pairs
            snp_name: SNP identifier
            output_file: Output file path
            variant_sites: List of variant sites for annotation
        """
        if variant_sites is None:
            variant_sites = []
        
        try:
            with open(output_file, 'w') as f:
                # Write header
                header = [
                    "index", "product_size", "type", "start", "end", 
                    "variation_number", "3'diffall", "length", "Tm", 
                    "GC_content", "any", "3'", "end_stability", "hairpin", 
                    "primer_seq", "reverse_complement", "penalty", 
                    "compl_any", "compl_end", "score"
                ]
                f.write("\t".join(header) + "\n")
                
                # Write primer pairs
                for i, pair in enumerate(primer_pairs):
                    index = f"{snp_name}-{i}"
                    
                    # Format left primer
                    left_line = self._format_primer_line(
                        index, pair.product_size, pair.left, 
                        pair.penalty, pair.compl_any, pair.compl_end, pair.score
                    )
                    f.write(left_line + "\n")
                    
                    # Format right primer
                    right_line = self._format_primer_line(
                        index, pair.product_size, pair.right,
                        pair.penalty, pair.compl_any, pair.compl_end, pair.score
                    )
                    f.write(right_line + "\n")
                
                # Write variant sites information
                if variant_sites:
                    f.write(f"\n\nSites that can differ all for {snp_name}\n")
                    f.write(", ".join(str(x + 1) for x in variant_sites))  # Convert to 1-based
                    f.write("\n\n")
                    
        except IOError as e:
            raise PrimerDesignError(f"Failed to write output file: {e}") from e
    
    def _format_primer_line(
        self,
        index: str,
        product_size: int,
        primer: Primer,
        penalty: float,
        compl_any: float,
        compl_end: float,
        score: float
    ) -> str:
        """Format a single primer line for output."""
        return "\t".join([
            index,
            str(product_size),
            primer.direction,
            str(primer.start),
            str(primer.end),
            str(primer.diff_num),
            "YES" if primer.diff_three_all else "NO",
            str(primer.length),
            f"{primer.tm:.2f}",
            f"{primer.gc_percent:.2f}",
            f"{primer.self_any:.2f}",
            f"{primer.self_end:.2f}",
            f"{primer.end_stability:.2f}",
            f"{primer.hairpin:.2f}",
            primer.sequence,
            self._reverse_complement(primer.sequence),
            f"{penalty:.2f}",
            f"{compl_any:.2f}",
            f"{compl_end:.2f}",
            f"{score:.2f}"
        ])