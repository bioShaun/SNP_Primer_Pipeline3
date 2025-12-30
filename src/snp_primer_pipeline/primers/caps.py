#!/usr/bin/env python3
"""
CAPS/dCAPS primer design module for SNP Primer Pipeline.

This module handles CAPS and dCAPS primer design using restriction enzymes.
"""

import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field

from ..models import Primer, PrimerPair, RestrictionEnzyme
from ..core.primer3_parser import Primer3Input, Primer3Runner, Primer3OutputParser
from ..core.alignment import MultipleSequenceAlignment
from ..exceptions import PrimerDesignError


class CAPSDesigner:
    """CAPS/dCAPS primer designer."""
    
    # IUPAC code mapping
    IUPAC_MAP = {
        "R": "AG", "Y": "TC", "S": "GC", 
        "W": "AT", "K": "TG", "M": "AC"
    }
    
    # IUPAC to regex pattern mapping
    IUPAC_PATTERNS = {
        "A": "A", "T": "T", "G": "G", "C": "C",
        "B": "[CGT]", "D": "[AGT]", "H": "[ACT]",
        "K": "[GT]", "M": "[AC]", "N": "[ACGT]",
        "R": "[AG]", "S": "[CG]", "V": "[ACG]",
        "W": "[AT]", "Y": "[CT]"
    }
    
    def __init__(
        self,
        primer3_path: Optional[Path] = None,
        config_path: Optional[Path] = None,
        enzyme_file: Optional[Path] = None,
        max_tm: float = 63.0,
        max_size: int = 25,
        pick_anyway: bool = False,
    ):
        """
        Initialize CAPS designer.
        
        Args:
            primer3_path: Path to primer3_core executable
            config_path: Path to Primer3 configuration directory
            enzyme_file: Path to restriction enzyme database file
            max_tm: Maximum primer Tm
            max_size: Maximum primer size
            pick_anyway: Pick primers anyway even if constraints violated
        """
        self.primer3_runner = Primer3Runner(primer3_path, config_path)
        self.primer3_parser = Primer3OutputParser()
        self.max_tm = max_tm
        self.max_size = max_size
        self.pick_anyway = pick_anyway
        self.enzyme_file = enzyme_file
        self.enzymes: Dict[str, RestrictionEnzyme] = {}
        
        if enzyme_file:
            self.load_enzymes()
    
    def load_enzymes(self) -> None:
        """Load restriction enzyme database from file."""
        if not self.enzyme_file or not Path(self.enzyme_file).exists():
            raise PrimerDesignError(f"Enzyme file not found: {self.enzyme_file}")
        
        try:
            with open(self.enzyme_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) != 2:
                        continue
                    
                    enzyme_info, sequence = parts
                    
                    # Parse enzyme name and price
                    # Format: "EnzymeName,price" or "EnzymeName1,EnzymeName2,price"
                    name_parts = enzyme_info.split(',')
                    price = int(name_parts[-1])
                    enzyme_name = name_parts[0]  # Use first name as primary
                    
                    enzyme = RestrictionEnzyme(
                        name=enzyme_name,
                        sequence=sequence.upper(),
                        price=price
                    )
                    
                    self.enzymes[enzyme_name] = enzyme
                    
        except (IOError, ValueError) as e:
            raise PrimerDesignError(f"Failed to load enzyme file: {e}") from e
    
    def find_usable_enzymes(
        self,
        template_sequence: str,
        snp_position: int,
        snp_alleles: Tuple[str, str],
        max_price: int = 200,
    ) -> Tuple[List[RestrictionEnzyme], List[RestrictionEnzyme]]:
        """
        Find enzymes that can be used for CAPS or dCAPS.
        
        Args:
            template_sequence: Template DNA sequence
            snp_position: SNP position in template (0-based)
            snp_alleles: Tuple of (allele_a, allele_b)
            max_price: Maximum enzyme price
            
        Returns:
            Tuple of (caps_enzymes, dcaps_enzymes)
        """
        allele_a, allele_b = snp_alleles
        
        # Create sequences with each allele
        wild_seq = (
            template_sequence[:snp_position] + 
            allele_a + 
            template_sequence[snp_position + 1:]
        ).lower()
        
        mut_seq = (
            template_sequence[:snp_position] + 
            allele_b + 
            template_sequence[snp_position + 1:]
        ).lower()
        
        caps_enzymes = []
        dcaps_enzymes = []
        
        for enzyme_name, enzyme in self.enzymes.items():
            if enzyme.price > max_price:
                continue
            
            # Test enzyme for CAPS/dCAPS capability
            tested_enzyme = self._test_enzyme(enzyme, wild_seq, mut_seq, snp_position)
            
            if tested_enzyme.caps:
                caps_enzymes.append(tested_enzyme)
            elif tested_enzyme.dcaps:
                dcaps_enzymes.append(tested_enzyme)
        
        return caps_enzymes, dcaps_enzymes
    
    def _test_enzyme(
        self,
        enzyme: RestrictionEnzyme,
        wild_seq: str,
        mut_seq: str,
        snp_position: int
    ) -> RestrictionEnzyme:
        """Test if enzyme can be used for CAPS or dCAPS."""
        # Create a copy to avoid modifying original
        test_enzyme = RestrictionEnzyme(
            name=enzyme.name,
            sequence=enzyme.sequence,
            price=enzyme.price
        )
        
        enzyme_seq = enzyme.sequence.lower()
        enzyme_seq_rc = self._reverse_complement(enzyme_seq)
        
        # Find all cutting positions in both sequences
        wild_positions = self._find_enzyme_sites(enzyme_seq, wild_seq)
        mut_positions = self._find_enzyme_sites(enzyme_seq, mut_seq)
        
        # Also check reverse complement
        wild_positions.extend(self._find_enzyme_sites(enzyme_seq_rc, wild_seq))
        mut_positions.extend(self._find_enzyme_sites(enzyme_seq_rc, mut_seq))
        
        test_enzyme.all_positions = list(set(wild_positions))
        
        # Check for CAPS (different number of cuts)
        if len(wild_positions) != len(mut_positions):
            test_enzyme.caps = True
            test_enzyme.template_seq = wild_seq
            return test_enzyme
        
        # Check for dCAPS (one base change can create/remove cut site)
        test_enzyme = self._check_dcaps(test_enzyme, wild_seq, mut_seq, snp_position)
        
        return test_enzyme
    
    def _find_enzyme_sites(self, enzyme_seq: str, sequence: str) -> List[int]:
        """Find all positions where enzyme cuts in sequence."""
        pattern = self._seq_to_pattern(enzyme_seq)
        return [m.start() for m in re.finditer(pattern, sequence)]
    
    def _seq_to_pattern(self, seq: str) -> str:
        """Convert enzyme sequence with IUPAC codes to regex pattern."""
        pattern = ""
        for base in seq.upper():
            pattern += self.IUPAC_PATTERNS.get(base, base)
        return pattern.lower()
    
    def _check_dcaps(
        self,
        enzyme: RestrictionEnzyme,
        wild_seq: str,
        mut_seq: str,
        snp_position: int
    ) -> RestrictionEnzyme:
        """Check if enzyme can be used for dCAPS."""
        enzyme_seq = enzyme.sequence.lower()
        
        # Try changing each position in enzyme sequence
        for i in range(len(enzyme_seq)):
            # Create pattern with one variable position
            pattern = self._seq_to_pattern(
                enzyme_seq[:i] + "N" + enzyme_seq[i+1:]
            )
            
            # Find matches in wild sequence
            for match in re.finditer(pattern, wild_seq):
                match_start = match.start()
                match_end = match.end()
                
                # Check if SNP is within match and change would affect cutting
                if (snp_position >= match_start and snp_position < match_end):
                    change_pos = match_start + i
                    
                    # Ensure change is at least 2 bases from SNP
                    if abs(snp_position - change_pos) > 1:
                        # Check if mutation prevents cutting
                        mut_region = mut_seq[match_start:match_end]
                        if not re.search(pattern, mut_region):
                            enzyme.dcaps = True
                            enzyme.template_seq = (
                                wild_seq[:change_pos] + 
                                enzyme_seq[i].upper() + 
                                wild_seq[change_pos + 1:]
                            )
                            enzyme.change_pos = change_pos + 1  # 1-based
                            
                            # Set primer end positions
                            if change_pos < snp_position:
                                enzyme.primer_end_positions = list(
                                    range(change_pos + 1, snp_position)
                                )
                            else:
                                enzyme.primer_end_positions = list(
                                    range(snp_position + 1, change_pos)
                                )
                            
                            return enzyme
        
        # Also try reverse complement
        enzyme_seq_rc = self._reverse_complement(enzyme_seq)
        if enzyme_seq_rc != enzyme_seq:
            enzyme.sequence = enzyme_seq_rc.upper()
            return self._check_dcaps(enzyme, wild_seq, mut_seq, snp_position)
        
        return enzyme
    
    def design_caps_primers(
        self,
        template_sequence: str,
        snp_position: int,
        snp_alleles: Tuple[str, str],
        enzyme: RestrictionEnzyme,
        product_size_range: Tuple[int, int] = (300, 900),
        variant_sites: Optional[List[int]] = None,
        output_dir: Optional[Path] = None,
    ) -> List[PrimerPair]:
        """
        Design CAPS primers for a specific enzyme.
        
        Args:
            template_sequence: Template DNA sequence
            snp_position: SNP position in template (0-based)
            snp_alleles: Tuple of (allele_a, allele_b)
            enzyme: Restriction enzyme to use
            product_size_range: Min and max product sizes
            variant_sites: List of variant positions for primer specificity
            output_dir: Output directory for intermediate files
            
        Returns:
            List of designed PrimerPair objects
        """
        if variant_sites is None:
            variant_sites = []
        
        product_min, product_max = product_size_range
        
        # Generate Primer3 inputs
        primer3_inputs = []
        
        if enzyme.dcaps:
            # dCAPS primers - force one primer to end at specific positions
            for primer_end_pos in enzyme.primer_end_positions:
                for var_site in variant_sites:
                    if primer_end_pos > snp_position:
                        left_end = var_site + 1
                        right_end = primer_end_pos + 1
                    else:
                        left_end = primer_end_pos + 1
                        right_end = var_site + 1
                    
                    # Check product size
                    if (right_end - left_end < product_min - 35 or 
                        right_end - left_end > product_max - 35):
                        continue
                    
                    p3_input = self._create_primer3_input(
                        enzyme.template_seq,
                        left_end,
                        right_end,
                        (product_min, product_max),
                        f"dCAPS-{enzyme.name}-{var_site+1}-{primer_end_pos+1}"
                    )
                    primer3_inputs.append(p3_input)
        
        elif enzyme.caps:
            # CAPS primers - target SNP region
            for var_site in variant_sites:
                if var_site < snp_position:
                    left_end = var_site + 1
                    right_end = -1000000  # Let Primer3 choose
                else:
                    left_end = -1000000  # Let Primer3 choose
                    right_end = var_site + 1
                
                p3_input = self._create_primer3_input(
                    enzyme.template_seq,
                    left_end,
                    right_end,
                    (product_min, product_max),
                    f"CAPS-{enzyme.name}-{var_site+1}",
                    target_region=(snp_position - 20, 40)
                )
                primer3_inputs.append(p3_input)
        
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
        
        return primer_pairs
    
    def _create_primer3_input(
        self,
        template: str,
        left_end: int,
        right_end: int,
        product_size_range: Tuple[int, int],
        sequence_id: str,
        target_region: Optional[Tuple[int, int]] = None
    ) -> str:
        """Create Primer3 input string."""
        p3_input = Primer3Input()
        p3_input.set_template(template)
        p3_input.set_product_size_range([product_size_range])
        
        if left_end > 0:
            p3_input.set_force_left_end(left_end)
        if right_end > 0:
            p3_input.set_force_right_end(right_end)
        
        if target_region:
            start, length = target_region
            p3_input.set_target(start, length)
        
        p3_input.set_setting("PRIMER_MAX_SIZE", self.max_size)
        p3_input.set_setting("PRIMER_MAX_TM", self.max_tm)
        p3_input.set_setting("PRIMER_PICK_ANYWAY", 1 if self.pick_anyway else 0)
        p3_input.set_setting("PRIMER_NUM_RETURN", 5)
        
        return p3_input.generate(sequence_id)
    
    def _reverse_complement(self, seq: str) -> str:
        """Get reverse complement of DNA sequence."""
        complement = {
            "a": "t", "t": "a", "g": "c", "c": "g",
            "A": "T", "T": "A", "G": "C", "C": "G",
            "r": "y", "y": "r", "s": "s", "w": "w",
            "k": "m", "m": "k", "b": "v", "v": "b",
            "d": "h", "h": "d", "n": "n",
            "R": "Y", "Y": "R", "S": "S", "W": "W",
            "K": "M", "M": "K", "B": "V", "V": "B",
            "D": "H", "H": "D", "N": "N"
        }
        return "".join(complement.get(base, base) for base in reversed(seq))
    
    def format_output(
        self,
        primer_pairs: List[PrimerPair],
        snp_name: str,
        output_file: Path,
        caps_enzymes: List[RestrictionEnzyme],
        dcaps_enzymes: List[RestrictionEnzyme],
        variant_sites: Optional[List[int]] = None
    ) -> None:
        """
        Format and write CAPS primer results to file.
        
        Args:
            primer_pairs: List of designed primer pairs
            snp_name: SNP identifier
            output_file: Output file path
            caps_enzymes: List of CAPS enzymes found
            dcaps_enzymes: List of dCAPS enzymes found
            variant_sites: List of variant sites for annotation
        """
        if variant_sites is None:
            variant_sites = []
        
        try:
            with open(output_file, 'w') as f:
                # Write header
                header = [
                    "index", "product_size", "type", "start", "end", 
                    "diff_number", "3'differall", "length", "Tm", 
                    "GC_content", "any", "3'", "end_stability", "hairpin", 
                    "primer_seq", "reverse_complement", "penalty", 
                    "compl_any", "compl_end"
                ]
                f.write("\t".join(header) + "\n")
                
                # Write primer pairs
                for i, pair in enumerate(primer_pairs):
                    index = f"{snp_name}-{i}"
                    
                    # Format left primer
                    left_line = self._format_primer_line(
                        index, pair.product_size, pair.left, 
                        pair.penalty, pair.compl_any, pair.compl_end
                    )
                    f.write(left_line + "\n")
                    
                    # Format right primer
                    right_line = self._format_primer_line(
                        index, pair.product_size, pair.right,
                        pair.penalty, pair.compl_any, pair.compl_end
                    )
                    f.write(right_line + "\n")
                
                # Write variant sites information
                if variant_sites:
                    f.write(f"\n\nSites that can differ all for {snp_name}\n")
                    f.write(", ".join(str(x + 1) for x in variant_sites))  # Convert to 1-based
                    f.write("\n")
                
                # Write CAPS enzyme information
                f.write(f"\nCAPS cut information for SNP {snp_name}\n")
                f.write("Enzyme\tEnzyme_seq\tCut_positions\n")
                for enzyme in caps_enzymes:
                    positions = ", ".join(str(x + 1) for x in enzyme.all_positions)
                    f.write(f"{enzyme.name}\t{enzyme.sequence}\t{positions}\n")
                
                # Write dCAPS enzyme information
                f.write(f"\ndCAPS cut information for SNP {snp_name}\n")
                f.write("Enzyme\tEnzyme_seq\tChange_pos\tOther_cut_pos\tPotential_primer\n")
                for enzyme in dcaps_enzymes:
                    positions = ", ".join(str(x + 1) for x in enzyme.all_positions)
                    potential_primer = getattr(enzyme, 'potential_primer', '')
                    f.write(f"{enzyme.name}\t{enzyme.sequence}\t{enzyme.change_pos}\t{positions}\t{potential_primer}\n")
                
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
        compl_end: float
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
            f"{compl_end:.2f}"
        ])