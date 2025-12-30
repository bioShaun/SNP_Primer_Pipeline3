#!/usr/bin/env python3
"""
KASP primer design module for SNP Primer Pipeline.

This module handles KASP primer design using Primer3 and multiple sequence alignment.
Refactored to match V2 behavior exactly.
"""

import copy
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import re

from ..models import Primer, PrimerPair
from ..core.primer3_parser import Primer3Input, Primer3Runner, Primer3OutputParser
from ..core.alignment import MultipleSequenceAlignment
from ..exceptions import PrimerDesignError


class KASPDesigner:
    """KASP primer designer matching V2 behavior."""
    
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
        self.config_path = config_path
        self.max_tm = max_tm
        self.max_size = max_size
        self.pick_anyway = pick_anyway
        
        # Store diffarray for V2 compatibility (position -> list of diffs from each homeolog)
        self.diffarray: Dict[int, List[int]] = {}
    
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
    ) -> List[Dict[str, Any]]:
        """
        Design KASP primers for a SNP (V2-compatible output).
        
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
            List of KASP primer dictionaries (V2-compatible format)
        """
        if variant_sites is None:
            variant_sites = []
        
        # Store variant sites for later use
        self._variant_sites = variant_sites
        self._snp_position = snp_position
        self._snp_alleles = snp_alleles
        self._template_sequence = template_sequence
        
        # Build diffarray from alignment if available (V2 compatibility)
        if alignment and target_name:
            self._build_diffarray(alignment, target_name, snp_position)
        
        # Determine which allele to use as template (prefer A/T for lower Tm)
        allele_a, allele_b = snp_alleles
        template_allele = allele_a if allele_a in "AT" else allele_b
        
        # Modify template sequence with preferred allele
        modified_template = (
            template_sequence[:snp_position] + 
            template_allele + 
            template_sequence[snp_position + 1:]
        )
        
        # If no homeologs found, use SNP-forced primers only
        if not variant_sites:
            return self._design_primers_no_homeologs(
                modified_template, snp_position, snp_alleles, product_size_range, output_dir
            )
        
        # Generate Primer3 inputs for each variant site (V2 style)
        primer3_inputs = self._generate_primer3_inputs_v2(
            modified_template,
            snp_position,
            variant_sites,
            product_size_range
        )
        
        if not primer3_inputs:
            return []
        
        # Run Primer3 and collect results
        raw_primer_pairs = {}
        for var_site, input_data in primer3_inputs:
            try:
                output_string = self.primer3_runner.run_string(input_data)
                results = self.primer3_parser.parse_string(output_string, num_pairs=1)
                
                for seq_id, pairs in results.items():
                    for pair in pairs:
                        if pair.product_size != 0:
                            # Use var_site + 1 for 1-based naming (V2 style)
                            key = f"{var_site + 1}-0"
                            raw_primer_pairs[key] = (pair, var_site)
            except PrimerDesignError:
                continue
        
        # Filter primer pairs (V2 style - check if common primer can differ all)
        final_primers = self._filter_primers_v2(
            raw_primer_pairs, snp_position, variant_sites, len(template_sequence)
        )
        
        # Convert to V2 output format
        kasp_output = self._convert_to_v2_format(
            final_primers, snp_position, snp_alleles, variant_sites
        )
        
        return kasp_output
    
    def _build_diffarray(
        self,
        alignment: MultipleSequenceAlignment,
        target_name: str,
        snp_position: int
    ) -> None:
        """
        Build diffarray for V2-compatible filtering.
        diffarray[pos] = list of 0/1 indicating if each homeolog differs at this position.
        """
        self.diffarray = {}
        
        target_seq = alignment.get_sequence_by_name(target_name)
        if not target_seq:
            return
        
        other_seqs = [s for s in alignment.sequences if s.name != target_name]
        if not other_seqs:
            return
        
        # Get clean template sequence
        template = target_seq.clean_sequence
        
        for pos in range(len(template)):
            align_pos = alignment.template_to_alignment(pos)
            if align_pos is None:
                continue
            
            target_base = target_seq.sequence[align_pos]
            if target_base == '-':
                continue
            
            diff_list = []
            for other_seq in other_seqs:
                other_base = other_seq.sequence[align_pos]
                if other_base == '-':
                    diff_list.append(0)  # Gap treated as no difference for V2
                elif target_base != other_base:
                    diff_list.append(1)
                else:
                    diff_list.append(0)
            
            self.diffarray[pos] = diff_list
    
    def _generate_primer3_inputs_v2(
        self,
        template: str,
        snp_position: int,
        variant_sites: List[int],
        product_size_range: Tuple[int, int]
    ) -> List[Tuple[int, str]]:
        """
        Generate Primer3 input strings for variant sites (V2 style).
        Returns list of (var_site, input_string) tuples.
        """
        inputs = []
        product_min, product_max = product_size_range
        
        for var_site in variant_sites:
            if var_site == snp_position:
                continue
            
            # Determine primer positions (V2 style)
            if var_site < snp_position:
                left_end = var_site
                right_end = snp_position
            else:
                left_end = snp_position
                right_end = var_site
            
            # Skip if product would be too large (V2: product_max - 35)
            if right_end - left_end > product_max - 35:
                continue
            
            # Create Primer3 input (V2 style settings)
            p3_input = Primer3Input()
            p3_input.set_template(template)
            p3_input.set_product_size_range([(50, 100), (100, 150), (150, 250)])
            p3_input.set_force_left_end(left_end + 1)  # Convert to 1-based
            p3_input.set_force_right_end(right_end + 1)
            p3_input.set_setting("PRIMER_MAX_SIZE", self.max_size)
            p3_input.set_setting("PRIMER_MIN_TM", 57.0)
            p3_input.set_setting("PRIMER_OPT_TM", 60.0)
            p3_input.set_setting("PRIMER_MAX_TM", self.max_tm)
            p3_input.set_setting("PRIMER_PAIR_MAX_DIFF_TM", 6.0)
            p3_input.set_setting("PRIMER_FIRST_BASE_INDEX", 1)
            p3_input.set_setting("PRIMER_LIBERAL_BASE", 1)
            p3_input.set_setting("PRIMER_NUM_RETURN", 5)
            p3_input.set_setting("PRIMER_EXPLAIN_FLAG", 1)
            p3_input.set_setting("PRIMER_PICK_ANYWAY", 1 if self.pick_anyway else 0)
            
            sequence_id = f"var_{var_site + 1}"  # 1-based for output
            inputs.append((var_site, p3_input.generate(sequence_id)))
        
        return inputs
    
    def _design_primers_no_homeologs(
        self,
        template: str,
        snp_position: int,
        snp_alleles: Tuple[str, str],
        product_size_range: Tuple[int, int],
        output_dir: Optional[Path] = None
    ) -> List[Dict[str, Any]]:
        """
        Design primers when no homeologs found (V2 style).
        Force primer ends at SNP position.
        """
        results = []
        allele_a, allele_b = snp_alleles
        
        # Generate inputs for left and right primers ending at SNP
        for direction in ['left', 'right']:
            p3_input = Primer3Input()
            p3_input.set_template(template)
            p3_input.set_product_size_range([(50, 100), (100, 150), (150, 250)])
            
            if direction == 'left':
                p3_input.set_force_left_end(snp_position + 1)  # 1-based
            else:
                p3_input.set_force_right_end(snp_position + 1)  # 1-based
            
            p3_input.set_setting("PRIMER_MAX_SIZE", self.max_size)
            p3_input.set_setting("PRIMER_MIN_TM", 57.0)
            p3_input.set_setting("PRIMER_OPT_TM", 60.0)
            p3_input.set_setting("PRIMER_MAX_TM", self.max_tm)
            p3_input.set_setting("PRIMER_PAIR_MAX_DIFF_TM", 6.0)
            p3_input.set_setting("PRIMER_FIRST_BASE_INDEX", 1)
            p3_input.set_setting("PRIMER_LIBERAL_BASE", 1)
            p3_input.set_setting("PRIMER_NUM_RETURN", 5)
            p3_input.set_setting("PRIMER_EXPLAIN_FLAG", 1)
            p3_input.set_setting("PRIMER_PICK_ANYWAY", 1 if self.pick_anyway else 0)
            
            try:
                output_string = self.primer3_runner.run_string(p3_input.generate(direction))
                parsed = self.primer3_parser.parse_string(output_string, num_pairs=5)
                
                for seq_id, pairs in parsed.items():
                    for idx, pair in enumerate(pairs):
                        if pair.product_size == 0:
                            continue
                        
                        # Format output (V2 style - no homeolog filtering)
                        result = self._format_primer_pair_v2(
                            pair, snp_position, snp_alleles, [], 
                            f"{direction}-{idx}"
                        )
                        results.extend(result)
            except PrimerDesignError:
                continue
        
        return results
    
    def _filter_primers_v2(
        self,
        raw_primer_pairs: Dict[str, Tuple[PrimerPair, int]],
        snp_position: int,
        variant_sites: List[int],
        template_length: int
    ) -> Dict[str, Tuple[PrimerPair, int]]:
        """
        Filter primer pairs using V2 logic:
        Only keep pairs where common primer can differ all homeologs in 3' 10bp region.
        
        In V3:
        - LEFT primer: start is 5' end, end is 3' end
        - RIGHT primer: start is 3' end, end is 5' end
        """
        final_primers = {}
        gap_left = 20  # V2 default
        
        for key, (pp, var_site) in raw_primer_pairs.items():
            # Check if var_site is in the variant list for 3'diffall
            dif3all = 1 if var_site in variant_sites else 0
            
            # Determine which is the common primer (V2 logic)
            if var_site < snp_position:
                # Left primer is common (ends at var_site, which is before SNP)
                pc = pp.left
                # For LEFT primer, 3' end is at pc.end
                # Range to check: 10 bases from 3' end
                three_prime_pos = pc.end
                rr = list(range(max(three_prime_pos - 9, gap_left), three_prime_pos + 1))
            else:
                # Right primer is common (ends at var_site, which is after SNP)
                pc = pp.right
                # For RIGHT primer in V3, 3' end is at pc.start (lower position)
                three_prime_pos = pc.start
                rr = list(range(three_prime_pos, min(three_prime_pos + 10, template_length - 20)))
            
            # Calculate score (V2 style)
            pp.score = dif3all * 5.0 + 150.0 / pp.product_size + pc.diff_num / 10.0 - abs(pp.left.tm - pp.right.tm) / 10.0
            
            # Check if common primer can differ all homeologs in its 3' region
            if self.diffarray:
                # Sum differences for each homeolog across the range
                try:
                    aa = []
                    for k in rr:
                        if k in self.diffarray:
                            aa.append(self.diffarray[k])
                    
                    if aa:
                        # Sum for each homeolog
                        sums = [sum(x) for x in zip(*aa)]
                        # Only keep if common primer can differ from all homeologs
                        if min(sums) > 0:
                            final_primers[key] = (pp, var_site)
                    else:
                        # No diffarray data, keep primer
                        final_primers[key] = (pp, var_site)
                except Exception:
                    # On any error, keep the primer
                    final_primers[key] = (pp, var_site)
            else:
                # No diffarray, keep all primers
                final_primers[key] = (pp, var_site)
        
        return final_primers
    
    def _format_primer_seq_v2(self, primer: Primer, variant_sites: List[int]) -> str:
        """
        Format primer sequence with variant sites highlighted (V2 style).
        Lowercase for normal bases, uppercase for variant positions.
        
        V2 convention:
        - LEFT primer: start < end, start is 5' end, end is 3' end
        - RIGHT primer (in V2): start > end, start is 5' end, end is 3' end
        - In V3: start < end for both, so for RIGHT: start is 3' end, end is 5' end
        """
        # Get primer range in template coordinates (0-based)
        if primer.direction == "LEFT":
            # LEFT: start is 5' end, end is 3' end
            range_start = primer.start
            range_end = primer.end
            seq = primer.sequence.lower()
        else:
            # RIGHT: in v3, start < end where start is 3' end, end is 5' end
            range_start = primer.start
            range_end = primer.end
            seq = self._reverse_complement(primer.sequence).lower()
        
        # Find variant sites within primer range
        primer_range = set(range(range_start, range_end + 1))
        var_sites_in_primer = set(variant_sites).intersection(primer_range)
        
        # Calculate relative positions and uppercase them
        for var_pos in var_sites_in_primer:
            relative_pos = var_pos - range_start
            if 0 <= relative_pos < len(seq):
                seq = seq[:relative_pos] + seq[relative_pos].upper() + seq[relative_pos + 1:]
        
        # Return in correct orientation for output
        if primer.direction == "LEFT":
            return seq
        else:
            return self._reverse_complement(seq)
    
    def _count_variants_in_primer(self, primer: Primer, variant_sites: List[int]) -> int:
        """Count number of variant sites within primer region."""
        if primer.direction == "LEFT":
            primer_range = set(range(primer.start, primer.end + 1))
        else:
            primer_range = set(range(primer.end, primer.start + 1))
        
        return len(primer_range.intersection(variant_sites))
    
    def _convert_to_v2_format(
        self,
        final_primers: Dict[str, Tuple[PrimerPair, int]],
        snp_position: int,
        snp_alleles: Tuple[str, str],
        variant_sites: List[int]
    ) -> List[Dict[str, Any]]:
        """
        Convert filtered primers to V2 output format.
        Returns list of dictionaries with V2-compatible fields.
        """
        results = []
        allele_a, allele_b = snp_alleles
        
        for key, (pp, var_site) in final_primers.items():
            result = self._format_primer_pair_v2(
                pp, snp_position, snp_alleles, variant_sites, key
            )
            results.extend(result)
        
        return results
    
    def _format_primer_pair_v2(
        self,
        pp: PrimerPair,
        snp_position: int,
        snp_alleles: Tuple[str, str],
        variant_sites: List[int],
        key: str
    ) -> List[Dict[str, Any]]:
        """
        Format a single primer pair in V2 output format.
        Returns list of 3 dictionaries: Allele-A, Allele-B, Common
        """
        allele_a, allele_b = snp_alleles
        results = []
        
        pl = pp.left
        pr = pp.right
        
        # Annotate variants in primers
        pl.diff_num = self._count_variants_in_primer(pl, variant_sites)
        pr.diff_num = self._count_variants_in_primer(pr, variant_sites)
        
        # Parse var_site from key (format: "var_site-index")
        # In V2, both primers get the same difthreeall based on whether varsite is in variation
        try:
            var_site_from_key = int(key.split('-')[0]) - 1  # Convert to 0-based
        except (ValueError, IndexError):
            var_site_from_key = -1
        
        # Check 3'diffall (V2 style: based on whether var_site is in variation list)
        dif3all = var_site_from_key in variant_sites
        pl.diff_three_all = dif3all
        pr.diff_three_all = dif3all
        
        # Format primer sequences with highlighted variants
        pl_seq_formatted = self._format_primer_seq_v2(pl, variant_sites)
        pr_seq_formatted = self._format_primer_seq_v2(pr, variant_sites)
        
        # Determine allele-specific and common primers (V2 logic)
        # If left primer ends at SNP: left is allele-specific
        # Otherwise: right is allele-specific
        if pl.end == snp_position:
            # Left is allele-specific, right is common
            pA_seq = pl_seq_formatted[:-1] + allele_a
            pB_seq = pl_seq_formatted[:-1] + allele_b
            pC_seq = pr_seq_formatted
            pA_primer = pl
            pB_primer = pl
            pC_primer = pr
            pA_dir = "LEFT"
            pB_dir = "LEFT"
            pC_dir = "RIGHT"
        else:
            # Right is allele-specific, left is common
            pA_seq = pr_seq_formatted[:-1] + self._reverse_complement(allele_a)
            pB_seq = pr_seq_formatted[:-1] + self._reverse_complement(allele_b)
            pC_seq = pl_seq_formatted
            pA_primer = pr
            pB_primer = pr
            pC_primer = pl
            pA_dir = "RIGHT"
            pB_dir = "RIGHT"
            pC_dir = "LEFT"
        
        # Create output dictionaries (V2 format)
        for suffix, p_seq, primer, direction in [
            (f"Allele-{allele_a}", pA_seq, pA_primer, pA_dir),
            (f"Allele-{allele_b}", pB_seq, pB_primer, pB_dir),
            ("Common", pC_seq, pC_primer, pC_dir),
        ]:
            # Determine start/end for output (V2 uses 1-based)
            # V2 format: LEFT start < end, RIGHT start > end
            if direction == "LEFT":
                start_out = primer.start + 1
                end_out = primer.end + 1
            else:
                # For RIGHT, swap to match V2 format (start > end)
                start_out = primer.end + 1   # 5' end (higher position) 
                end_out = primer.start + 1   # 3' end (lower position)
            
            results.append({
                'index': f"{key}-{suffix}",
                'product_size': pp.product_size,
                'direction': direction,
                'start': start_out,
                'end': end_out,
                'diff_num': primer.diff_num,
                'diff_three_all': "YES" if primer.diff_three_all else "NO",
                'length': primer.length,
                'tm': primer.tm,
                'gc_percent': primer.gc_percent,
                'self_any': primer.self_any,
                'self_end': primer.self_end,
                'end_stability': primer.end_stability,
                'hairpin': primer.hairpin,
                'primer_seq': p_seq,
                'reverse_complement': self._reverse_complement(p_seq),
                'penalty': pp.penalty,
                'compl_any': pp.compl_any,
                'compl_end': pp.compl_end,
                'score': pp.score,
            })
        
        return results
    
    def _reverse_complement(self, seq: str) -> str:
        """Get reverse complement of DNA sequence."""
        complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
        return "".join(complement.get(base, base) for base in reversed(seq.upper()))
    
    def format_output(
        self,
        kasp_results: List[Dict[str, Any]],
        snp_name: str,
        output_file: Path,
        variant_sites: Optional[List[int]] = None
    ) -> None:
        """
        Format and write KASP primer results to file (V2 compatible format).
        
        Args:
            kasp_results: List of KASP primer dictionaries from design_primers
            snp_name: SNP identifier
            output_file: Output file path
            variant_sites: List of variant sites for annotation
        """
        if variant_sites is None:
            variant_sites = []
        
        try:
            with open(output_file, 'w') as f:
                # Write header (V2 format)
                header = [
                    "index", "product_size", "type", "start", "end", 
                    "variation number", "3'diffall", "length", "Tm", 
                    "GCcontent", "any", "3'", "end_stability", "hairpin", 
                    "primer_seq", "ReverseComplement", "penalty", 
                    "compl_any", "compl_end", "score"
                ]
                f.write("\t".join(header) + "\n")
                
                # Write primer results
                for result in kasp_results:
                    line = "\t".join([
                        f"{snp_name}-{result['index']}",
                        str(result['product_size']),
                        result['direction'],
                        str(result['start']),
                        str(result['end']),
                        str(result['diff_num']),
                        result['diff_three_all'],
                        str(result['length']),
                        f"{result['tm']:.2f}",
                        f"{result['gc_percent']:.2f}",
                        f"{result['self_any']:.2f}",
                        f"{result['self_end']:.2f}",
                        f"{result['end_stability']:.2f}",
                        f"{result['hairpin']:.2f}",
                        result['primer_seq'],
                        result['reverse_complement'],
                        f"{result['penalty']:.2f}" if isinstance(result['penalty'], float) else str(result['penalty']),
                        f"{result['compl_any']:.2f}" if isinstance(result['compl_any'], float) else str(result['compl_any']),
                        f"{result['compl_end']:.2f}" if isinstance(result['compl_end'], float) else str(result['compl_end']),
                        f"{result['score']:.2f}" if isinstance(result['score'], float) else str(result['score']),
                    ])
                    f.write(line + "\n")
                
                # Write variant sites information (V2 format)
                if variant_sites:
                    f.write(f"\n\nSites that can differ all for {snp_name}\n")
                    f.write(", ".join(str(x + 1) for x in variant_sites))  # Convert to 1-based
                    f.write("\n\n\n")
                    
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
        """Format a single primer line for output (legacy method)."""
        # Convert to 1-based coordinates and match V2 format
        if primer.direction == "LEFT":
            start_out = primer.start + 1
            end_out = primer.end + 1
        else:
            start_out = primer.end + 1
            end_out = primer.start + 1
            
        return "\t".join([
            index,
            str(product_size),
            primer.direction,
            str(start_out),
            str(end_out),
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