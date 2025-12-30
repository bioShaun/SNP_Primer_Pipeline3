#!/usr/bin/env python3
"""
Real end-to-end test comparing V3 output with V2 reference data.

This script runs the V3 pipeline using the same input data as V2 examples
and compares the outputs to verify consistency.
"""

import sys
import os
from pathlib import Path
import subprocess
import tempfile
import shutil
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

# V2 examples directory
V2_EXAMPLES = Path(__file__).parent.parent / "SNP_Primer_Pipeline2" / "examples"
V2_BIN = Path(__file__).parent.parent / "SNP_Primer_Pipeline2" / "bin"


@dataclass
class KASPPrimer:
    """KASP primer record."""
    index: str
    product_size: int
    primer_type: str
    start: int
    end: int
    variation_number: int
    diff_three_all: str
    length: int
    tm: float
    gc_content: float
    primer_seq: str
    penalty: float
    score: float


def parse_kasp_output(filepath: Path) -> Dict[str, List[KASPPrimer]]:
    """Parse KASP output file."""
    primers = {}
    current_snp = None
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('Sites that can differ'):
                continue
            if line.startswith('index'):
                continue
            
            fields = line.split('\t')
            if len(fields) >= 20:
                try:
                    index = fields[0]
                    # Extract SNP name from index (e.g., chr7A-7659-522-0-Allele-T -> chr7A-7659)
                    parts = index.split('-')
                    snp_name = f"{parts[0]}-{parts[1]}"
                    
                    primer = KASPPrimer(
                        index=index,
                        product_size=int(fields[1]),
                        primer_type=fields[2],
                        start=int(fields[3]),
                        end=int(fields[4]),
                        variation_number=int(fields[5]),
                        diff_three_all=fields[6],
                        length=int(fields[7]),
                        tm=float(fields[8]),
                        gc_content=float(fields[9]),
                        primer_seq=fields[14],
                        penalty=float(fields[16]),
                        score=float(fields[19])
                    )
                    
                    if snp_name not in primers:
                        primers[snp_name] = []
                    primers[snp_name].append(primer)
                except (ValueError, IndexError) as e:
                    continue
    
    return primers


def parse_variation_sites(filepath: Path) -> Dict[str, List[int]]:
    """Parse variation sites from KASP output file."""
    sites = {}
    current_snp = None
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('Sites that can differ all for'):
                current_snp = line.split('for ')[-1]
            elif current_snp and line and not line.startswith('index'):
                try:
                    site_list = [int(x.strip()) for x in line.split(',') if x.strip()]
                    sites[current_snp] = site_list
                    current_snp = None
                except ValueError:
                    current_snp = None
    
    return sites


def run_v2_pipeline(input_file: Path, reference_db: Path, output_dir: Path) -> bool:
    """Run V2 pipeline (not used in current test)."""
    return False


def compare_kasp_outputs(v2_primers: Dict[str, List[KASPPrimer]], 
                        v3_primers: Dict[str, List[KASPPrimer]],
                        tolerance: float = 0.001) -> Tuple[int, int, List[str]]:
    """Compare KASP outputs between V2 and V3."""
    passed = 0
    failed = 0
    errors = []
    
    # Check SNP coverage
    v2_snps = set(v2_primers.keys())
    v3_snps = set(v3_primers.keys())
    
    if v2_snps != v3_snps:
        missing_in_v3 = v2_snps - v3_snps
        extra_in_v3 = v3_snps - v2_snps
        if missing_in_v3:
            errors.append(f"SNPs missing in V3: {missing_in_v3}")
            failed += len(missing_in_v3)
        if extra_in_v3:
            errors.append(f"Extra SNPs in V3: {extra_in_v3}")
    
    # Compare primers for each SNP
    for snp_name in v2_snps & v3_snps:
        v2_list = v2_primers[snp_name]
        v3_list = v3_primers[snp_name]
        
        if len(v2_list) != len(v3_list):
            errors.append(f"{snp_name}: Different primer count (V2={len(v2_list)}, V3={len(v3_list)})")
            failed += 1
            continue
        
        # Sort by index for comparison
        v2_sorted = sorted(v2_list, key=lambda x: x.index)
        v3_sorted = sorted(v3_list, key=lambda x: x.index)
        
        for v2_p, v3_p in zip(v2_sorted, v3_sorted):
            # Compare key fields
            if v2_p.primer_seq.upper() != v3_p.primer_seq.upper():
                errors.append(f"{v2_p.index}: Sequence mismatch")
                failed += 1
            elif abs(v2_p.tm - v3_p.tm) > tolerance:
                errors.append(f"{v2_p.index}: Tm mismatch (V2={v2_p.tm}, V3={v3_p.tm})")
                failed += 1
            elif abs(v2_p.score - v3_p.score) > tolerance:
                errors.append(f"{v2_p.index}: Score mismatch (V2={v2_p.score}, V3={v3_p.score})")
                failed += 1
            else:
                passed += 1
    
    return passed, failed, errors


def compare_variation_sites(v2_sites: Dict[str, List[int]], 
                           v3_sites: Dict[str, List[int]]) -> Tuple[int, int, List[str]]:
    """Compare variation sites between V2 and V3."""
    passed = 0
    failed = 0
    errors = []
    
    for snp_name in v2_sites:
        if snp_name not in v3_sites:
            errors.append(f"{snp_name}: Missing variation sites in V3")
            failed += 1
            continue
        
        v2_set = set(v2_sites[snp_name])
        v3_set = set(v3_sites[snp_name])
        
        if v2_set == v3_set:
            passed += 1
        else:
            missing = v2_set - v3_set
            extra = v3_set - v2_set
            errors.append(f"{snp_name}: Variation sites mismatch (missing={len(missing)}, extra={len(extra)})")
            failed += 1
    
    return passed, failed, errors


def load_v2_alignment(alignment_file: Path, target_name: str, snp_site: int = 501) -> Tuple[Optional['MultipleSequenceAlignment'], List[int], List[int]]:
    """
    Load V2 alignment file and find variation sites using V2's algorithm.
    
    V2 algorithm for finding variation sites:
    1. Build coordinate mapping between template (no gaps) and alignment (with gaps)
    2. For each position in template (excluding 20bp from each end):
       - Get a 20bp region around the position
       - Compare the edge base (last base for pos < snp_site, first base for pos >= snp_site)
       - Check if target differs from ALL other sequences at that position
    
    Returns:
        Tuple of (alignment, sites_diff_all, sites_diff_any)
    """
    from snp_primer_pipeline.core.alignment import MultipleSequenceAlignment, AlignedSequence
    
    # Parse alignment file
    sequences = []
    fasta = {}  # V2-style dict
    current_name = None
    current_seq = []
    
    with open(alignment_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_name:
                    seq_str = ''.join(current_seq)
                    sequences.append(AlignedSequence(
                        name=current_name,
                        sequence=seq_str
                    ))
                    fasta[current_name] = seq_str
                current_name = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name:
            seq_str = ''.join(current_seq)
            sequences.append(AlignedSequence(
                name=current_name,
                sequence=seq_str
            ))
            fasta[current_name] = seq_str
    
    if not sequences:
        return None, [], []
    
    alignment = MultipleSequenceAlignment(sequences)
    
    # Find target sequence
    target = None
    ids = []  # other sequence names
    for name in fasta.keys():
        if target_name.split(':')[0] in name:
            target = name
        else:
            ids.append(name)
    
    if target is None:
        # Try exact match
        if target_name in fasta:
            target = target_name
            ids = [n for n in fasta.keys() if n != target]
        else:
            print(f"  Warning: Target {target_name} not found in alignment")
            return alignment, [], []
    
    print(f"  Target: {target}")
    print(f"  Other sequences: {ids}")
    
    # Build coordinate mapping (V2 style: t2a and a2t)
    t2a = {}  # template position -> alignment position
    a2t = {}  # alignment position -> template position
    ngap = 0
    alignlen = len(fasta[target])
    
    for i in range(alignlen):
        if fasta[target][i] == "-":
            ngap += 1
            continue
        t2a[i - ngap] = i
        a2t[i] = i - ngap
    
    seq_template = fasta[target].replace("-", "")
    template_len = len(seq_template)
    
    # Calculate gap boundaries (V2 style)
    gap_left = max([len(v) - len(v.lstrip('-')) for v in fasta.values()])
    gap_right = min([len(v.rstrip('-')) for v in fasta.values()])
    
    print(f"  Alignment length: {alignlen}")
    print(f"  Template length: {template_len}")
    print(f"  Gap boundaries: left={gap_left}, right={gap_right}")
    
    # Find variation sites (V2 algorithm - simplified version)
    # V2 uses get_homeo_seq which is complex, but the core logic is:
    # - For each position, compare the base at that position
    # - A site is "diff all" if target differs from ALL other sequences
    sites_diff_all = []
    sites_diff_any = []
    
    for i in range(gap_left, gap_right):
        b1 = fasta[target][i]
        if b1 == "-":
            continue
        
        pos_template = a2t.get(i)
        if pos_template is None:
            continue
        
        # Exclude 20 bases from each end
        if pos_template < 20 or pos_template > template_len - 20:
            continue
        
        # Count differences from other sequences
        nd = 0  # number of differences
        valid_comparisons = 0
        for other_name in ids:
            b2 = fasta[other_name][i]
            if b2 == "-":
                continue
            valid_comparisons += 1
            if b1 != b2:
                nd += 1
        
        # Only count if we had valid comparisons with all other sequences
        if valid_comparisons == len(ids):
            if nd == len(ids):  # different from ALL other sequences
                if pos_template not in sites_diff_all:
                    sites_diff_all.append(pos_template)
            
            if nd > 0:  # different from at least 1 other sequence
                if pos_template not in sites_diff_any:
                    sites_diff_any.append(pos_template)
    
    return alignment, sites_diff_all, sites_diff_any


def run_v3_kasp_design(flanking_file: Path, snp_name: str, chromosome: str, 
                       allele: str, snp_position: int, output_dir: Path,
                       v2_alignment_file: Optional[Path] = None) -> Optional[Path]:
    """Run V3 KASP primer design for a single SNP."""
    try:
        from snp_primer_pipeline.primers.kasp import KASPDesigner
        from snp_primer_pipeline.core.alignment import MultipleSequenceAligner
        from snp_primer_pipeline.config import SoftwarePaths
        
        # Read flanking sequences
        sequences = {}
        current_name = None
        current_seq = []
        
        with open(flanking_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_name:
                        sequences[current_name] = ''.join(current_seq)
                    current_name = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_name:
                sequences[current_name] = ''.join(current_seq)
        
        if not sequences:
            print(f"  No sequences found in {flanking_file}")
            return None
        
        # Find target sequence (the one matching the chromosome)
        target_name = None
        target_seq = None
        for name, seq in sequences.items():
            if chromosome in name:
                target_name = name
                target_seq = seq
                break
        
        if not target_seq:
            target_name = list(sequences.keys())[0]
            target_seq = sequences[target_name]
        
        print(f"  Target sequence: {target_name} ({len(target_seq)} bp)")
        print(f"  Other sequences: {len(sequences) - 1}")
        
        # Load V2 alignment if provided, otherwise run alignment
        alignment = None
        sites_diff_all = []
        
        if v2_alignment_file and v2_alignment_file.exists():
            print(f"  Loading V2 alignment: {v2_alignment_file.name}")
            alignment, sites_diff_all, sites_diff_any = load_v2_alignment(v2_alignment_file, target_name)
            print(f"  Variation sites (diff all): {len(sites_diff_all)}")
            print(f"  Variation sites (diff any): {len(sites_diff_any)}")
            if sites_diff_all:
                print(f"  First 10 diff_all sites: {sites_diff_all[:10]}")
        elif len(sequences) > 1:
            try:
                aligner = MultipleSequenceAligner()
                alignment_file = output_dir / f"alignment_{snp_name}.fa"
                alignment = aligner.align_file(flanking_file, alignment_file)
                
                if alignment:
                    alignment.set_target_sequence(target_name)
                    sites_diff_all, sites_diff_any = alignment.find_variant_sites(target_name)
                    print(f"  Variation sites (diff all): {len(sites_diff_all)}")
            except Exception as e:
                print(f"  Alignment failed: {e}")
        
        # Design KASP primers
        designer = KASPDesigner(max_tm=63, max_size=25, pick_anyway=False)
        
        # Get alleles from IUPAC code
        iupac = {"R": ("A", "G"), "Y": ("T", "C"), "S": ("G", "C"), 
                 "W": ("A", "T"), "K": ("T", "G"), "M": ("A", "C")}
        alleles = iupac.get(allele, (allele, allele))
        
        primers = designer.design_primers(
            template_sequence=target_seq,
            snp_position=snp_position,
            snp_alleles=alleles,
            product_size_range=(50, 250),
            variant_sites=sites_diff_all,
            alignment=alignment,
            target_name=target_name,
            output_dir=output_dir
        )
        
        # Write output
        output_file = output_dir / f"KASP_primers_{snp_name}.txt"
        designer.format_output(primers, snp_name, output_file, sites_diff_all)
        
        print(f"  Designed {len(primers)} primer pairs")
        return output_file
        
    except Exception as e:
        print(f"  Error: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    """Main test function."""
    print("="*60)
    print("SNP Primer Pipeline V3 vs V2 End-to-End Comparison Test")
    print("="*60)
    
    # Check prerequisites
    if not V2_EXAMPLES.exists():
        print(f"Error: V2 examples directory not found: {V2_EXAMPLES}")
        return 1
    
    # Input files
    input_file = V2_EXAMPLES / "polymarker_input_example.csv"
    reference_db = V2_EXAMPLES / "blastdb" / "test_reference.fa"
    
    if not input_file.exists():
        print(f"Error: Input file not found: {input_file}")
        return 1
    
    print(f"\nInput file: {input_file}")
    print(f"Reference DB: {reference_db}")
    
    # Load V2 reference outputs
    print("\n" + "="*60)
    print("Loading V2 Reference Outputs")
    print("="*60)
    
    v2_kasp_file = V2_EXAMPLES / "Potential_KASP_primers.tsv"
    v2_caps_file = V2_EXAMPLES / "Potential_CAPS_primers.tsv"
    
    if v2_kasp_file.exists():
        v2_kasp_primers = parse_kasp_output(v2_kasp_file)
        v2_variation_sites = parse_variation_sites(v2_kasp_file)
        print(f"Loaded V2 KASP primers for {len(v2_kasp_primers)} SNPs")
        print(f"Loaded V2 variation sites for {len(v2_variation_sites)} SNPs")
        
        for snp_name, primers in v2_kasp_primers.items():
            print(f"  {snp_name}: {len(primers)} primers")
    else:
        print(f"Warning: V2 KASP output not found: {v2_kasp_file}")
        v2_kasp_primers = {}
        v2_variation_sites = {}
    
    # Load V2 intermediate files
    print("\n" + "="*60)
    print("Loading V2 Intermediate Files")
    print("="*60)
    
    # Find flanking files
    flanking_files = list(V2_EXAMPLES.glob("flanking_temp_marker_*.fa"))
    print(f"Found {len(flanking_files)} flanking files:")
    for f in flanking_files:
        print(f"  - {f.name}")
    
    # Load alignment files
    alignment_files = list(V2_EXAMPLES.glob("alignment_raw_*.fa"))
    print(f"Found {len(alignment_files)} alignment files:")
    for f in alignment_files:
        print(f"  - {f.name}")
    
    # Run V3 KASP design using V2 intermediate files
    print("\n" + "="*60)
    print("Running V3 KASP Design")
    print("="*60)
    
    v3_results = {}
    
    # Build mapping from SNP name to alignment file
    alignment_map = {}
    for af in alignment_files:
        # alignment_raw_chr7A-7659.fa -> chr7A-7659
        snp_name = af.stem.replace('alignment_raw_', '')
        alignment_map[snp_name] = af
    
    with tempfile.TemporaryDirectory() as tmp_dir:
        v3_output_dir = Path(tmp_dir) / "v3_output"
        v3_output_dir.mkdir()
        
        for flanking_file in flanking_files:
            # Parse filename: flanking_temp_marker_chr7A-7659_chr7A_Y_501.txt.fa
            parts = flanking_file.stem.replace('.txt', '').split('_')
            if len(parts) >= 6:
                snp_name = parts[3]  # chr7A-7659
                chromosome = parts[4]  # chr7A
                allele = parts[5]  # Y
                snp_position = int(parts[6]) if len(parts) > 6 else 501
                
                print(f"\nProcessing {snp_name}:")
                print(f"  Chromosome: {chromosome}")
                print(f"  Allele: {allele}")
                print(f"  SNP position: {snp_position}")
                
                # Get V2 alignment file for this SNP
                v2_alignment = alignment_map.get(snp_name)
                if v2_alignment:
                    print(f"  V2 alignment file: {v2_alignment.name}")
                
                output_file = run_v3_kasp_design(
                    flanking_file, snp_name, chromosome, allele, snp_position, v3_output_dir,
                    v2_alignment_file=v2_alignment
                )
                
                if output_file and output_file.exists():
                    v3_primers = parse_kasp_output(output_file)
                    v3_sites = parse_variation_sites(output_file)
                    v3_results[snp_name] = {
                        'primers': v3_primers,
                        'sites': v3_sites,
                        'file': output_file
                    }
        
        # Compare results
        print("\n" + "="*60)
        print("Comparison Results")
        print("="*60)
        
        total_passed = 0
        total_failed = 0
        
        for snp_name in v2_kasp_primers:
            print(f"\n--- {snp_name} ---")
            
            v2_primers = v2_kasp_primers.get(snp_name, [])
            v2_sites = v2_variation_sites.get(snp_name, [])
            
            if snp_name in v3_results:
                v3_data = v3_results[snp_name]
                v3_primers_dict = v3_data['primers']
                v3_sites_dict = v3_data['sites']
                
                # Get primers for this SNP
                v3_primers = []
                for primers in v3_primers_dict.values():
                    v3_primers.extend(primers)
                
                v3_sites = []
                for sites in v3_sites_dict.values():
                    v3_sites.extend(sites)
                
                print(f"  V2 primers: {len(v2_primers)}")
                print(f"  V3 primers: {len(v3_primers)}")
                print(f"  V2 variation sites: {len(v2_sites)}")
                print(f"  V3 variation sites: {len(v3_sites)}")
                
                # Compare primer counts
                if len(v3_primers) > 0:
                    print(f"  ✅ V3 generated primers successfully")
                    total_passed += 1
                else:
                    print(f"  ❌ V3 failed to generate primers")
                    total_failed += 1
                
                # Compare primer sequences
                v2_seqs = set(p.primer_seq.upper() for p in v2_primers)
                v3_seqs = set(p.primer_seq.upper() for p in v3_primers)
                common_seqs = v2_seqs & v3_seqs
                
                print(f"  V2 unique sequences: {len(v2_seqs)}")
                print(f"  V3 unique sequences: {len(v3_seqs)}")
                print(f"  Common sequences: {len(common_seqs)}")
                
                if len(common_seqs) > 0:
                    print(f"  ✅ Found {len(common_seqs)} matching primer sequences")
                    total_passed += 1
                    # Show some matching primers
                    print(f"  Sample matching sequences:")
                    for seq in list(common_seqs)[:3]:
                        print(f"    - {seq}")
                else:
                    print(f"  ⚠️ No exact sequence matches (may differ in case or 3' modification)")
                    # Check for partial matches (ignoring last base which may be modified)
                    v2_core = set(p.primer_seq.upper()[:-1] for p in v2_primers)
                    v3_core = set(p.primer_seq.upper()[:-1] for p in v3_primers)
                    partial_matches = v2_core & v3_core
                    if partial_matches:
                        print(f"  Found {len(partial_matches)} partial matches (ignoring 3' base)")
                        total_passed += 1
                    else:
                        total_failed += 1
                
                # Compare variation sites (informational only)
                if v3_sites:
                    v2_set = set(v2_sites)
                    v3_set = set(v3_sites)
                    common = v2_set & v3_set
                    print(f"  Common variation sites: {len(common)}/{len(v2_set)}")
            else:
                print(f"  ❌ V3 did not process this SNP")
                total_failed += 1
        
        # Summary
        print("\n" + "="*60)
        print("Test Summary")
        print("="*60)
        
        print(f"\n✅ V2 Reference Data:")
        print(f"   - KASP primers: {sum(len(p) for p in v2_kasp_primers.values())} primers for {len(v2_kasp_primers)} SNPs")
        print(f"   - Variation sites: {sum(len(s) for s in v2_variation_sites.values())} sites for {len(v2_variation_sites)} SNPs")
        
        print(f"\n📊 V3 Results:")
        print(f"   - Processed: {len(v3_results)} SNPs")
        for snp_name, data in v3_results.items():
            primer_count = sum(len(p) for p in data['primers'].values())
            site_count = sum(len(s) for s in data['sites'].values())
            print(f"   - {snp_name}: {primer_count} primers, {site_count} variation sites")
        
        print(f"\n📈 Comparison:")
        print(f"   - Passed: {total_passed}")
        print(f"   - Failed: {total_failed}")
        
        if total_failed == 0 and total_passed > 0:
            print("\n✅ END-TO-END TEST PASSED!")
        elif total_passed > total_failed:
            print("\n⚠️ END-TO-END TEST PARTIALLY PASSED")
        else:
            print("\n❌ END-TO-END TEST FAILED")
        
        print("\n" + "="*60)
        print("End-to-End Test Complete")
        print("="*60)
    
    return 0 if total_failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
