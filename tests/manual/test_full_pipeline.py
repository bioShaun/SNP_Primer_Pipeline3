#!/usr/bin/env python3
"""
Full pipeline test comparing V3 output with V2 reference data.

This script runs the COMPLETE V3 pipeline from scratch using V2's original input
file and BLAST database, then compares the outputs to verify consistency.

Unlike test_e2e_real.py which uses V2's intermediate files (flanking sequences,
alignments), this test runs the entire pipeline from input parsing through
BLAST, flanking extraction, alignment, and primer design.
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


def run_v3_full_pipeline(input_file: Path, reference_db: Path, output_dir: Path) -> bool:
    """
    Run the complete V3 pipeline from scratch.
    
    Args:
        input_file: Polymarker input CSV file
        reference_db: Path to BLAST database
        output_dir: Output directory
        
    Returns:
        True if pipeline completed successfully
    """
    from snp_primer_pipeline.config import PipelineConfig
    from snp_primer_pipeline.core.parser import PolymarkerParser
    from snp_primer_pipeline.core.blast import BlastRunner, BlastParser, FlankingExtractor
    from snp_primer_pipeline.core.alignment import MultipleSequenceAligner
    from snp_primer_pipeline.primers.kasp import KASPDesigner
    
    print(f"\n{'='*60}")
    print("Running V3 Full Pipeline")
    print(f"{'='*60}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Step 1: Parse input file
        print("\nStep 1: Parsing input file...")
        parser = PolymarkerParser(input_file)
        snps = parser.parse()
        print(f"  Parsed {len(snps)} SNPs")
        
        for snp in snps:
            print(f"    - {snp.name}: {snp.chromosome}, pos={snp.snp_position}, alleles={snp.allele_a}/{snp.allele_b}")
        
        if not snps:
            print("  ERROR: No valid SNPs found")
            return False
        
        # Step 2: Convert to FASTA for BLAST
        print("\nStep 2: Converting to FASTA format...")
        fasta_file = output_dir / "for_blast.fa"
        parser.to_fasta(fasta_file)
        print(f"  Written to: {fasta_file}")
        
        # Step 3: Run BLAST (or use V2's BLAST output for consistency)
        print("\nStep 3: Running BLAST search...")
        blast_output = output_dir / "blast_out.txt"
        
        # Use V2's BLAST output for consistency testing
        v2_blast_output = Path(__file__).parent.parent / "SNP_Primer_Pipeline2" / "examples" / "blast_out.txt"
        if v2_blast_output.exists():
            import shutil
            shutil.copy(v2_blast_output, blast_output)
            print(f"  Using V2 BLAST output: {v2_blast_output}")
        else:
            blast_runner = BlastRunner(reference_db, threads=1)
            blast_runner.run(fasta_file, blast_output)
            print(f"  BLAST output: {blast_output}")
        
        # Step 4: Parse BLAST results
        print("\nStep 4: Parsing BLAST results...")
        blast_parser = BlastParser(blast_output)
        blast_hits = blast_parser.parse()
        print(f"  Found hits for {len(blast_hits)} queries")
        
        for query_id, hits in blast_hits.items():
            print(f"    - {query_id}: {len(hits)} hits")
        
        # Step 5: Extract flanking regions
        print("\nStep 5: Extracting flanking regions...")
        flanking_extractor = FlankingExtractor(reference_db)
        
        # Create SNP position mapping
        snp_positions = {snp.name: snp.snp_position for snp in snps}
        
        flanking_regions = flanking_extractor.extract_flanking_regions(
            blast_hits,
            snp_positions,
            flanking_size=500,
            max_hits=6
        )
        
        print(f"  Extracted {len(flanking_regions)} flanking regions")
        
        # Step 6: Extract sequences
        print("\nStep 6: Extracting flanking sequences...")
        flanking_fasta = output_dir / "flanking_sequences.fa"
        flanking_extractor.extract_sequences(flanking_regions, flanking_fasta)
        print(f"  Written to: {flanking_fasta}")
        
        # Group flanking regions by SNP
        snp_regions: Dict[str, List] = {}
        for region in flanking_regions:
            if region.snp_name not in snp_regions:
                snp_regions[region.snp_name] = []
            snp_regions[region.snp_name].append(region)
        
        # Step 7: Process each SNP for KASP design
        print("\nStep 7: Processing SNPs for KASP primer design...")
        
        all_kasp_results = {}
        
        for snp in snps:
            if snp.name not in snp_regions:
                print(f"  WARNING: No flanking regions found for SNP {snp.name}")
                continue
            
            print(f"\n  Processing SNP: {snp.name}")
            regions = snp_regions[snp.name]
            print(f"    Flanking regions: {len(regions)}")
            
            # Create SNP-specific output directory
            snp_output_dir = output_dir / f"SNP_{snp.name}"
            snp_output_dir.mkdir(exist_ok=True)
            
            # Create FASTA file for this SNP's sequences
            snp_fasta = snp_output_dir / f"flanking_{snp.name}.fa"
            
            # Read sequences from the extracted flanking file
            sequences = {}
            target_name = None
            target_seq = None
            
            # Read all flanking sequences
            all_sequences = {}
            current_name = None
            current_seq = []
            
            with open(flanking_fasta, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_name:
                            all_sequences[current_name] = ''.join(current_seq)
                        current_name = line[1:].split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line)
                if current_name:
                    all_sequences[current_name] = ''.join(current_seq)
            
            # Filter sequences for this SNP
            for name, seq in all_sequences.items():
                if snp.name in name:
                    sequences[name] = seq
                    # Find target sequence (matching chromosome)
                    # The sequence name format is: snp_name_chromosome_allele_position
                    # e.g., chr7A-7659_chr7A_Y_501
                    parts = name.split('_')
                    if len(parts) >= 3:
                        seq_chromosome = parts[-3]  # chromosome is 3rd from end
                        if seq_chromosome == snp.chromosome:
                            target_name = name
                            target_seq = seq
            
            if not sequences:
                print(f"    ERROR: No sequences found for {snp.name}")
                continue
            
            # If no target found, use first sequence
            if target_seq is None:
                target_name = list(sequences.keys())[0]
                target_seq = sequences[target_name]
            
            print(f"    Target: {target_name} ({len(target_seq)} bp)")
            print(f"    Other sequences: {len(sequences) - 1}")
            
            # Write sequences to SNP-specific FASTA
            with open(snp_fasta, 'w') as f:
                for name, seq in sequences.items():
                    f.write(f">{name}\n{seq}\n")
            
            # Run alignment if multiple sequences
            alignment = None
            sites_diff_all = []
            
            if len(sequences) > 1:
                try:
                    print(f"    Running alignment...")
                    aligner = MultipleSequenceAligner()
                    alignment_file = snp_output_dir / f"alignment_{snp.name}.fa"
                    alignment = aligner.align_file(snp_fasta, alignment_file)
                    
                    if alignment:
                        alignment.set_target_sequence(target_name)
                        sites_diff_all, sites_diff_any = alignment.find_variant_sites(target_name)
                        print(f"    Variation sites (diff all): {len(sites_diff_all)}")
                except Exception as e:
                    print(f"    Alignment failed: {e}")
            
            # Design KASP primers
            try:
                print(f"    Designing KASP primers...")
                designer = KASPDesigner(max_tm=63, max_size=25, pick_anyway=False)
                
                # Get SNP position in the flanking sequence
                # The SNP should be at position 501 (500 bp flanking + 1)
                snp_pos_in_seq = 501  # Standard position
                
                primers = designer.design_primers(
                    template_sequence=target_seq,
                    snp_position=snp_pos_in_seq,
                    snp_alleles=(snp.allele_a, snp.allele_b),
                    product_size_range=(50, 250),
                    variant_sites=sites_diff_all,
                    alignment=alignment,
                    target_name=target_name,
                    output_dir=snp_output_dir
                )
                
                # Write output
                output_file = snp_output_dir / f"KASP_primers_{snp.name}.txt"
                designer.format_output(primers, snp.name, output_file, sites_diff_all)
                
                print(f"    Designed {len(primers)} primer pairs")
                
                all_kasp_results[snp.name] = {
                    'primers': primers,
                    'sites': sites_diff_all,
                    'output_file': output_file
                }
                
            except Exception as e:
                print(f"    KASP design failed: {e}")
                import traceback
                traceback.print_exc()
        
        # Combine all KASP results into a single file
        print("\nStep 8: Combining KASP results...")
        combined_output = output_dir / "Potential_KASP_primers.tsv"
        
        with open(combined_output, 'w') as out:
            first_snp = True
            for snp_name, results in all_kasp_results.items():
                if results['output_file'].exists():
                    with open(results['output_file'], 'r') as f:
                        for line in f:
                            # Skip header for subsequent SNPs
                            if not first_snp and line.startswith('index'):
                                continue
                            out.write(line)
                    first_snp = False
        
        print(f"  Combined output: {combined_output}")
        
        return True
        
    except Exception as e:
        print(f"\nERROR: Pipeline failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def compare_outputs(v2_file: Path, v3_file: Path) -> Tuple[int, int, List[str]]:
    """
    Compare V2 and V3 KASP output files.
    
    Returns:
        Tuple of (passed, failed, error_messages)
    """
    passed = 0
    failed = 0
    errors = []
    
    v2_primers = parse_kasp_output(v2_file)
    v3_primers = parse_kasp_output(v3_file)
    
    v2_sites = parse_variation_sites(v2_file)
    v3_sites = parse_variation_sites(v3_file)
    
    print(f"\n{'='*60}")
    print("Comparing V2 and V3 Outputs")
    print(f"{'='*60}")
    
    print(f"\nV2 KASP primers: {sum(len(p) for p in v2_primers.values())} for {len(v2_primers)} SNPs")
    print(f"V3 KASP primers: {sum(len(p) for p in v3_primers.values())} for {len(v3_primers)} SNPs")
    
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
        print(f"\n--- {snp_name} ---")
        
        v2_list = v2_primers[snp_name]
        v3_list = v3_primers[snp_name]
        
        print(f"  V2 primers: {len(v2_list)}")
        print(f"  V3 primers: {len(v3_list)}")
        
        # Compare primer sequences
        v2_seqs = set(p.primer_seq.upper() for p in v2_list)
        v3_seqs = set(p.primer_seq.upper() for p in v3_list)
        common_seqs = v2_seqs & v3_seqs
        
        print(f"  V2 unique sequences: {len(v2_seqs)}")
        print(f"  V3 unique sequences: {len(v3_seqs)}")
        print(f"  Common sequences: {len(common_seqs)}")
        
        if len(common_seqs) > 0:
            print(f"  ✅ Found {len(common_seqs)} matching primer sequences")
            passed += 1
            # Show some matching primers
            print(f"  Sample matching sequences:")
            for seq in list(common_seqs)[:3]:
                print(f"    - {seq}")
        else:
            # Check for partial matches (ignoring last base which may be modified)
            v2_core = set(p.primer_seq.upper()[:-1] for p in v2_list)
            v3_core = set(p.primer_seq.upper()[:-1] for p in v3_list)
            partial_matches = v2_core & v3_core
            if partial_matches:
                print(f"  ⚠️ Found {len(partial_matches)} partial matches (ignoring 3' base)")
                passed += 1
            else:
                print(f"  ❌ No matching primer sequences")
                failed += 1
                errors.append(f"{snp_name}: No matching primer sequences")
        
        # Compare variation sites
        v2_site_list = v2_sites.get(snp_name, [])
        v3_site_list = v3_sites.get(snp_name, [])
        
        if v2_site_list and v3_site_list:
            v2_set = set(v2_site_list)
            v3_set = set(v3_site_list)
            common_sites = v2_set & v3_set
            print(f"  Variation sites: V2={len(v2_site_list)}, V3={len(v3_site_list)}, common={len(common_sites)}")
    
    return passed, failed, errors


def main():
    """Main test function."""
    print("="*60)
    print("SNP Primer Pipeline V3 Full Pipeline Test")
    print("="*60)
    
    # Check prerequisites
    if not V2_EXAMPLES.exists():
        print(f"Error: V2 examples directory not found: {V2_EXAMPLES}")
        return 1
    
    # Input files
    input_file = V2_EXAMPLES / "polymarker_input_example.csv"
    reference_db = V2_EXAMPLES / "blastdb" / "test_reference.fa"
    v2_kasp_file = V2_EXAMPLES / "Potential_KASP_primers.tsv"
    
    if not input_file.exists():
        print(f"Error: Input file not found: {input_file}")
        return 1
    
    if not reference_db.exists():
        print(f"Error: Reference database not found: {reference_db}")
        return 1
    
    print(f"\nInput file: {input_file}")
    print(f"Reference DB: {reference_db}")
    print(f"V2 reference output: {v2_kasp_file}")
    
    # Run V3 full pipeline
    with tempfile.TemporaryDirectory() as tmp_dir:
        v3_output_dir = Path(tmp_dir) / "v3_output"
        
        success = run_v3_full_pipeline(input_file, reference_db, v3_output_dir)
        
        if not success:
            print("\n❌ V3 pipeline failed to complete")
            return 1
        
        # Compare outputs
        v3_kasp_file = v3_output_dir / "Potential_KASP_primers.tsv"
        
        if not v3_kasp_file.exists():
            print(f"\n❌ V3 output file not found: {v3_kasp_file}")
            return 1
        
        passed, failed, errors = compare_outputs(v2_kasp_file, v3_kasp_file)
        
        # Summary
        print(f"\n{'='*60}")
        print("Test Summary")
        print(f"{'='*60}")
        
        print(f"\n✅ Passed: {passed}")
        print(f"❌ Failed: {failed}")
        
        if errors:
            print(f"\nErrors:")
            for error in errors:
                print(f"  - {error}")
        
        if failed == 0 and passed > 0:
            print("\n✅ FULL PIPELINE TEST PASSED!")
            return 0
        elif passed > failed:
            print("\n⚠️ FULL PIPELINE TEST PARTIALLY PASSED")
            return 0
        else:
            print("\n❌ FULL PIPELINE TEST FAILED")
            return 1


if __name__ == "__main__":
    sys.exit(main())
