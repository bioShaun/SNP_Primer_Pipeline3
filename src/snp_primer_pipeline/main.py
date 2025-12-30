#!/usr/bin/env python3
"""
Main pipeline module for SNP Primer Pipeline.

This module provides the main entry point and orchestrates the entire pipeline.
"""

import logging
import sys
from pathlib import Path
from typing import List, Optional, Dict, Any
import argparse

from .config import PipelineConfig, SoftwarePaths
from .core.parser import PolymarkerParser
from .core.blast import BlastRunner, BlastParser, FlankingExtractor
from .core.alignment import MultipleSequenceAligner
from .primers.kasp import KASPDesigner
from .primers.caps import CAPSDesigner
from .exceptions import PipelineError, ParseError, BlastError, AlignmentError, PrimerDesignError


def setup_logging(log_level: str = "INFO") -> None:
    """Setup logging configuration."""
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
        ]
    )


def run_pipeline(config: PipelineConfig) -> None:
    """
    Run the complete SNP primer design pipeline.
    
    Args:
        config: Pipeline configuration
    """
    logger = logging.getLogger(__name__)
    
    # Create output directory
    config.output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Starting SNP Primer Pipeline")
    logger.info(f"Input file: {config.input_file}")
    logger.info(f"Reference: {config.reference_file}")
    logger.info(f"Output directory: {config.output_dir}")
    
    try:
        # Step 1: Parse input file
        logger.info("Step 1: Parsing input file...")
        parser = PolymarkerParser(config.input_file)
        snps = parser.parse()
        logger.info(f"Parsed {len(snps)} SNPs")
        
        if not snps:
            logger.warning("No valid SNPs found in input file")
            return
        
        # Step 2: Convert to FASTA for BLAST
        logger.info("Step 2: Converting to FASTA format...")
        fasta_file = config.output_dir / "for_blast.fa"
        parser.to_fasta(fasta_file)
        
        # Step 3: Run BLAST
        logger.info("Step 3: Running BLAST search...")
        blast_runner = BlastRunner(config.reference_file, config.threads)
        blast_output = config.output_dir / "blast_out.txt"
        blast_runner.run(fasta_file, blast_output)
        
        # Step 4: Parse BLAST results
        logger.info("Step 4: Parsing BLAST results...")
        blast_parser = BlastParser(blast_output)
        blast_hits = blast_parser.parse()
        
        # Step 5: Extract flanking regions
        logger.info("Step 5: Extracting flanking regions...")
        flanking_extractor = FlankingExtractor(config.reference_file)
        
        # Create SNP position mapping
        snp_positions = {snp.name: snp.snp_position for snp in snps}
        
        flanking_regions = flanking_extractor.extract_flanking_regions(
            blast_hits,
            snp_positions,
            config.flanking_size,
            config.max_hits
        )
        
        logger.info(f"Extracted {len(flanking_regions)} flanking regions")
        
        # Step 6: Extract sequences
        logger.info("Step 6: Extracting flanking sequences...")
        flanking_fasta = config.output_dir / "flanking_sequences.fa"
        flanking_extractor.extract_sequences(flanking_regions, flanking_fasta)
        
        # Load extracted sequences
        extracted_sequences = {}
        current_id = None
        current_seq = []
        
        if flanking_fasta.exists():
            with open(flanking_fasta, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_id:
                            extracted_sequences[current_id] = ''.join(current_seq)
                        current_id = line[1:].split()[0]  # Take ID up to first whitespace
                        current_seq = []
                    else:
                        current_seq.append(line)
                if current_id:
                    extracted_sequences[current_id] = ''.join(current_seq)
        
        # Step 7: Process each SNP
        logger.info("Step 7: Processing SNPs for primer design...")
        
        # Group flanking regions by SNP
        snp_regions: Dict[str, List] = {}
        for region in flanking_regions:
            if region.snp_name not in snp_regions:
                snp_regions[region.snp_name] = []
            snp_regions[region.snp_name].append(region)
        
        # Process each SNP
        for snp in snps:
            if snp.name not in snp_regions:
                logger.warning(f"No flanking regions found for SNP {snp.name}")
                continue
            
            logger.info(f"Processing SNP: {snp.name}")
            
            try:
                process_snp(
                    snp,
                    snp_regions[snp.name],
                    extracted_sequences,
                    config
                )
            except Exception as e:
                logger.error(f"Failed to process SNP {snp.name}: {e}")
                continue
        
        logger.info("Pipeline completed successfully")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        raise


def process_snp(
    snp,
    flanking_regions: List,
    extracted_sequences: Dict[str, str],
    config: PipelineConfig
) -> None:
    """
    Process a single SNP for primer design.
    
    Args:
        snp: SNP object
        flanking_regions: List of flanking regions for this SNP
        extracted_sequences: Dictionary of extracted sequences by ID
        config: Pipeline configuration
    """
    logger = logging.getLogger(__name__)
    
    # Create SNP-specific output directory
    snp_output_dir = config.output_dir / f"SNP_{snp.name}"
    snp_output_dir.mkdir(exist_ok=True)
    
    # Create FASTA file for this SNP's sequences
    snp_fasta = snp_output_dir / f"flanking_{snp.name}.fa"
    
    # Write sequences to FASTA
    sequences = {}
    target_sequence = None
    target_snp_position = None
    genomic_start = None
    genomic_strand = None
    
    with open(snp_fasta, 'w') as f:
        for region in flanking_regions:
            seq_id = f"{region.snp_name}_{region.chromosome}_{region.allele}_{region.snp_position_in_region}"
            
            if seq_id not in extracted_sequences:
                logger.warning(f"Sequence not found for {seq_id}, skipping")
                continue
                
            sequence = extracted_sequences[seq_id]
            sequences[seq_id] = sequence
            f.write(f">{seq_id}\n{sequence}\n")
            
            # Identify target sequence (first one for now)
            if target_sequence is None:
                target_sequence = sequence
                target_snp_position = region.snp_position_in_region - 1
                genomic_start = region.start
                genomic_strand = region.strand
    
    if len(sequences) < 2:
        logger.warning(f"Not enough sequences for alignment for SNP {snp.name}")
        return
    
    # Perform multiple sequence alignment
    alignment = None
    target_name = None
    
    try:
        aligner = MultipleSequenceAligner()
        alignment = aligner.align_file(
            snp_fasta,
            snp_output_dir / f"alignment_{snp.name}.fa"
        )
        
        # Find variant sites
        target_name = list(sequences.keys())[0]
        alignment.set_target_sequence(target_name)
        
        # Use V2-style variant site detection with SNP position awareness
        sites_diff_all, sites_diff_any, diffarray = alignment.find_variant_sites_v2(
            target_name, target_snp_position
        )
        
        logger.info(f"Found {len(sites_diff_all)} sites that differ from all homeologs")
        
    except AlignmentError as e:
        logger.warning(f"Alignment failed for SNP {snp.name}: {e}")
        sites_diff_all = []
        sites_diff_any = []
        diffarray = {}
    
    # Design KASP primers
    if config.design_kasp:
        try:
            logger.info(f"Designing KASP primers for {snp.name}")
            
            kasp_designer = KASPDesigner(
                max_tm=config.max_tm,
                max_size=config.max_primer_size,
                pick_anyway=config.pick_anyway
            )
            
            kasp_primers = kasp_designer.design_primers(
                target_sequence,
                target_snp_position,
                (snp.allele_a, snp.allele_b),
                config.primer_product_size_range,
                sites_diff_all,
                alignment if len(sequences) > 1 else None,
                target_name if len(sequences) > 1 else None,
                snp_output_dir,
                diffarray,  # Pass diffarray for V2-style filtering
                genomic_start=genomic_start,
                genomic_strand=genomic_strand
            )
            
            # Write KASP results
            kasp_output = snp_output_dir / f"KASP_primers_{snp.name}.txt"
            kasp_designer.format_output(
                kasp_primers, snp.name, kasp_output, sites_diff_all,
                show_variant_sites=config.show_variant_sites
            )
            
            # Write simplified summary
            kasp_summary = snp_output_dir / f"KASP_primers_{snp.name}_summary.txt"
            kasp_designer.format_simple_output(kasp_primers, snp.name, kasp_summary)
            
            logger.info(f"Designed {len(kasp_primers) // 3} KASP primer pairs for {snp.name}")
            
        except PrimerDesignError as e:
            logger.warning(f"KASP design failed for SNP {snp.name}: {e}")
    
    # Design CAPS primers
    if config.design_caps:
        try:
            logger.info(f"Designing CAPS primers for {snp.name}")
            
            # Load enzyme database
            enzyme_file = Path(__file__).parent.parent.parent / "resources" / "NEB_parsed_REs.txt"
            
            caps_designer = CAPSDesigner(
                enzyme_file=enzyme_file,
                max_tm=config.max_tm,
                max_size=config.max_primer_size,
                pick_anyway=config.pick_anyway
            )
            
            # Find usable enzymes
            caps_enzymes, dcaps_enzymes = caps_designer.find_usable_enzymes(
                target_sequence,
                target_snp_position,
                (snp.allele_a, snp.allele_b),
                config.max_price
            )
            
            logger.info(f"Found {len(caps_enzymes)} CAPS and {len(dcaps_enzymes)} dCAPS enzymes")
            
            # Design primers for each enzyme
            all_caps_primers = []
            
            for enzyme in caps_enzymes + dcaps_enzymes:
                try:
                    primers = caps_designer.design_caps_primers(
                        target_sequence,
                        target_snp_position,
                        (snp.allele_a, snp.allele_b),
                        enzyme,
                        (300, 900),  # CAPS product size range
                        sites_diff_all,
                        snp_output_dir
                    )
                    all_caps_primers.extend(primers)
                except PrimerDesignError:
                    continue
            
            # Write CAPS results
            caps_output = snp_output_dir / f"CAPS_primers_{snp.name}.txt"
            caps_designer.format_output(
                all_caps_primers, 
                snp.name, 
                caps_output,
                caps_enzymes,
                dcaps_enzymes,
                sites_diff_all
            )
            
            logger.info(f"Designed {len(all_caps_primers)} CAPS primer pairs for {snp.name}")
            
        except PrimerDesignError as e:
            logger.warning(f"CAPS design failed for SNP {snp.name}: {e}")


def main() -> None:
    """Main entry point for command line interface."""
    parser = argparse.ArgumentParser(
        description="SNP Primer Pipeline - Design KASP and CAPS primers for SNPs"
    )
    
    parser.add_argument(
        "input_file",
        type=Path,
        help="Input polymarker CSV file"
    )
    
    parser.add_argument(
        "reference",
        type=Path,
        help="Reference genome BLAST database"
    )
    
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("output"),
        help="Output directory (default: output)"
    )
    
    parser.add_argument(
        "--no-kasp",
        action="store_true",
        help="Skip KASP primer design"
    )
    
    parser.add_argument(
        "--no-caps",
        action="store_true",
        help="Skip CAPS primer design"
    )
    
    parser.add_argument(
        "--max-price",
        type=int,
        default=200,
        help="Maximum enzyme price for CAPS (default: 200)"
    )
    
    parser.add_argument(
        "--max-tm",
        type=float,
        default=63.0,
        help="Maximum primer Tm (default: 63.0)"
    )
    
    parser.add_argument(
        "--max-size",
        type=int,
        default=25,
        help="Maximum primer size (default: 25)"
    )
    
    parser.add_argument(
        "--pick-anyway",
        action="store_true",
        help="Pick primers even if constraints are violated"
    )
    
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads for BLAST (default: 1)"
    )
    
    parser.add_argument(
        "--show-variant-sites",
        action="store_true",
        help="Show variant sites in KASP output (default: hidden for cleaner output)"
    )
    
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level (default: INFO)"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.log_level)
    
    # Create configuration
    config = PipelineConfig(
        input_file=args.input_file,
        reference_file=args.reference,
        output_dir=args.output,
        design_kasp=not args.no_kasp,
        design_caps=not args.no_caps,
        max_price=args.max_price,
        max_tm=args.max_tm,
        max_primer_size=args.max_size,
        pick_anyway=args.pick_anyway,
        threads=args.threads,
        show_variant_sites=args.show_variant_sites
    )
    
    try:
        run_pipeline(config)
    except PipelineError as e:
        logging.error(f"Pipeline failed: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        logging.info("Pipeline interrupted by user")
        sys.exit(1)


if __name__ == "__main__":
    main()