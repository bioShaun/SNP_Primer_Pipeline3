#!/usr/bin/env python3
"""
Inspect the output of the pipeline to verify correctness.
"""

import sys
import tempfile
from pathlib import Path
import logging

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

from snp_primer_pipeline.config import PipelineConfig
from snp_primer_pipeline.main import run_pipeline, setup_logging


def inspect_pipeline_output():
    """Run pipeline and inspect the output files."""
    
    # Setup logging
    setup_logging("INFO")
    logger = logging.getLogger(__name__)
    
    # Input files
    input_file = Path("test_data/polymarker_input_example.csv")
    reference_db = Path("test_data/blastdb/test_reference.fa")
    
    # Create output directory (not temporary so we can inspect)
    output_dir = Path("test_output")
    output_dir.mkdir(exist_ok=True)
    
    # Create configuration
    config = PipelineConfig(
        input_file=input_file,
        reference_file=reference_db,
        output_dir=output_dir,
        design_kasp=True,
        design_caps=False,  # Skip CAPS to avoid recursion issue for now
        max_price=200,
        max_tm=63.0,
        max_primer_size=25,
        pick_anyway=True,
        threads=1
    )
    
    try:
        logger.info("Running pipeline for output inspection...")
        run_pipeline(config)
        
        # Inspect generated files
        logger.info("\n=== PIPELINE OUTPUT INSPECTION ===")
        
        # 1. Check input parsing
        logger.info("\n1. Input file content:")
        with open(input_file, 'r') as f:
            for line in f:
                logger.info(f"  {line.strip()}")
        
        # 2. Check FASTA conversion
        fasta_file = output_dir / "for_blast.fa"
        if fasta_file.exists():
            logger.info("\n2. Generated FASTA for BLAST:")
            with open(fasta_file, 'r') as f:
                content = f.read()
                logger.info(f"  Content (first 200 chars): {content[:200]}...")
        
        # 3. Check BLAST output
        blast_file = output_dir / "blast_out.txt"
        if blast_file.exists():
            logger.info("\n3. BLAST output:")
            with open(blast_file, 'r') as f:
                lines = f.readlines()
                logger.info(f"  Total lines: {len(lines)}")
                if lines:
                    logger.info(f"  First few lines:")
                    for i, line in enumerate(lines[:5]):
                        logger.info(f"    {i+1}: {line.strip()}")
        
        # 4. Check flanking sequences
        flanking_file = output_dir / "flanking_sequences.fa"
        if flanking_file.exists():
            logger.info("\n4. Flanking sequences:")
            with open(flanking_file, 'r') as f:
                content = f.read()
                seq_count = content.count('>')
                logger.info(f"  Number of sequences: {seq_count}")
                logger.info(f"  Content (first 300 chars): {content[:300]}...")
        
        # 5. Check SNP directories
        snp_dirs = list(output_dir.glob("SNP_*"))
        logger.info(f"\n5. SNP directories ({len(snp_dirs)} found):")
        
        for snp_dir in snp_dirs:
            logger.info(f"\n  Directory: {snp_dir.name}")
            files = list(snp_dir.glob("*"))
            for file in files:
                logger.info(f"    - {file.name} ({file.stat().st_size} bytes)")
                
                # Show content of small text files
                if file.suffix in ['.txt', '.fa'] and file.stat().st_size < 1000:
                    with open(file, 'r') as f:
                        content = f.read().strip()
                        logger.info(f"      Content: {content[:200]}...")
        
        logger.info("\n=== INSPECTION COMPLETE ===")
        return True
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = inspect_pipeline_output()
    sys.exit(0 if success else 1)