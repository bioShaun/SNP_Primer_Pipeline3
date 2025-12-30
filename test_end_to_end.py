#!/usr/bin/env python3
"""
End-to-end test for SNP Primer Pipeline using example data from SNP_Primer_Pipeline2.
"""

import sys
import tempfile
from pathlib import Path
import logging

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

from snp_primer_pipeline.config import PipelineConfig
from snp_primer_pipeline.main import run_pipeline, setup_logging


def test_end_to_end():
    """Test the complete pipeline with example data."""
    
    # Setup logging
    setup_logging("INFO")
    logger = logging.getLogger(__name__)
    
    # Input files
    input_file = Path("test_data/polymarker_input_example.csv")
    reference_db = Path("test_data/blastdb/test_reference.fa")
    
    # Check input files exist
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        return False
        
    if not reference_db.exists():
        logger.error(f"Reference database not found: {reference_db}")
        return False
    
    # Create temporary output directory
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir) / "pipeline_output"
        
        # Create configuration
        config = PipelineConfig(
            input_file=input_file,
            reference_file=reference_db,
            output_dir=output_dir,
            design_kasp=True,
            design_caps=True,
            max_price=200,
            max_tm=63.0,
            max_primer_size=25,
            pick_anyway=True,  # Allow primers even if constraints violated
            threads=1
        )
        
        try:
            logger.info("Starting end-to-end pipeline test...")
            run_pipeline(config)
            
            # Check if output files were created
            expected_files = [
                "for_blast.fa",
                "blast_out.txt",
                "flanking_sequences.fa"
            ]
            
            success = True
            for expected_file in expected_files:
                file_path = output_dir / expected_file
                if file_path.exists():
                    logger.info(f"✓ Created: {expected_file}")
                else:
                    logger.error(f"✗ Missing: {expected_file}")
                    success = False
            
            # Check for SNP-specific directories
            snp_dirs = list(output_dir.glob("SNP_*"))
            if snp_dirs:
                logger.info(f"✓ Created {len(snp_dirs)} SNP directories")
                
                # Check contents of first SNP directory
                first_snp_dir = snp_dirs[0]
                snp_files = list(first_snp_dir.glob("*"))
                logger.info(f"✓ SNP directory contains {len(snp_files)} files")
                for snp_file in snp_files:
                    logger.info(f"  - {snp_file.name}")
            else:
                logger.error("✗ No SNP directories created")
                success = False
            
            if success:
                logger.info("🎉 End-to-end test PASSED!")
                return True
            else:
                logger.error("❌ End-to-end test FAILED!")
                return False
                
        except Exception as e:
            logger.error(f"❌ Pipeline failed with error: {e}")
            import traceback
            traceback.print_exc()
            return False


if __name__ == "__main__":
    success = test_end_to_end()
    sys.exit(0 if success else 1)