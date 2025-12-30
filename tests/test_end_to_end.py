#!/usr/bin/env python3
"""
End-to-end integration tests for SNP Primer Pipeline.

Tests the complete pipeline using example data.
"""

import pytest
import tempfile
import shutil
from pathlib import Path
import logging

from snp_primer_pipeline.config import PipelineConfig
from snp_primer_pipeline.main import run_pipeline, setup_logging


class TestEndToEnd:
    """End-to-end integration tests."""
    
    @pytest.fixture(autouse=True)
    def setup_logging_fixture(self):
        """Setup logging for tests."""
        setup_logging("WARNING")  # Reduce log noise in tests
    
    @pytest.fixture
    def test_data_dir(self):
        """Get test data directory."""
        return Path(__file__).parent.parent / "test_data"
    
    @pytest.fixture
    def input_file(self, test_data_dir):
        """Get input polymarker file."""
        return test_data_dir / "polymarker_input_example.csv"
    
    @pytest.fixture
    def reference_db(self, test_data_dir):
        """Get reference BLAST database."""
        return test_data_dir / "blastdb" / "test_reference.fa"
    
    @pytest.fixture
    def temp_output_dir(self):
        """Create temporary output directory."""
        with tempfile.TemporaryDirectory() as temp_dir:
            yield Path(temp_dir)
    
    def test_pipeline_basic_functionality(self, input_file, reference_db, temp_output_dir):
        """Test basic pipeline functionality with example data."""
        
        # Skip if input files don't exist
        if not input_file.exists():
            pytest.skip(f"Input file not found: {input_file}")
        if not reference_db.exists():
            pytest.skip(f"Reference database not found: {reference_db}")
        
        # Create configuration
        config = PipelineConfig(
            input_file=input_file,
            reference_file=reference_db,
            output_dir=temp_output_dir,
            design_kasp=True,
            design_caps=False,  # Skip CAPS to avoid potential issues
            max_price=200,
            max_tm=63.0,
            max_primer_size=25,
            pick_anyway=True,
            threads=1
        )
        
        # Run pipeline
        run_pipeline(config)
        
        # Verify expected output files exist
        expected_files = [
            "for_blast.fa",
            "blast_out.txt", 
            "flanking_sequences.fa"
        ]
        
        for expected_file in expected_files:
            file_path = temp_output_dir / expected_file
            assert file_path.exists(), f"Expected output file not found: {expected_file}"
            assert file_path.stat().st_size > 0, f"Output file is empty: {expected_file}"
    
    def test_input_parsing(self, input_file, reference_db, temp_output_dir):
        """Test that input parsing works correctly."""
        
        if not input_file.exists():
            pytest.skip(f"Input file not found: {input_file}")
        if not reference_db.exists():
            pytest.skip(f"Reference database not found: {reference_db}")
        
        config = PipelineConfig(
            input_file=input_file,
            reference_file=reference_db,
            output_dir=temp_output_dir,
            design_kasp=True,
            design_caps=False,
            pick_anyway=True,
            threads=1
        )
        
        run_pipeline(config)
        
        # Check FASTA conversion
        fasta_file = temp_output_dir / "for_blast.fa"
        assert fasta_file.exists()
        
        with open(fasta_file, 'r') as f:
            content = f.read()
            
        # Should contain 2 sequences (from example file)
        assert content.count('>') == 2, "Expected 2 sequences in FASTA output"
        
        # Should contain IUPAC codes converted to degenerate bases
        assert 'Y' in content or 'R' in content, "Expected IUPAC codes in FASTA output"
    
    def test_blast_execution(self, input_file, reference_db, temp_output_dir):
        """Test that BLAST execution produces results."""
        
        if not input_file.exists():
            pytest.skip(f"Input file not found: {input_file}")
        if not reference_db.exists():
            pytest.skip(f"Reference database not found: {reference_db}")
        
        config = PipelineConfig(
            input_file=input_file,
            reference_file=reference_db,
            output_dir=temp_output_dir,
            design_kasp=True,
            design_caps=False,
            pick_anyway=True,
            threads=1
        )
        
        run_pipeline(config)
        
        # Check BLAST output
        blast_file = temp_output_dir / "blast_out.txt"
        assert blast_file.exists()
        
        with open(blast_file, 'r') as f:
            lines = f.readlines()
        
        # Should have BLAST hits
        assert len(lines) > 0, "Expected BLAST results"
        
        # Check format (tab-delimited)
        for line in lines:
            if line.strip():
                fields = line.strip().split('\t')
                assert len(fields) >= 12, "Expected at least 12 fields in BLAST output"
    
    def test_flanking_sequence_extraction(self, input_file, reference_db, temp_output_dir):
        """Test that flanking sequences are extracted correctly."""
        
        if not input_file.exists():
            pytest.skip(f"Input file not found: {input_file}")
        if not reference_db.exists():
            pytest.skip(f"Reference database not found: {reference_db}")
        
        config = PipelineConfig(
            input_file=input_file,
            reference_file=reference_db,
            output_dir=temp_output_dir,
            design_kasp=True,
            design_caps=False,
            pick_anyway=True,
            threads=1
        )
        
        run_pipeline(config)
        
        # Check flanking sequences
        flanking_file = temp_output_dir / "flanking_sequences.fa"
        assert flanking_file.exists()
        
        with open(flanking_file, 'r') as f:
            content = f.read()
        
        # Should contain extracted sequences
        seq_count = content.count('>')
        assert seq_count > 0, "Expected flanking sequences"
        
        # Sequences should be reasonably long (flanking regions)
        sequences = content.split('>')[1:]  # Skip first empty element
        for seq in sequences:
            lines = seq.strip().split('\n')
            if len(lines) > 1:
                seq_data = ''.join(lines[1:])
                assert len(seq_data) > 100, "Expected flanking sequences to be substantial length"
    
    def test_snp_processing_directories(self, input_file, reference_db, temp_output_dir):
        """Test that SNP-specific directories are created."""
        
        if not input_file.exists():
            pytest.skip(f"Input file not found: {input_file}")
        if not reference_db.exists():
            pytest.skip(f"Reference database not found: {reference_db}")
        
        config = PipelineConfig(
            input_file=input_file,
            reference_file=reference_db,
            output_dir=temp_output_dir,
            design_kasp=True,
            design_caps=False,
            pick_anyway=True,
            threads=1
        )
        
        run_pipeline(config)
        
        # Check for SNP directories
        snp_dirs = list(temp_output_dir.glob("SNP_*"))
        assert len(snp_dirs) > 0, "Expected SNP-specific directories"
        
        # Check directory contents
        for snp_dir in snp_dirs:
            assert snp_dir.is_dir(), f"Expected {snp_dir} to be a directory"
            
            # Should contain flanking sequence file
            flanking_files = list(snp_dir.glob("flanking_*.fa"))
            assert len(flanking_files) > 0, f"Expected flanking file in {snp_dir}"
    
    def test_kasp_primer_design(self, input_file, reference_db, temp_output_dir):
        """Test KASP primer design functionality."""
        
        if not input_file.exists():
            pytest.skip(f"Input file not found: {input_file}")
        if not reference_db.exists():
            pytest.skip(f"Reference database not found: {reference_db}")
        
        config = PipelineConfig(
            input_file=input_file,
            reference_file=reference_db,
            output_dir=temp_output_dir,
            design_kasp=True,
            design_caps=False,
            pick_anyway=True,
            threads=1
        )
        
        run_pipeline(config)
        
        # Check for KASP output files
        snp_dirs = list(temp_output_dir.glob("SNP_*"))
        
        kasp_files_found = False
        for snp_dir in snp_dirs:
            kasp_files = list(snp_dir.glob("KASP_primers_*.txt"))
            if kasp_files:
                kasp_files_found = True
                
                # Check file content
                for kasp_file in kasp_files:
                    assert kasp_file.stat().st_size > 0, f"KASP file should not be empty: {kasp_file}"
        
        # At least one KASP file should be created (even if empty due to no primers found)
        assert kasp_files_found, "Expected at least one KASP primer file"
    
    def test_pipeline_with_missing_input(self, reference_db, temp_output_dir):
        """Test pipeline behavior with missing input file."""
        
        if not reference_db.exists():
            pytest.skip(f"Reference database not found: {reference_db}")
        
        missing_input = Path("nonexistent_file.csv")
        
        # Should raise an exception during config creation for missing input
        from snp_primer_pipeline.exceptions import ConfigurationError
        with pytest.raises(ConfigurationError):
            config = PipelineConfig(
                input_file=missing_input,
                reference_file=reference_db,
                output_dir=temp_output_dir,
                design_kasp=True,
                design_caps=False,
                pick_anyway=True,
                threads=1
            )
    
    def test_pipeline_with_missing_reference(self, input_file, temp_output_dir):
        """Test pipeline behavior with missing reference database."""
        
        if not input_file.exists():
            pytest.skip(f"Input file not found: {input_file}")
        
        missing_ref = Path("nonexistent_reference.fa")
        
        config = PipelineConfig(
            input_file=input_file,
            reference_file=missing_ref,
            output_dir=temp_output_dir,
            design_kasp=True,
            design_caps=False,
            pick_anyway=True,
            threads=1
        )
        
        # Should raise an exception during pipeline execution for missing reference
        with pytest.raises(Exception):
            run_pipeline(config)


if __name__ == "__main__":
    # Run tests directly
    pytest.main([__file__, "-v"])