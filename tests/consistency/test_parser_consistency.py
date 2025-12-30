#!/usr/bin/env python3
"""
Parser consistency tests.

Tests that V3 parser produces the same output as V2 parser.
"""

import pytest
from pathlib import Path
from typing import Dict, List

from ..utils.output_comparator import ComparisonResult


class TestParserConsistency:
    """Test parser consistency between V2 and V3."""
    
    def test_snp_name_extraction(self, reference_polymarker_input, comparator, reporter, tmp_path):
        """Test SNP name extraction consistency."""
        results = []
        
        try:
            from snp_primer_pipeline.core.parser import PolymarkerParser
        except ImportError:
            pytest.skip("V3 parser not available")
        
        # Create temporary input file with test data
        test_input_file = tmp_path / "test_input.csv"
        with open(test_input_file, 'w') as f:
            for snp_data in reference_polymarker_input:
                f.write(f"{snp_data['name']},{snp_data['chromosome']},{snp_data['sequence']}\n")
        
        # Parse with V3
        parser = PolymarkerParser(test_input_file)
        parsed_snps = parser.parse()
        
        # Test each SNP from reference input
        for i, snp_data in enumerate(reference_polymarker_input):
            expected_name = snp_data['name']
            
            if i < len(parsed_snps):
                actual_name = parsed_snps[i].name
                
                results.append(comparator.compare_sequences(
                    expected=expected_name,
                    actual=actual_name,
                    name=f"SNP name extraction for {expected_name}"
                ))
            else:
                results.append(ComparisonResult(
                    field_name=f"SNP name extraction for {expected_name}",
                    expected_value=expected_name,
                    actual_value="<missing>",
                    is_match=False,
                    message="SNP not parsed"
                ))
        
        reporter.add_section_results("Parser - SNP Name Extraction", results)
        
        # Assert all passed
        failed = [r for r in results if not r.is_match]
        if failed:
            pytest.fail(f"SNP name extraction failed for {len(failed)} cases: {[f.field_name for f in failed[:3]]}")
    
    def test_chromosome_extraction(self, reference_polymarker_input, comparator, reporter, tmp_path):
        """Test chromosome information extraction consistency."""
        results = []
        
        try:
            from snp_primer_pipeline.core.parser import PolymarkerParser
        except ImportError:
            pytest.skip("V3 parser not available")
        
        # Create temporary input file with test data
        test_input_file = tmp_path / "test_input.csv"
        with open(test_input_file, 'w') as f:
            for snp_data in reference_polymarker_input:
                f.write(f"{snp_data['name']},{snp_data['chromosome']},{snp_data['sequence']}\n")
        
        # Parse with V3
        parser = PolymarkerParser(test_input_file)
        parsed_snps = parser.parse()
        
        # Test each SNP from reference input
        for i, snp_data in enumerate(reference_polymarker_input):
            expected_chr = snp_data['chromosome']
            
            if i < len(parsed_snps):
                actual_chr = parsed_snps[i].chromosome
                
                results.append(comparator.compare_sequences(
                    expected=expected_chr,
                    actual=actual_chr,
                    name=f"Chromosome extraction for {snp_data['name']}"
                ))
            else:
                results.append(ComparisonResult(
                    field_name=f"Chromosome extraction for {snp_data['name']}",
                    expected_value=expected_chr,
                    actual_value="<missing>",
                    is_match=False,
                    message="SNP not parsed"
                ))
        
        reporter.add_section_results("Parser - Chromosome Extraction", results)
        
        # Assert all passed
        failed = [r for r in results if not r.is_match]
        if failed:
            pytest.fail(f"Chromosome extraction failed for {len(failed)} cases: {[f.field_name for f in failed[:3]]}")
    
    def test_iupac_conversion(self, reference_polymarker_input, comparator, reporter, tmp_path):
        """Test IUPAC code conversion consistency."""
        results = []
        
        try:
            from snp_primer_pipeline.core.parser import PolymarkerParser
        except ImportError:
            pytest.skip("V3 parser not available")
        
        # Expected IUPAC conversions based on V2 behavior
        expected_conversions = {
            '[T/C]': 'Y',
            '[A/G]': 'R',
            '[A/T]': 'W',
            '[C/G]': 'S',
            '[A/C]': 'M',
            '[G/T]': 'K'
        }
        
        # Test each conversion
        for snp_notation, expected_iupac in expected_conversions.items():
            test_sequence = f"ATCGATCGATCG{snp_notation}GATCGATCGATC"
            
            # Create temporary input file
            test_input_file = tmp_path / f"test_iupac_{snp_notation.replace('/', '_').replace('[', '').replace(']', '')}.csv"
            with open(test_input_file, 'w') as f:
                f.write(f"test_snp,chr1,{test_sequence}\n")
            
            # Parse with V3
            parser = PolymarkerParser(test_input_file)
            parsed_snps = parser.parse()
            
            if parsed_snps:
                actual_iupac = parsed_snps[0].iupac_code
                
                results.append(comparator.compare_sequences(
                    expected=expected_iupac,
                    actual=actual_iupac,
                    name=f"IUPAC conversion for {snp_notation}"
                ))
            else:
                results.append(ComparisonResult(
                    field_name=f"IUPAC conversion for {snp_notation}",
                    expected_value=expected_iupac,
                    actual_value="<no parse result>",
                    is_match=False,
                    message="Failed to parse SNP"
                ))
        
        reporter.add_section_results("Parser - IUPAC Conversion", results)
        
        # Assert all passed
        failed = [r for r in results if not r.is_match]
        if failed:
            pytest.fail(f"IUPAC conversion failed for {len(failed)} cases: {[f.field_name for f in failed[:3]]}")
    
    def test_fasta_output_format(self, reference_polymarker_input, comparator, reporter, tmp_path):
        """Test FASTA output format consistency."""
        results = []
        
        try:
            from snp_primer_pipeline.core.parser import PolymarkerParser
        except ImportError:
            pytest.skip("V3 parser not available")
        
        # Create temporary input file with test data
        test_input_file = tmp_path / "test_input.csv"
        with open(test_input_file, 'w') as f:
            for snp_data in reference_polymarker_input:
                f.write(f"{snp_data['name']},{snp_data['chromosome']},{snp_data['sequence']}\n")
        
        # Parse with V3
        parser = PolymarkerParser(test_input_file)
        parsed_snps = parser.parse()
        
        # Generate FASTA output
        fasta_output_file = tmp_path / "test_output.fa"
        parser.to_fasta(fasta_output_file)
        
        # Read and check FASTA output
        with open(fasta_output_file, 'r') as f:
            fasta_content = f.read()
        
        lines = fasta_content.strip().split('\n')
        
        # Test FASTA format for each SNP
        for i, snp_data in enumerate(reference_polymarker_input):
            if i < len(parsed_snps):
                # Each SNP should have 2 lines (header + sequence)
                line_start = i * 2
                
                if line_start + 1 < len(lines):
                    header_line = lines[line_start]
                    sequence_line = lines[line_start + 1]
                    
                    # Header should start with >
                    header_ok = header_line.startswith('>')
                    results.append(ComparisonResult(
                        field_name=f"FASTA header format for {snp_data['name']}",
                        expected_value="starts with >",
                        actual_value=f"starts with {header_line[:1] if header_line else 'empty'}",
                        is_match=header_ok
                    ))
                    
                    # Header should contain SNP name
                    name_in_header = snp_data['name'] in header_line
                    results.append(ComparisonResult(
                        field_name=f"SNP name in FASTA header for {snp_data['name']}",
                        expected_value=f"contains {snp_data['name']}",
                        actual_value=header_line,
                        is_match=name_in_header
                    ))
                    
                    # Sequence should not be empty
                    seq_not_empty = len(sequence_line.strip()) > 0
                    results.append(ComparisonResult(
                        field_name=f"FASTA sequence not empty for {snp_data['name']}",
                        expected_value="non-empty sequence",
                        actual_value=f"length {len(sequence_line.strip())}",
                        is_match=seq_not_empty
                    ))
        
        reporter.add_section_results("Parser - FASTA Output Format", results)
        
        # Assert all passed
        failed = [r for r in results if not r.is_match]
        if failed:
            pytest.fail(f"FASTA output format failed for {len(failed)} cases: {[f.field_name for f in failed[:3]]}")
    
    def test_sequence_processing(self, reference_polymarker_input, comparator, reporter, tmp_path):
        """Test sequence processing consistency."""
        results = []
        
        try:
            from snp_primer_pipeline.core.parser import PolymarkerParser
        except ImportError:
            pytest.skip("V3 parser not available")
        
        # Create temporary input file with test data
        test_input_file = tmp_path / "test_input.csv"
        with open(test_input_file, 'w') as f:
            for snp_data in reference_polymarker_input:
                f.write(f"{snp_data['name']},{snp_data['chromosome']},{snp_data['sequence']}\n")
        
        # Parse with V3
        parser = PolymarkerParser(test_input_file)
        parsed_snps = parser.parse()
        
        # Test sequence processing for each SNP
        for i, snp_data in enumerate(reference_polymarker_input):
            if i < len(parsed_snps):
                original_seq = snp_data['sequence']
                processed_seq = parsed_snps[i].flanking_sequence
                
                # Check that brackets are removed/converted
                has_brackets = '[' in processed_seq or ']' in processed_seq
                results.append(ComparisonResult(
                    field_name=f"Brackets removed from sequence for {snp_data['name']}",
                    expected_value="no brackets",
                    actual_value=f"brackets present: {has_brackets}",
                    is_match=not has_brackets
                ))
                
                # Check sequence length is reasonable (should be similar to original minus brackets)
                original_clean = original_seq.replace('[', '').replace(']', '').replace('/', '')
                length_diff = abs(len(processed_seq) - len(original_clean))
                length_ok = length_diff <= 2  # Allow small differences for IUPAC conversion
                
                results.append(ComparisonResult(
                    field_name=f"Processed sequence length for {snp_data['name']}",
                    expected_value=f"~{len(original_clean)} (±2)",
                    actual_value=len(processed_seq),
                    is_match=length_ok
                ))
        
        reporter.add_section_results("Parser - Sequence Processing", results)
        
        # Assert all passed
        failed = [r for r in results if not r.is_match]
        if failed:
            pytest.fail(f"Sequence processing failed for {len(failed)} cases: {[f.field_name for f in failed[:3]]}")
    
    def test_error_handling(self, comparator, reporter, tmp_path):
        """Test parser error handling consistency."""
        results = []
        
        try:
            from snp_primer_pipeline.core.parser import PolymarkerParser
        except ImportError:
            pytest.skip("V3 parser not available")
        
        # Test invalid input formats
        invalid_inputs = [
            "",  # Empty string
            "invalid_format",  # No commas
            "name,chr",  # Missing sequence
            "name,chr,",  # Empty sequence
            "name,chr,ATCG[X/Y]ATCG",  # Invalid IUPAC
        ]
        
        for i, invalid_input in enumerate(invalid_inputs):
            # Create temporary input file with invalid data
            test_input_file = tmp_path / f"invalid_test_{i}.csv"
            with open(test_input_file, 'w') as f:
                if invalid_input:  # Don't write empty string
                    f.write(invalid_input + "\n")
            
            try:
                parser = PolymarkerParser(test_input_file)
                parsed_snps = parser.parse()
                # Should either return empty list or raise exception
                error_handled = len(parsed_snps) == 0
            except Exception:
                # Exception is acceptable error handling
                error_handled = True
            
            results.append(ComparisonResult(
                field_name=f"Error handling for '{invalid_input[:20]}...'",
                expected_value="error handled (empty result or exception)",
                actual_value="error handled" if error_handled else "no error handling",
                is_match=error_handled
            ))
        
        reporter.add_section_results("Parser - Error Handling", results)
        
        # Assert all passed
        failed = [r for r in results if not r.is_match]
        if failed:
            pytest.fail(f"Error handling failed for {len(failed)} cases: {[f.field_name for f in failed[:3]]}")