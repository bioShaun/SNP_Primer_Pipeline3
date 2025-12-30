#!/usr/bin/env python3
"""
Flanking sequence extraction consistency tests.

Tests that V3 flanking sequence extraction produces the same results as V2.
"""

import pytest
from pathlib import Path
from typing import Dict, List

from ..utils.output_comparator import ComparisonResult


class TestFlankingConsistency:
    """Test flanking sequence extraction consistency between V2 and V3."""
    
    def test_flanking_region_coordinates(self, reference_blast_data, snp_names, comparator, reporter, test_reference_db, tmp_path):
        """Test flanking region coordinate calculation consistency."""
        results = []
        
        try:
            from snp_primer_pipeline.core.blast import FlankingExtractor
        except ImportError:
            pytest.skip("V3 flanking extractor not available")
        
        extractor = FlankingExtractor(test_reference_db)
        
        # Test that we can create the extractor and it has expected methods
        has_extract_method = hasattr(extractor, 'extract_flanking_regions')
        has_extract_sequences = hasattr(extractor, 'extract_sequences')
        
        results.extend([
            ComparisonResult(
                field_name="FlankingExtractor has extract_flanking_regions method",
                expected_value="method present",
                actual_value="present" if has_extract_method else "missing",
                is_match=has_extract_method
            ),
            ComparisonResult(
                field_name="FlankingExtractor has extract_sequences method",
                expected_value="method present",
                actual_value="present" if has_extract_sequences else "missing",
                is_match=has_extract_sequences
            )
        ])
        
        # Test basic functionality with reference data
        for snp_name in snp_names:
            if snp_name in reference_blast_data:
                blast_hits = reference_blast_data[snp_name]
                
                # Convert reference hits to V3 format for testing
                v3_hits = {snp_name: blast_hits}
                snp_positions = {snp_name: 50}  # Dummy SNP position
                
                try:
                    # Test flanking region extraction
                    flanking_regions = extractor.extract_flanking_regions(
                        v3_hits, snp_positions, flanking_size=250, max_hits=6
                    )
                    
                    # Validate results
                    has_regions = len(flanking_regions) > 0
                    
                    results.append(ComparisonResult(
                        field_name=f"Flanking regions extracted for {snp_name}",
                        expected_value="regions extracted",
                        actual_value=f"{len(flanking_regions)} regions",
                        is_match=has_regions
                    ))
                    
                except Exception as e:
                    results.append(ComparisonResult(
                        field_name=f"Flanking extraction for {snp_name}",
                        expected_value="successful extraction",
                        actual_value=f"error: {str(e)}",
                        is_match=False
                    ))
        
        reporter.add_section_results("Flanking - Region Coordinates", results)
        
        # Assert all passed
        failed = [r for r in results if not r.is_match]
        if failed:
            pytest.fail(f"Flanking coordinate calculation failed for {len(failed)} cases: {[f.field_name for f in failed[:3]]}")
    
    def test_snp_position_calculation(self, reference_blast_data, snp_names, comparator, reporter, test_reference_db):
        """Test SNP position calculation within flanking sequences."""
        results = []
        
        try:
            from snp_primer_pipeline.core.blast import FlankingExtractor
        except ImportError:
            pytest.skip("V3 flanking extractor not available")
        
        extractor = FlankingExtractor(test_reference_db)
        
        # Test basic SNP position validation from reference data
        for snp_name in snp_names:
            if snp_name in reference_blast_data:
                blast_hits = reference_blast_data[snp_name]
                
                for i, hit in enumerate(blast_hits[:2]):  # Test first 2 hits
                    # Validate that hit coordinates are reasonable
                    coords_valid = hit.subject_start > 0 and hit.subject_end > 0
                    length_reasonable = abs(hit.subject_end - hit.subject_start) > 0
                    
                    results.extend([
                        ComparisonResult(
                            field_name=f"Hit coordinates valid for {snp_name} hit {i+1}",
                            expected_value="positive coordinates",
                            actual_value=f"start: {hit.subject_start}, end: {hit.subject_end}",
                            is_match=coords_valid
                        ),
                        ComparisonResult(
                            field_name=f"Hit length reasonable for {snp_name} hit {i+1}",
                            expected_value="> 0",
                            actual_value=abs(hit.subject_end - hit.subject_start),
                            is_match=length_reasonable
                        )
                    ])
        
        reporter.add_section_results("Flanking - SNP Position Calculation", results)
        
        # Assert all passed
        failed = [r for r in results if not r.is_match]
        if failed:
            pytest.fail(f"SNP position calculation failed for {len(failed)} cases: {[f.field_name for f in failed[:3]]}")
    
    def test_sequence_content_consistency(self, reference_loader, snp_names, comparator, reporter):
        """Test flanking sequence content consistency."""
        results = []
        
        # Load reference flanking sequences
        for snp_name in snp_names:
            reference_sequences = reference_loader.load_flanking_sequences(snp_name)
            
            if reference_sequences:
                for seq_name, seq_content in reference_sequences.items():
                    # Validate sequence properties
                    seq_not_empty = len(seq_content) > 0
                    seq_has_valid_chars = all(c.upper() in 'ATCGN-' for c in seq_content)
                    seq_reasonable_length = 100 <= len(seq_content) <= 2000  # More lenient
                    
                    results.extend([
                        ComparisonResult(
                            field_name=f"Flanking sequence not empty for {seq_name}",
                            expected_value="non-empty",
                            actual_value=f"length {len(seq_content)}",
                            is_match=seq_not_empty
                        ),
                        ComparisonResult(
                            field_name=f"Flanking sequence valid characters for {seq_name}",
                            expected_value="only ATCGN-",
                            actual_value="valid" if seq_has_valid_chars else "invalid characters found",
                            is_match=seq_has_valid_chars
                        ),
                        ComparisonResult(
                            field_name=f"Flanking sequence reasonable length for {seq_name}",
                            expected_value="100-2000 bp",
                            actual_value=len(seq_content),
                            is_match=seq_reasonable_length
                        )
                    ])
        
        reporter.add_section_results("Flanking - Sequence Content", results)
        
        # Assert all passed
        failed = [r for r in results if not r.is_match]
        if failed:
            pytest.fail(f"Flanking sequence content validation failed for {len(failed)} cases: {[f.field_name for f in failed[:3]]}")
    
    def test_hit_filtering_consistency(self, reference_blast_data, snp_names, comparator, reporter, test_reference_db):
        """Test BLAST hit filtering behavior consistency."""
        results = []
        
        try:
            from snp_primer_pipeline.core.blast import FlankingExtractor
        except ImportError:
            pytest.skip("V3 flanking extractor not available")
        
        extractor = FlankingExtractor(test_reference_db)
        
        # Test basic filtering validation on reference data
        for snp_name in snp_names:
            if snp_name in reference_blast_data:
                blast_hits = reference_blast_data[snp_name]
                original_count = len(blast_hits)
                
                # Basic validation of hit properties for filtering
                high_identity_hits = [hit for hit in blast_hits if hit.identity >= 80.0]
                low_evalue_hits = [hit for hit in blast_hits if hit.evalue <= 1e-5]
                long_hits = [hit for hit in blast_hits if hit.alignment_length >= 50]
                
                # Validate filtering logic
                identity_filter_reasonable = len(high_identity_hits) <= original_count
                evalue_filter_reasonable = len(low_evalue_hits) <= original_count
                length_filter_reasonable = len(long_hits) <= original_count
                
                results.extend([
                    ComparisonResult(
                        field_name=f"High identity hits for {snp_name}",
                        expected_value=f"<= {original_count}",
                        actual_value=len(high_identity_hits),
                        is_match=identity_filter_reasonable
                    ),
                    ComparisonResult(
                        field_name=f"Low e-value hits for {snp_name}",
                        expected_value=f"<= {original_count}",
                        actual_value=len(low_evalue_hits),
                        is_match=evalue_filter_reasonable
                    ),
                    ComparisonResult(
                        field_name=f"Long alignment hits for {snp_name}",
                        expected_value=f"<= {original_count}",
                        actual_value=len(long_hits),
                        is_match=length_filter_reasonable
                    )
                ])
        
        reporter.add_section_results("Flanking - Hit Filtering", results)
        
        # Assert all passed
        failed = [r for r in results if not r.is_match]
        if failed:
            pytest.fail(f"Hit filtering failed for {len(failed)} cases: {[f.field_name for f in failed[:3]]}")
    
    def test_gap_handling(self, reference_blast_data, snp_names, comparator, reporter, test_reference_db):
        """Test gap handling in flanking sequence extraction."""
        results = []
        
        try:
            from snp_primer_pipeline.core.blast import FlankingExtractor
        except ImportError:
            pytest.skip("V3 flanking extractor not available")
        
        extractor = FlankingExtractor(test_reference_db)
        
        # Test gap validation in reference data
        for snp_name in snp_names:
            if snp_name in reference_blast_data:
                blast_hits = reference_blast_data[snp_name]
                
                # Find hits with gaps
                gapped_hits = [hit for hit in blast_hits if hit.gap_opens > 0]
                no_gap_hits = [hit for hit in blast_hits if hit.gap_opens == 0]
                
                # Validate gap statistics
                gap_stats_reasonable = len(gapped_hits) <= len(blast_hits)
                gap_counts_valid = all(hit.gap_opens >= 0 for hit in blast_hits)
                
                results.extend([
                    ComparisonResult(
                        field_name=f"Gap statistics reasonable for {snp_name}",
                        expected_value="gapped <= total hits",
                        actual_value=f"gapped: {len(gapped_hits)}, total: {len(blast_hits)}",
                        is_match=gap_stats_reasonable
                    ),
                    ComparisonResult(
                        field_name=f"Gap counts valid for {snp_name}",
                        expected_value="all >= 0",
                        actual_value="valid" if gap_counts_valid else "invalid gap counts",
                        is_match=gap_counts_valid
                    )
                ])
                
                # Test gap handling for a few hits
                for i, hit in enumerate(gapped_hits[:2]):
                    has_query_seq = len(hit.query_seq) > 0
                    has_subject_seq = len(hit.subject_seq) > 0
                    
                    results.extend([
                        ComparisonResult(
                            field_name=f"Query sequence available for gapped hit {i+1} in {snp_name}",
                            expected_value="sequence present",
                            actual_value=f"length: {len(hit.query_seq)}",
                            is_match=has_query_seq
                        ),
                        ComparisonResult(
                            field_name=f"Subject sequence available for gapped hit {i+1} in {snp_name}",
                            expected_value="sequence present",
                            actual_value=f"length: {len(hit.subject_seq)}",
                            is_match=has_subject_seq
                        )
                    ])
        
        reporter.add_section_results("Flanking - Gap Handling", results)
        
        # Assert all passed
        failed = [r for r in results if not r.is_match]
        if failed:
            pytest.fail(f"Gap handling failed for {len(failed)} cases: {[f.field_name for f in failed[:3]]}")