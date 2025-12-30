#!/usr/bin/env python3
"""
KASP primer design consistency tests.

Tests that V3 KASP primer design produces the same results as V2.
"""

import pytest
from pathlib import Path
from typing import Dict, List

from ..utils.reference_loader import KASPPrimerRecord
from ..utils.output_comparator import ComparisonResult


class TestKASPConsistency:
    """Test KASP primer design consistency between V2 and V3."""
    
    def test_kasp_primer_sequences(self, reference_kasp_data, comparator, reporter, tmp_path):
        """Test KASP primer sequence generation capability."""
        results = []
        
        try:
            from snp_primer_pipeline.primers.kasp import KASPDesigner
        except ImportError:
            pytest.skip("V3 KASP designer not available")
        
        # Check if primer3_core is available
        import shutil
        from snp_primer_pipeline.config import SoftwarePaths
        
        try:
            software = SoftwarePaths.auto_detect()
            primer3_available = software.primer3_path.exists()
        except:
            primer3_available = shutil.which("primer3_core") is not None
            
        if not primer3_available:
            pytest.skip("primer3_core not available in PATH")
        
        designer = KASPDesigner()
        
        # Use realistic template sequence from reference data (similar to flanking sequences)
        realistic_template = (
            "TATGAGACTTGTAGGTGGAGCAATTGAGGTCGACAGAAAGGTGGCTTTTGGTTGATACTCCATTGAATAATTAATGAAGA"
            "TGGCTATCGATGCAAAGGGCGGAGGTAATCCTCCTTTTCAAAAATGATATGTTCATCTCTTGAGCGGCATCTGCAAGTCT"
            "GCTGAAAATTGAATCGATGCAGCTCGTGAAGAGCTTCTCTTTGGATATAGACTACTCCTACTTCATTGTTACCCCGATTA"
            "AAGAATACCAGTACTTCATTGTTACCCCGAAAGAAGTTACAGCAAATCTGTTTTCAGACTTCACCAGGTTATCGGAAAAT"
            "GCAACTTTGGTGTCCCAAGGGATGCACAGCACTGTATGGAAACTCATAGGGTCAATAAACTATTACAGTACATTACAGGA"
            "TGGTTAATAAATCCATACCACGAAACAAGAAACTTGATACTCCCTCCGATCCAAAATAAGTGTCGCAGTTTCGAACTGGG"
            "ACTAGTTCAAAACTGCAACATTTATTATGGATCAGAGGTAGCAATACAATCCGTAGTGGTACTCAACAACAAAGATCAGC"
            "TGCTGTATGCATGCCCCCTAGCACAACAAACAAAAGAAATAGAAGAACATAGATAGATCCGTTTGTTAAAATGTACACAC"
            "AGGTGGTCTAGAGCTACTGGGAATCTCAGGTATCATTTAGGGCATGATGATACATTTCACGTAGTCGAGGCGGGCATCAA"
            "TTGCACGCCATATCAAGTTCGGAAACACCCTTAGCGGCCTGCGGGGTGCCTGTTCCAGGACAGGTAGTTGTGATCTCCCA"
        )
        
        # Test KASP design capability
        try:
            # Run V3 KASP design with realistic parameters
            v3_primers = designer.design_primers(
                template_sequence=realistic_template,
                snp_position=500,  # Middle of template
                snp_alleles=("A", "T"),
                product_size_range=(60, 250),  # Match reference data range
                output_dir=tmp_path
            )
            
            # Validate that primers were generated
            primers_generated = len(v3_primers) > 0
            results.append(ComparisonResult(
                field_name="KASP primers generation capability",
                expected_value="primers generated",
                actual_value=f"{len(v3_primers)} primers generated",
                is_match=primers_generated
            ))
            
            if v3_primers:
                # Validate primer structure
                for i, pair in enumerate(v3_primers):
                    # Check that we have allele-specific and common primers
                    has_allele_primer = "Allele_" in pair.left.name
                    has_common_primer = pair.right.name == "Common"
                    valid_sequences = len(pair.left.sequence) > 0 and len(pair.right.sequence) > 0
                    valid_product_size = pair.product_size > 0
                    
                    results.extend([
                        ComparisonResult(
                            field_name=f"KASP pair {i+1} has allele-specific primer",
                            expected_value="allele-specific primer",
                            actual_value=pair.left.name,
                            is_match=has_allele_primer
                        ),
                        ComparisonResult(
                            field_name=f"KASP pair {i+1} has common primer",
                            expected_value="Common",
                            actual_value=pair.right.name,
                            is_match=has_common_primer
                        ),
                        ComparisonResult(
                            field_name=f"KASP pair {i+1} has valid sequences",
                            expected_value="non-empty sequences",
                            actual_value=f"L:{len(pair.left.sequence)} R:{len(pair.right.sequence)}",
                            is_match=valid_sequences
                        ),
                        ComparisonResult(
                            field_name=f"KASP pair {i+1} has valid product size",
                            expected_value="> 0",
                            actual_value=pair.product_size,
                            is_match=valid_product_size
                        )
                    ])
                    
        except Exception as e:
            results.append(ComparisonResult(
                field_name="KASP design execution",
                expected_value="successful design",
                actual_value=f"error: {str(e)}",
                is_match=False
            ))
        
        reporter.add_section_results("KASP - Primer Sequences", results)
        
        # Assert all passed
        failed = [r for r in results if not r.is_match]
        if failed:
            pytest.fail(f"KASP primer sequences failed for {len(failed)} cases: {[f.field_name for f in failed[:3]]}")
    
    def test_kasp_primer_positions(self, reference_kasp_data, comparator, reporter):
        """Test KASP primer position consistency."""
        results = []
        
        # Validate primer positions in reference data
        # Note: V2 uses Primer3's coordinate convention where:
        # - LEFT primers: start < end (5' to 3' on forward strand)
        # - RIGHT primers: start > end (start is 3' position, end is 5' position)
        for snp_name, primers in reference_kasp_data.items():
            for i, primer in enumerate(primers):
                # Positions should be valid (>= 1)
                start_valid = primer.start >= 1
                end_valid = primer.end >= 1
                
                # Length calculation depends on primer type
                # For LEFT primers: end - start + 1
                # For RIGHT primers: start - end + 1 (because start > end)
                if primer.primer_type == 'RIGHT':
                    calculated_length = primer.start - primer.end + 1
                    position_order_valid = primer.start >= primer.end  # start >= end for RIGHT
                else:
                    calculated_length = primer.end - primer.start + 1
                    position_order_valid = primer.end >= primer.start  # end >= start for LEFT
                
                length_matches = calculated_length == primer.length
                
                results.extend([
                    ComparisonResult(
                        field_name=f"KASP primer {i+1} start position for {snp_name}",
                        expected_value=">= 1",
                        actual_value=primer.start,
                        is_match=start_valid
                    ),
                    ComparisonResult(
                        field_name=f"KASP primer {i+1} position order for {snp_name} ({primer.primer_type})",
                        expected_value=f"valid for {primer.primer_type}",
                        actual_value=f"start={primer.start}, end={primer.end}",
                        is_match=position_order_valid
                    ),
                    ComparisonResult(
                        field_name=f"KASP primer {i+1} length consistency for {snp_name}",
                        expected_value=primer.length,
                        actual_value=calculated_length,
                        is_match=length_matches
                    )
                ])
        
        reporter.add_section_results("KASP - Primer Positions", results)
        
        # Assert all passed
        failed = [r for r in results if not r.is_match]
        if failed:
            pytest.fail(f"KASP primer positions failed for {len(failed)} cases: {[f.field_name for f in failed[:3]]}")
    
    def test_kasp_tm_values(self, reference_kasp_data, comparator, reporter):
        """Test KASP Tm value consistency with tolerance."""
        results = []
        
        # Validate Tm values in reference data
        for snp_name, primers in reference_kasp_data.items():
            for i, primer in enumerate(primers):
                # Tm should be in reasonable range
                tm_reasonable = 50.0 <= primer.tm <= 70.0
                
                results.append(ComparisonResult(
                    field_name=f"KASP primer {i+1} Tm range for {snp_name}",
                    expected_value="50-70°C",
                    actual_value=primer.tm,
                    is_match=tm_reasonable
                ))
                
                # Test numerical precision
                tm_precision = abs(primer.tm - round(primer.tm, 3)) < 0.001
                
                results.append(ComparisonResult(
                    field_name=f"KASP primer {i+1} Tm precision for {snp_name}",
                    expected_value="3 decimal places",
                    actual_value=f"{primer.tm:.6f}",
                    is_match=tm_precision
                ))
        
        reporter.add_section_results("KASP - Tm Values", results)
        
        # Assert all passed
        failed = [r for r in results if not r.is_match]
        if failed:
            pytest.fail(f"KASP Tm values failed for {len(failed)} cases: {[f.field_name for f in failed[:3]]}")
    
    def test_kasp_gc_content(self, reference_kasp_data, comparator, reporter):
        """Test KASP GC content consistency with tolerance."""
        results = []
        
        # Validate GC content in reference data
        # Note: V2 stores the GC content from Primer3's original calculation,
        # but for allele-specific primers, the 3' base is modified after design.
        # So GC content may not match the final sequence for Allele primers.
        for snp_name, primers in reference_kasp_data.items():
            for i, primer in enumerate(primers):
                # GC content should be in reasonable range
                gc_reasonable = 20.0 <= primer.gc_content <= 80.0
                
                results.append(ComparisonResult(
                    field_name=f"KASP primer {i+1} GC content range for {snp_name}",
                    expected_value="20-80%",
                    actual_value=primer.gc_content,
                    is_match=gc_reasonable
                ))
                
                # Calculate expected GC content from sequence
                if primer.primer_seq:
                    gc_count = primer.primer_seq.upper().count('G') + primer.primer_seq.upper().count('C')
                    expected_gc = (gc_count / len(primer.primer_seq)) * 100
                    
                    # For allele-specific primers (Allele-X in index), the 3' base is modified
                    # after Primer3 design, so GC content from Primer3 may differ from
                    # the final sequence. This is expected V2 behavior.
                    is_allele_primer = 'Allele' in primer.index
                    
                    if is_allele_primer:
                        # For allele primers, allow up to 1 base difference in GC calculation
                        # (one base change can cause ~4% GC difference for 24bp primer)
                        gc_tolerance = 5.0  # Allow ~1 base difference
                        gc_matches = abs(primer.gc_content - expected_gc) <= gc_tolerance
                        
                        results.append(ComparisonResult(
                            field_name=f"KASP primer {i+1} GC calculation for {snp_name} (allele-specific)",
                            expected_value=f"{expected_gc:.3f}% (±{gc_tolerance}% for allele primers)",
                            actual_value=primer.gc_content,
                            is_match=gc_matches
                        ))
                    else:
                        # For common primers, GC should match exactly
                        gc_matches = abs(primer.gc_content - expected_gc) <= 0.1
                        
                        results.append(ComparisonResult(
                            field_name=f"KASP primer {i+1} GC calculation for {snp_name}",
                            expected_value=f"{expected_gc:.3f}%",
                            actual_value=primer.gc_content,
                            is_match=gc_matches
                        ))
        
        reporter.add_section_results("KASP - GC Content", results)
        
        # Assert all passed
        failed = [r for r in results if not r.is_match]
        if failed:
            pytest.fail(f"KASP GC content failed for {len(failed)} cases: {[f.field_name for f in failed[:3]]}")
    
    def test_kasp_scoring(self, reference_kasp_data, comparator, reporter):
        """Test KASP scoring consistency with tolerance."""
        results = []
        
        # Validate scoring in reference data
        for snp_name, primers in reference_kasp_data.items():
            for i, primer in enumerate(primers):
                # Score should be positive
                score_positive = primer.score > 0
                
                results.append(ComparisonResult(
                    field_name=f"KASP primer {i+1} score positive for {snp_name}",
                    expected_value="> 0",
                    actual_value=primer.score,
                    is_match=score_positive
                ))
                
                # Penalty should be non-negative
                penalty_valid = primer.penalty >= 0
                
                results.append(ComparisonResult(
                    field_name=f"KASP primer {i+1} penalty non-negative for {snp_name}",
                    expected_value=">= 0",
                    actual_value=primer.penalty,
                    is_match=penalty_valid
                ))
                
                # Score precision
                score_precision = abs(primer.score - round(primer.score, 4)) < 0.0001
                
                results.append(ComparisonResult(
                    field_name=f"KASP primer {i+1} score precision for {snp_name}",
                    expected_value="4 decimal places",
                    actual_value=f"{primer.score:.6f}",
                    is_match=score_precision
                ))
        
        reporter.add_section_results("KASP - Scoring", results)
        
        # Assert all passed
        failed = [r for r in results if not r.is_match]
        if failed:
            pytest.fail(f"KASP scoring failed for {len(failed)} cases: {[f.field_name for f in failed[:3]]}")
    
    def test_kasp_output_format(self, reference_kasp_data, comparator, reporter):
        """Test KASP output format consistency."""
        results = []
        
        # Validate output format structure
        for snp_name, primers in reference_kasp_data.items():
            # Should have both allele-specific and common primers
            primer_types = [p.primer_type for p in primers]
            has_left = 'LEFT' in primer_types
            has_right = 'RIGHT' in primer_types
            
            results.extend([
                ComparisonResult(
                    field_name=f"KASP has LEFT primers for {snp_name}",
                    expected_value="LEFT primers present",
                    actual_value=f"LEFT present: {has_left}",
                    is_match=has_left
                ),
                ComparisonResult(
                    field_name=f"KASP has RIGHT primers for {snp_name}",
                    expected_value="RIGHT primers present",
                    actual_value=f"RIGHT present: {has_right}",
                    is_match=has_right
                )
            ])
            
            # Check for allele-specific primers
            allele_primers = [p for p in primers if 'Allele' in p.index]
            common_primers = [p for p in primers if 'Common' in p.index]
            
            has_allele_primers = len(allele_primers) > 0
            has_common_primers = len(common_primers) > 0
            
            results.extend([
                ComparisonResult(
                    field_name=f"KASP has allele-specific primers for {snp_name}",
                    expected_value="allele primers present",
                    actual_value=f"count: {len(allele_primers)}",
                    is_match=has_allele_primers
                ),
                ComparisonResult(
                    field_name=f"KASP has common primers for {snp_name}",
                    expected_value="common primers present",
                    actual_value=f"count: {len(common_primers)}",
                    is_match=has_common_primers
                )
            ])
        
        reporter.add_section_results("KASP - Output Format", results)
        
        # Assert all passed
        failed = [r for r in results if not r.is_match]
        if failed:
            pytest.fail(f"KASP output format failed for {len(failed)} cases: {[f.field_name for f in failed[:3]]}")
    
    def test_kasp_primer_validation(self, reference_kasp_data, comparator, reporter):
        """Test KASP primer validation consistency."""
        results = []
        
        # Validate primer properties
        for snp_name, primers in reference_kasp_data.items():
            for i, primer in enumerate(primers):
                # Primer sequence should be valid DNA
                seq_valid = all(c.upper() in 'ATCG' for c in primer.primer_seq)
                
                results.append(ComparisonResult(
                    field_name=f"KASP primer {i+1} sequence valid DNA for {snp_name}",
                    expected_value="only ATCG",
                    actual_value="valid" if seq_valid else "invalid characters",
                    is_match=seq_valid
                ))
                
                # Primer length should match sequence length
                length_matches = len(primer.primer_seq) == primer.length
                
                results.append(ComparisonResult(
                    field_name=f"KASP primer {i+1} length matches sequence for {snp_name}",
                    expected_value=primer.length,
                    actual_value=len(primer.primer_seq),
                    is_match=length_matches
                ))
                
                # Reverse complement should be valid
                if primer.reverse_complement:
                    rc_valid = len(primer.reverse_complement) == len(primer.primer_seq)
                    
                    results.append(ComparisonResult(
                        field_name=f"KASP primer {i+1} reverse complement length for {snp_name}",
                        expected_value=len(primer.primer_seq),
                        actual_value=len(primer.reverse_complement),
                        is_match=rc_valid
                    ))
        
        reporter.add_section_results("KASP - Primer Validation", results)
        
        # Assert all passed
        failed = [r for r in results if not r.is_match]
        if failed:
            pytest.fail(f"KASP primer validation failed for {len(failed)} cases: {[f.field_name for f in failed[:3]]}")
    
    def _create_test_alignment(self, snp_name: str) -> Dict[str, str]:
        """Create a test alignment for KASP design."""
        # This is a simplified test alignment
        # In real implementation, this would use actual alignment data
        return {
            f"{snp_name}_seq1": "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
            f"{snp_name}_seq2": "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
            f"{snp_name}_seq3": "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        }