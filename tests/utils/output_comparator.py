#!/usr/bin/env python3
"""
Output comparator for consistency tests.

Compares V3 pipeline outputs with V2 reference data.
"""

from dataclasses import dataclass
from typing import List, Dict, Any, Optional, Union
import math
from .reference_loader import KASPPrimerRecord, CAPSPrimerRecord, BlastHitRecord


@dataclass
class ComparisonResult:
    """Result of comparing two values."""
    field_name: str
    expected_value: Any
    actual_value: Any
    is_match: bool
    tolerance: Optional[float] = None
    message: str = ""
    
    def __post_init__(self):
        if not self.message:
            if self.is_match:
                self.message = f"✓ {self.field_name} matches"
            else:
                if self.tolerance is not None:
                    self.message = f"✗ {self.field_name}: expected {self.expected_value}, got {self.actual_value} (tolerance: {self.tolerance})"
                else:
                    self.message = f"✗ {self.field_name}: expected {self.expected_value}, got {self.actual_value}"


class OutputComparator:
    """Compares V3 outputs with V2 reference data."""
    
    # Numerical comparison tolerances
    TM_TOLERANCE = 0.001
    GC_TOLERANCE = 0.001
    PENALTY_TOLERANCE = 0.0001
    SCORE_TOLERANCE = 0.0001
    EVALUE_TOLERANCE = 1e-10
    BIT_SCORE_TOLERANCE = 0.001
    
    def compare_kasp_primers(
        self,
        expected: List[KASPPrimerRecord],
        actual: List[KASPPrimerRecord],
        snp_name: str = ""
    ) -> List[ComparisonResult]:
        """Compare KASP primer lists."""
        results = []
        
        # Check count
        results.append(ComparisonResult(
            field_name=f"KASP primer count for {snp_name}",
            expected_value=len(expected),
            actual_value=len(actual),
            is_match=len(expected) == len(actual)
        ))
        
        if len(expected) != len(actual):
            return results
        
        # Sort both lists by index for consistent comparison
        expected_sorted = sorted(expected, key=lambda x: x.index)
        actual_sorted = sorted(actual, key=lambda x: x.index)
        
        for i, (exp, act) in enumerate(zip(expected_sorted, actual_sorted)):
            prefix = f"KASP primer {i+1}"
            
            # Compare all fields
            results.extend([
                ComparisonResult(
                    field_name=f"{prefix} index",
                    expected_value=exp.index,
                    actual_value=act.index,
                    is_match=exp.index == act.index
                ),
                ComparisonResult(
                    field_name=f"{prefix} product_size",
                    expected_value=exp.product_size,
                    actual_value=act.product_size,
                    is_match=exp.product_size == act.product_size
                ),
                ComparisonResult(
                    field_name=f"{prefix} primer_type",
                    expected_value=exp.primer_type,
                    actual_value=act.primer_type,
                    is_match=exp.primer_type == act.primer_type
                ),
                ComparisonResult(
                    field_name=f"{prefix} start",
                    expected_value=exp.start,
                    actual_value=act.start,
                    is_match=exp.start == act.start
                ),
                ComparisonResult(
                    field_name=f"{prefix} end",
                    expected_value=exp.end,
                    actual_value=act.end,
                    is_match=exp.end == act.end
                ),
                ComparisonResult(
                    field_name=f"{prefix} variation_number",
                    expected_value=exp.variation_number,
                    actual_value=act.variation_number,
                    is_match=exp.variation_number == act.variation_number
                ),
                ComparisonResult(
                    field_name=f"{prefix} diff_three_all",
                    expected_value=exp.diff_three_all,
                    actual_value=act.diff_three_all,
                    is_match=exp.diff_three_all == act.diff_three_all
                ),
                ComparisonResult(
                    field_name=f"{prefix} length",
                    expected_value=exp.length,
                    actual_value=act.length,
                    is_match=exp.length == act.length
                ),
                self.compare_numeric(
                    expected=exp.tm,
                    actual=act.tm,
                    field_name=f"{prefix} tm",
                    tolerance=self.TM_TOLERANCE
                ),
                self.compare_numeric(
                    expected=exp.gc_content,
                    actual=act.gc_content,
                    field_name=f"{prefix} gc_content",
                    tolerance=self.GC_TOLERANCE
                ),
                self.compare_numeric(
                    expected=exp.self_any,
                    actual=act.self_any,
                    field_name=f"{prefix} self_any",
                    tolerance=self.SCORE_TOLERANCE
                ),
                self.compare_numeric(
                    expected=exp.self_three,
                    actual=act.self_three,
                    field_name=f"{prefix} self_three",
                    tolerance=self.SCORE_TOLERANCE
                ),
                self.compare_numeric(
                    expected=exp.end_stability,
                    actual=act.end_stability,
                    field_name=f"{prefix} end_stability",
                    tolerance=self.SCORE_TOLERANCE
                ),
                self.compare_numeric(
                    expected=exp.hairpin,
                    actual=act.hairpin,
                    field_name=f"{prefix} hairpin",
                    tolerance=self.SCORE_TOLERANCE
                ),
                ComparisonResult(
                    field_name=f"{prefix} primer_seq",
                    expected_value=exp.primer_seq,
                    actual_value=act.primer_seq,
                    is_match=exp.primer_seq == act.primer_seq
                ),
                ComparisonResult(
                    field_name=f"{prefix} reverse_complement",
                    expected_value=exp.reverse_complement,
                    actual_value=act.reverse_complement,
                    is_match=exp.reverse_complement == act.reverse_complement
                ),
                self.compare_numeric(
                    expected=exp.penalty,
                    actual=act.penalty,
                    field_name=f"{prefix} penalty",
                    tolerance=self.PENALTY_TOLERANCE
                ),
                self.compare_numeric(
                    expected=exp.compl_any,
                    actual=act.compl_any,
                    field_name=f"{prefix} compl_any",
                    tolerance=self.SCORE_TOLERANCE
                ),
                self.compare_numeric(
                    expected=exp.compl_end,
                    actual=act.compl_end,
                    field_name=f"{prefix} compl_end",
                    tolerance=self.SCORE_TOLERANCE
                ),
                self.compare_numeric(
                    expected=exp.score,
                    actual=act.score,
                    field_name=f"{prefix} score",
                    tolerance=self.SCORE_TOLERANCE
                )
            ])
        
        return results
    
    def compare_caps_primers(
        self,
        expected: List[CAPSPrimerRecord],
        actual: List[CAPSPrimerRecord],
        snp_name: str = ""
    ) -> List[ComparisonResult]:
        """Compare CAPS primer lists."""
        results = []
        
        # Check count
        results.append(ComparisonResult(
            field_name=f"CAPS primer count for {snp_name}",
            expected_value=len(expected),
            actual_value=len(actual),
            is_match=len(expected) == len(actual)
        ))
        
        if len(expected) != len(actual):
            return results
        
        # Sort both lists by index for consistent comparison
        expected_sorted = sorted(expected, key=lambda x: x.index)
        actual_sorted = sorted(actual, key=lambda x: x.index)
        
        for i, (exp, act) in enumerate(zip(expected_sorted, actual_sorted)):
            prefix = f"CAPS primer {i+1}"
            
            # Compare all fields
            results.extend([
                ComparisonResult(
                    field_name=f"{prefix} index",
                    expected_value=exp.index,
                    actual_value=act.index,
                    is_match=exp.index == act.index
                ),
                ComparisonResult(
                    field_name=f"{prefix} product_size",
                    expected_value=exp.product_size,
                    actual_value=act.product_size,
                    is_match=exp.product_size == act.product_size
                ),
                ComparisonResult(
                    field_name=f"{prefix} primer_type",
                    expected_value=exp.primer_type,
                    actual_value=act.primer_type,
                    is_match=exp.primer_type == act.primer_type
                ),
                ComparisonResult(
                    field_name=f"{prefix} start",
                    expected_value=exp.start,
                    actual_value=act.start,
                    is_match=exp.start == act.start
                ),
                ComparisonResult(
                    field_name=f"{prefix} end",
                    expected_value=exp.end,
                    actual_value=act.end,
                    is_match=exp.end == act.end
                ),
                ComparisonResult(
                    field_name=f"{prefix} diff_number",
                    expected_value=exp.diff_number,
                    actual_value=act.diff_number,
                    is_match=exp.diff_number == act.diff_number
                ),
                ComparisonResult(
                    field_name=f"{prefix} diff_three_all",
                    expected_value=exp.diff_three_all,
                    actual_value=act.diff_three_all,
                    is_match=exp.diff_three_all == act.diff_three_all
                ),
                ComparisonResult(
                    field_name=f"{prefix} length",
                    expected_value=exp.length,
                    actual_value=act.length,
                    is_match=exp.length == act.length
                ),
                self.compare_numeric(
                    expected=exp.tm,
                    actual=act.tm,
                    field_name=f"{prefix} tm",
                    tolerance=self.TM_TOLERANCE
                ),
                self.compare_numeric(
                    expected=exp.gc_content,
                    actual=act.gc_content,
                    field_name=f"{prefix} gc_content",
                    tolerance=self.GC_TOLERANCE
                ),
                self.compare_numeric(
                    expected=exp.self_any,
                    actual=act.self_any,
                    field_name=f"{prefix} self_any",
                    tolerance=self.SCORE_TOLERANCE
                ),
                self.compare_numeric(
                    expected=exp.self_three,
                    actual=act.self_three,
                    field_name=f"{prefix} self_three",
                    tolerance=self.SCORE_TOLERANCE
                ),
                self.compare_numeric(
                    expected=exp.end_stability,
                    actual=act.end_stability,
                    field_name=f"{prefix} end_stability",
                    tolerance=self.SCORE_TOLERANCE
                ),
                self.compare_numeric(
                    expected=exp.hairpin,
                    actual=act.hairpin,
                    field_name=f"{prefix} hairpin",
                    tolerance=self.SCORE_TOLERANCE
                ),
                ComparisonResult(
                    field_name=f"{prefix} primer_seq",
                    expected_value=exp.primer_seq,
                    actual_value=act.primer_seq,
                    is_match=exp.primer_seq == act.primer_seq
                ),
                ComparisonResult(
                    field_name=f"{prefix} reverse_complement",
                    expected_value=exp.reverse_complement,
                    actual_value=act.reverse_complement,
                    is_match=exp.reverse_complement == act.reverse_complement
                ),
                self.compare_numeric(
                    expected=exp.penalty,
                    actual=act.penalty,
                    field_name=f"{prefix} penalty",
                    tolerance=self.PENALTY_TOLERANCE
                ),
                self.compare_numeric(
                    expected=exp.compl_any,
                    actual=act.compl_any,
                    field_name=f"{prefix} compl_any",
                    tolerance=self.SCORE_TOLERANCE
                ),
                self.compare_numeric(
                    expected=exp.compl_end,
                    actual=act.compl_end,
                    field_name=f"{prefix} compl_end",
                    tolerance=self.SCORE_TOLERANCE
                )
            ])
        
        return results
    
    def compare_blast_results(
        self,
        expected: List[BlastHitRecord],
        actual: List[BlastHitRecord],
        query_id: str = ""
    ) -> List[ComparisonResult]:
        """Compare BLAST result lists."""
        results = []
        
        # Check count
        results.append(ComparisonResult(
            field_name=f"BLAST hit count for {query_id}",
            expected_value=len(expected),
            actual_value=len(actual),
            is_match=len(expected) == len(actual)
        ))
        
        if len(expected) != len(actual):
            return results
        
        # Sort both lists by subject_id and coordinates for consistent comparison
        expected_sorted = sorted(expected, key=lambda x: (x.subject_id, x.subject_start, x.subject_end))
        actual_sorted = sorted(actual, key=lambda x: (x.subject_id, x.subject_start, x.subject_end))
        
        for i, (exp, act) in enumerate(zip(expected_sorted, actual_sorted)):
            prefix = f"BLAST hit {i+1}"
            
            # Compare all fields
            results.extend([
                ComparisonResult(
                    field_name=f"{prefix} query_id",
                    expected_value=exp.query_id,
                    actual_value=act.query_id,
                    is_match=exp.query_id == act.query_id
                ),
                ComparisonResult(
                    field_name=f"{prefix} subject_id",
                    expected_value=exp.subject_id,
                    actual_value=act.subject_id,
                    is_match=exp.subject_id == act.subject_id
                ),
                self.compare_numeric(
                    expected=exp.identity,
                    actual=act.identity,
                    field_name=f"{prefix} identity",
                    tolerance=0.01  # 1% tolerance for identity
                ),
                ComparisonResult(
                    field_name=f"{prefix} alignment_length",
                    expected_value=exp.alignment_length,
                    actual_value=act.alignment_length,
                    is_match=exp.alignment_length == act.alignment_length
                ),
                ComparisonResult(
                    field_name=f"{prefix} mismatches",
                    expected_value=exp.mismatches,
                    actual_value=act.mismatches,
                    is_match=exp.mismatches == act.mismatches
                ),
                ComparisonResult(
                    field_name=f"{prefix} gap_opens",
                    expected_value=exp.gap_opens,
                    actual_value=act.gap_opens,
                    is_match=exp.gap_opens == act.gap_opens
                ),
                ComparisonResult(
                    field_name=f"{prefix} query_start",
                    expected_value=exp.query_start,
                    actual_value=act.query_start,
                    is_match=exp.query_start == act.query_start
                ),
                ComparisonResult(
                    field_name=f"{prefix} query_end",
                    expected_value=exp.query_end,
                    actual_value=act.query_end,
                    is_match=exp.query_end == act.query_end
                ),
                ComparisonResult(
                    field_name=f"{prefix} subject_start",
                    expected_value=exp.subject_start,
                    actual_value=act.subject_start,
                    is_match=exp.subject_start == act.subject_start
                ),
                ComparisonResult(
                    field_name=f"{prefix} subject_end",
                    expected_value=exp.subject_end,
                    actual_value=act.subject_end,
                    is_match=exp.subject_end == act.subject_end
                ),
                self.compare_numeric(
                    expected=exp.evalue,
                    actual=act.evalue,
                    field_name=f"{prefix} evalue",
                    tolerance=self.EVALUE_TOLERANCE
                ),
                self.compare_numeric(
                    expected=exp.bit_score,
                    actual=act.bit_score,
                    field_name=f"{prefix} bit_score",
                    tolerance=self.BIT_SCORE_TOLERANCE
                )
            ])
        
        return results
    
    def compare_alignments(
        self,
        expected: Dict[str, str],
        actual: Dict[str, str],
        snp_name: str = ""
    ) -> List[ComparisonResult]:
        """Compare multiple sequence alignments."""
        results = []
        
        # Check sequence count
        results.append(ComparisonResult(
            field_name=f"Alignment sequence count for {snp_name}",
            expected_value=len(expected),
            actual_value=len(actual),
            is_match=len(expected) == len(actual)
        ))
        
        # Check sequence names
        expected_names = set(expected.keys())
        actual_names = set(actual.keys())
        
        results.append(ComparisonResult(
            field_name=f"Alignment sequence names for {snp_name}",
            expected_value=sorted(expected_names),
            actual_value=sorted(actual_names),
            is_match=expected_names == actual_names
        ))
        
        # Compare sequences for common names
        common_names = expected_names & actual_names
        for name in sorted(common_names):
            results.append(self.compare_sequences(
                expected=expected[name],
                actual=actual[name],
                name=f"Alignment sequence {name}"
            ))
        
        return results
    
    def compare_variation_sites(
        self,
        expected: List[int],
        actual: List[int],
        snp_name: str = ""
    ) -> List[ComparisonResult]:
        """Compare variation site lists."""
        results = []
        
        # Sort both lists for comparison
        expected_sorted = sorted(expected)
        actual_sorted = sorted(actual)
        
        results.append(ComparisonResult(
            field_name=f"Variation sites for {snp_name}",
            expected_value=expected_sorted,
            actual_value=actual_sorted,
            is_match=expected_sorted == actual_sorted
        ))
        
        return results
    
    def compare_sequences(
        self,
        expected: str,
        actual: str,
        name: str
    ) -> ComparisonResult:
        """Compare two sequences."""
        return ComparisonResult(
            field_name=name,
            expected_value=expected,
            actual_value=actual,
            is_match=expected == actual
        )
    
    def compare_numeric(
        self,
        expected: Union[int, float],
        actual: Union[int, float],
        field_name: str,
        tolerance: float
    ) -> ComparisonResult:
        """Compare numeric values with tolerance."""
        if isinstance(expected, (int, float)) and isinstance(actual, (int, float)):
            diff = abs(expected - actual)
            is_match = diff <= tolerance
        else:
            is_match = expected == actual
        
        return ComparisonResult(
            field_name=field_name,
            expected_value=expected,
            actual_value=actual,
            is_match=is_match,
            tolerance=tolerance
        )
    
    def get_failed_results(self, results: List[ComparisonResult]) -> List[ComparisonResult]:
        """Get only the failed comparison results."""
        return [r for r in results if not r.is_match]
    
    def get_success_rate(self, results: List[ComparisonResult]) -> float:
        """Calculate success rate as percentage."""
        if not results:
            return 100.0
        
        passed = sum(1 for r in results if r.is_match)
        return (passed / len(results)) * 100.0