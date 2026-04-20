"""Unit tests for MultipleSequenceAlignment.find_variant_sites (V1 logic)."""

from __future__ import annotations

import pytest

from snp_primer_pipeline.core.alignment import AlignedSequence, MultipleSequenceAlignment
from snp_primer_pipeline.exceptions import AlignmentError

pytestmark = pytest.mark.unit


def _pad(core: str, total: int = 100) -> str:
    """Pad a core sequence with 'A' to reach total length (variant at center)."""
    left = (total - len(core)) // 2
    right = total - len(core) - left
    return "A" * left + core + "A" * right


# ── basic: no other sequences ────────────────────────────────────────


def test_single_sequence_returns_empty() -> None:
    msa = MultipleSequenceAlignment([AlignedSequence("t", "A" * 100)])
    msa.set_target_sequence("t")
    diff_all, diff_any = msa.find_variant_sites("t")
    assert diff_all == []
    assert diff_any == []


# ── target not found ─────────────────────────────────────────────────


def test_target_not_found_raises() -> None:
    msa = MultipleSequenceAlignment([AlignedSequence("a", "A" * 60)])
    with pytest.raises(AlignmentError, match="not found"):
        msa.find_variant_sites("missing")


# ── two sequences, one difference ────────────────────────────────────


def test_two_sequences_one_diff() -> None:
    """Difference at center (pos 50) should appear in both diff_all and diff_any."""
    t_seq = "A" * 50 + "C" + "A" * 49  # pos 50 = C
    h_seq = "A" * 50 + "G" + "A" * 49  # pos 50 = G
    msa = MultipleSequenceAlignment(
        [AlignedSequence("target", t_seq), AlignedSequence("homeo", h_seq)]
    )
    msa.set_target_sequence("target")
    diff_all, diff_any = msa.find_variant_sites("target")
    assert 50 in diff_all
    assert 50 in diff_any


# ── boundary exclusion (pos < 20 or pos > len-20) ───────────────────


def test_boundary_positions_excluded() -> None:
    """Variants within 20bp of either end are excluded by V2-style boundary check."""
    # Diff at pos 5 (< 20) and pos 95 (> 100-20=80): both excluded
    t_seq = "A" * 5 + "C" + "A" * 89 + "C" + "A" * 4  # diffs at 5 and 95
    h_seq = "A" * 5 + "G" + "A" * 89 + "G" + "A" * 4
    msa = MultipleSequenceAlignment(
        [AlignedSequence("target", t_seq), AlignedSequence("homeo", h_seq)]
    )
    msa.set_target_sequence("target")
    _diff_all, diff_any = msa.find_variant_sites("target")
    assert 5 not in diff_any
    assert 95 not in diff_any


# ── three sequences: diff_all vs diff_any ────────────────────────────


def test_three_sequences_diff_all_vs_diff_any() -> None:
    """Position where target differs from one but matches another:
    in diff_any but not diff_all."""
    t_seq = "A" * 50 + "C" + "A" * 49
    h1 = "A" * 50 + "G" + "A" * 49  # differs from target
    h2 = "A" * 50 + "C" + "A" * 49  # same as target
    msa = MultipleSequenceAlignment(
        [
            AlignedSequence("target", t_seq),
            AlignedSequence("h1", h1),
            AlignedSequence("h2", h2),
        ]
    )
    msa.set_target_sequence("target")
    diff_all, diff_any = msa.find_variant_sites("target")
    assert 50 in diff_any
    assert 50 not in diff_all


def test_three_sequences_diff_from_all() -> None:
    """Position where target differs from ALL others: in both lists."""
    t_seq = "A" * 50 + "C" + "A" * 49
    h1 = "A" * 50 + "G" + "A" * 49
    h2 = "A" * 50 + "T" + "A" * 49
    msa = MultipleSequenceAlignment(
        [
            AlignedSequence("target", t_seq),
            AlignedSequence("h1", h1),
            AlignedSequence("h2", h2),
        ]
    )
    msa.set_target_sequence("target")
    diff_all, diff_any = msa.find_variant_sites("target")
    assert 50 in diff_all
    assert 50 in diff_any


# ── gaps in target are skipped ───────────────────────────────────────


def test_gaps_in_target_skipped() -> None:
    """Gap in target at alignment position → that position is skipped."""
    # Insert a gap in target at alignment pos 50; homeolog has 'G' there.
    t_seq = "A" * 50 + "-" + "A" * 49  # 101 chars aligned, 100 template
    h_seq = "A" * 50 + "G" + "A" * 49
    msa = MultipleSequenceAlignment(
        [AlignedSequence("target", t_seq), AlignedSequence("homeo", h_seq)]
    )
    msa.set_target_sequence("target")
    _diff_all, diff_any = msa.find_variant_sites("target")
    # The gap position itself produces no template coordinate → no variant
    # Template pos 50 maps to alignment pos 51 (shifted by gap), both are 'A' → no diff
    assert len(diff_any) == 0


# ── gaps in other sequences are skipped ──────────────────────────────


def test_gaps_in_other_sequence_skipped() -> None:
    """Gap in homeolog at comparison position → that comparison is skipped."""
    t_seq = "A" * 50 + "C" + "A" * 49
    h_seq = "A" * 50 + "-" + "A" * 49  # gap at alignment pos 50
    msa = MultipleSequenceAlignment(
        [AlignedSequence("target", t_seq), AlignedSequence("homeo", h_seq)]
    )
    msa.set_target_sequence("target")
    diff_all, diff_any = msa.find_variant_sites("target")
    # Gap in other seq → comparison skipped → no diff recorded
    # diff_from_all stays True (no non-gap other matched), diff_from_any stays False
    # Actually: the loop skips gaps with `continue`, so b1!=b2 never fires,
    # but diff_from_all remains True (default) → pos in diff_all but not diff_any.
    # This is debatable V2 behavior; let's verify what actually happens.
    # If only one other seq and it's a gap, diff_from_any=False, diff_from_all=True
    assert 50 in diff_all
    assert 50 not in diff_any


# ── leading gaps affect gap_left boundary ────────────────────────────


def test_leading_gaps_shift_gap_left() -> None:
    """Leading gaps in sequences shift gap_left, excluding early positions."""
    # homeolog has 30 leading gaps → gap_left = 30
    t_seq = "A" * 100
    h_seq = "-" * 30 + "A" * 20 + "C" + "A" * 49  # diff at aligned pos 50
    msa = MultipleSequenceAlignment(
        [AlignedSequence("target", t_seq), AlignedSequence("homeo", h_seq)]
    )
    msa.set_target_sequence("target")
    _diff_all, diff_any = msa.find_variant_sites("target")
    # The diff at alignment pos 50 → template pos 50 → should be in range
    assert 50 in diff_any


# ── results are sorted ───────────────────────────────────────────────


def test_results_are_sorted() -> None:
    """Output lists are sorted by template position."""
    t_seq = "A" * 40 + "C" + "A" * 9 + "C" + "A" * 49  # diffs at 40, 50
    h_seq = "A" * 40 + "G" + "A" * 9 + "G" + "A" * 49
    msa = MultipleSequenceAlignment(
        [AlignedSequence("target", t_seq), AlignedSequence("homeo", h_seq)]
    )
    msa.set_target_sequence("target")
    diff_all, diff_any = msa.find_variant_sites("target")
    assert diff_all == sorted(diff_all)
    assert diff_any == sorted(diff_any)
    assert 40 in diff_any
    assert 50 in diff_any
