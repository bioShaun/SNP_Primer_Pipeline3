"""Unit tests for MultipleSequenceAlignment variant detection (no MUSCLE/MAFFT)."""

from __future__ import annotations

import pytest

from snp_primer_pipeline.core.alignment import AlignedSequence, MultipleSequenceAlignment
from snp_primer_pipeline.exceptions import AlignmentError

pytestmark = pytest.mark.unit


def test_find_variant_sites_v2_no_homeologs_returns_empty() -> None:
    """Single target sequence: no comparison, empty diffarray."""
    msa = MultipleSequenceAlignment([AlignedSequence("t", "A" * 100)])
    msa.set_target_sequence("t")
    all_sites, any_sites, diff = msa.find_variant_sites_v2("t", 50)
    assert all_sites == []
    assert any_sites == []
    assert diff == {}


def test_find_variant_sites_v2_raises_if_target_missing() -> None:
    msa = MultipleSequenceAlignment(
        [
            AlignedSequence("a", "A" * 50),
            AlignedSequence("b", "A" * 50),
        ]
    )
    with pytest.raises(AlignmentError, match="not found"):
        msa.find_variant_sites_v2("missing", 10)


def test_find_variant_sites_v2_two_sequences_no_gaps() -> None:
    """Two aligned sequences with one difference at template position 50."""
    t = "A" * 50 + "C" + "A" * 49
    h = "A" * 50 + "G" + "A" * 49
    msa = MultipleSequenceAlignment(
        [
            AlignedSequence("target", t),
            AlignedSequence("homeo", h),
        ]
    )
    msa.set_target_sequence("target")
    _all_sites, any_sites, diff = msa.find_variant_sites_v2("target", 50)
    assert 50 in any_sites
    assert diff[50] == [1]


def test_find_variant_sites_gap_left_boundary() -> None:
    """Leading gaps in one sequence: gap_left > 0."""
    lead = "----" + "A" * 96
    msa = MultipleSequenceAlignment(
        [
            AlignedSequence("target", lead),
            AlignedSequence("homeo", "A" * 100),
        ]
    )
    msa.set_target_sequence("target")
    all_sites, any_sites, _ = msa.find_variant_sites_v2("target", 50)
    assert isinstance(all_sites, list)
    assert isinstance(any_sites, list)
