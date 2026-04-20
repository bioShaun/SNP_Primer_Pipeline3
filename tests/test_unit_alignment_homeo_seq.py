"""Unit tests for _get_homeo_seq and its helper functions."""

from __future__ import annotations

import pytest

from snp_primer_pipeline.core.alignment import (
    _find_longest_substring,
    _gap_diff,
    _get_homeo_seq,
    _score_pairwise,
)

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# _score_pairwise
# ---------------------------------------------------------------------------


class TestScorePairwise:
    def test_identical_sequences_positive_score(self) -> None:
        score = _score_pairwise("ATCGATCG", "ATCGATCG")
        assert score > 0
        assert score == 8.0  # 8 matches * 1.0

    def test_all_mismatches_negative_score(self) -> None:
        score = _score_pairwise("AAAA", "TTTT")
        assert score < 0
        assert score == -4.0  # 4 mismatches * -1.0

    def test_gap_penalties(self) -> None:
        # One gap opening + one gap extension
        score = _score_pairwise("AATT--CC", "AATTGGCC")
        # positions: 4 matches + gapopen + gapext + 2 matches
        # = 4*1.0 + (-4.0) + (-1.0) + 2*1.0 = 1.0
        assert score == 1.0

    def test_single_gap_opening(self) -> None:
        score = _score_pairwise("AAT-GG", "AATTGG")
        # 3 matches + gapopen + 2 matches = 3 + (-4) + 2 = 1.0
        assert score == 1.0


# ---------------------------------------------------------------------------
# _gap_diff
# ---------------------------------------------------------------------------


class TestGapDiff:
    def test_no_gaps_returns_zero(self) -> None:
        assert _gap_diff("ATCG", "ATCG") == 0

    def test_one_sided_gap_counted(self) -> None:
        assert _gap_diff("AT-G", "ATCG") == 1

    def test_multiple_one_sided_gaps(self) -> None:
        assert _gap_diff("A--G", "ATCG") == 2

    def test_both_sides_gap_not_counted(self) -> None:
        # pair "--" has count("-") == 2, not 1, so not counted
        assert _gap_diff("A-CG", "A-CG") == 0

    def test_mixed_gaps(self) -> None:
        # pos 1: both gap (0), pos 2: one-sided (1), pos 3: one-sided (1)
        assert _gap_diff("A--G", "A-CG") == 1


# ---------------------------------------------------------------------------
# _find_longest_substring
# ---------------------------------------------------------------------------


class TestFindLongestSubstring:
    def test_no_gaps_returns_whole_sequence(self) -> None:
        start, end, n_l, n_r = _find_longest_substring("ATCGATCG", "ATCGATCG")
        assert start == 0
        assert end == 8
        assert n_l == 0
        assert n_r == 0

    def test_one_gap_returns_longest_segment(self) -> None:
        # s1 = "AT-GATCG", gap at position 2
        # segments: [0:2] = "AT" Tm=4, [3:8] = "GATCG" Tm=14
        start, end, n_l, n_r = _find_longest_substring("AT-GATCG", "ATCGATCG")
        assert start == 3
        assert end == 8
        # n_left = chars before position 3 in s1 without gaps = "AT" → 2
        assert n_l == 2
        # n_right = chars from position 8 onward = "" → 0
        assert n_r == 0

    def test_gap_in_s2_also_detected(self) -> None:
        # gap at position 2 (s2 has "-")
        start, end, n_l, n_r = _find_longest_substring("ATCGATCG", "AT-GATCG")
        assert start == 3
        assert end == 8
        # n_left = len("ATC".replace("-", "")) = 3 (s1 chars before longest_start)
        assert n_l == 3
        assert n_r == 0


# ---------------------------------------------------------------------------
# _get_homeo_seq
# ---------------------------------------------------------------------------


class TestGetHomeoSeq:
    def test_single_homeolog_no_gaps(self) -> None:
        """No gaps → gap-free comparison returns the aligned region directly."""
        fasta = {
            "target": "AAAAATTTTTCCCCC",
            "homeo1": "AAAAAGGGGGCCCCC",
        }
        result = _get_homeo_seq(fasta, "target", ["homeo1"], 0, 14)
        assert len(result) == 1
        assert result[0] == "AAAAAGGGGGCCCCC"

    def test_with_gaps_and_good_alignment_score(self) -> None:
        """score1 > score2 and gap_diff < 4 → use gap-removed version."""
        # Build sequences where aligned score beats reconstructed score.
        # target and homeo are identical except for one gap in homeo.
        # This forces score1 (with gap penalty) to still beat score2.
        target = "ATCGATCGATCG"
        homeo = "ATC-ATCGATCG"
        fasta = {"target": target, "homeo1": homeo}
        result = _get_homeo_seq(fasta, "target", ["homeo1"], 0, 11)
        assert len(result) == 1
        # Verify the branch condition holds
        target_seq = target
        homeo_seq = homeo
        from snp_primer_pipeline.core.alignment import _find_longest_substring

        score1 = _score_pairwise(target_seq, homeo_seq)
        idx_l, idx_r, n_l, n_r = _find_longest_substring(target_seq, homeo_seq)
        s2 = homeo
        seq_l = s2[:idx_l].replace("-", "")
        seq_r = s2[idx_r:].replace("-", "")
        if len(seq_l) < n_l:
            seq_l = "-" * (n_l - len(seq_l)) + seq_l
        if len(seq_r) < n_r:
            seq_r = seq_r + "-" * (n_r - len(seq_r))
        seqk = seq_l[::-1][:n_l][::-1] + s2[idx_l:idx_r] + seq_r[:n_r]
        score2 = _score_pairwise(target_seq.replace("-", ""), seqk)
        assert score1 > score2, "test setup: score1 must beat score2"
        assert _gap_diff(target_seq, homeo_seq) < 4, "test setup: gap_diff must be < 4"
        # gap-removed: take homeo[i] where target[i] != "-"
        expected = "".join(c for i, c in enumerate(homeo_seq) if target_seq[i] != "-")
        assert result[0] == expected

    def test_with_gaps_and_poor_score(self) -> None:
        """score1 <= score2 → uses the reconstructed (seqk) version."""
        # Construct a scenario where aligned score is worse than reconstructed
        # Target with many gaps → poor alignment score, gap_diff >= 4
        target = "A---A---A---TTTT"
        homeo = "ACCCAGGGACCCTTTT"
        fasta = {"target": target, "homeo1": homeo}
        result = _get_homeo_seq(fasta, "target", ["homeo1"], 0, 15)
        assert len(result) == 1
        # gap_diff >= 4, so the gap-removed branch is NOT taken
        # Result is the reconstructed seqk (length matches target without gaps)
        assert len(result[0]) == len(target.replace("-", ""))

    def test_multiple_homeologs(self) -> None:
        """Returns one entry per homeolog."""
        fasta = {
            "target": "AAAATTTTCCCC",
            "h1": "AAAATTTTCCCC",
            "h2": "AAAATTTTGGGG",
        }
        result = _get_homeo_seq(fasta, "target", ["h1", "h2"], 0, 11)
        assert len(result) == 2
