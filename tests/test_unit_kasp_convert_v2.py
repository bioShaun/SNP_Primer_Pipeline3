"""Unit tests for KASPDesigner._convert_to_v2_format / _format_primer_pair_v2."""

from __future__ import annotations

from pathlib import Path

import pytest

from dataclasses import asdict

from snp_primer_pipeline.models import Primer, PrimerPair, Strand
from snp_primer_pipeline.primers.kasp import KASPDesigner

pytestmark = pytest.mark.unit


@pytest.fixture
def designer(tmp_path: Path) -> KASPDesigner:
    fake = tmp_path / "primer3_core"
    fake.write_text("#!/bin/sh\necho\n")
    fake.chmod(0o755)
    d = KASPDesigner(primer3_path=fake, config_path=None)
    # Defaults: no genomic info
    d._genomic_start = None
    d._genomic_strand = None
    d._template_sequence = "a" * 200
    return d


def _make_pair(
    left_start: int = 10,
    left_end: int = 30,
    right_start: int = 150,
    right_end: int = 170,
) -> PrimerPair:
    left = Primer(
        name="L",
        sequence="a" * (left_end - left_start + 1),
        start=left_start,
        end=left_end,
        length=left_end - left_start + 1,
        direction="LEFT",
    )
    right = Primer(
        name="R",
        sequence="t" * (right_end - right_start + 1),
        start=right_start,
        end=right_end,
        length=right_end - right_start + 1,
        direction="RIGHT",
    )
    return PrimerPair(
        left=left,
        right=right,
        product_size=right_end - left_start,
        penalty=1.5,
        compl_any=0.1,
        compl_end=0.2,
        score=2.0,
    )


# ── _convert_to_v2_format ────────────────────────────────────────────


def test_convert_to_v2_format_returns_three_rows_per_pair(designer: KASPDesigner) -> None:
    """Each primer pair produces 3 output rows: Allele-A, Allele-B, Common."""
    pp = _make_pair()
    final = {"31-0": (pp, 30)}
    results = designer._convert_to_v2_format(
        final, snp_position=50, snp_alleles=("A", "G"), variant_sites=[30]
    )
    assert len(results) == 3


def test_convert_to_v2_format_multiple_pairs(designer: KASPDesigner) -> None:
    """Two pairs produce 6 rows."""
    pp1 = _make_pair()
    pp2 = _make_pair(left_start=20, left_end=40, right_start=160, right_end=180)
    final = {"31-0": (pp1, 30), "41-1": (pp2, 40)}
    results = designer._convert_to_v2_format(
        final, snp_position=50, snp_alleles=("A", "G"), variant_sites=[30, 40]
    )
    assert len(results) == 6


# ── _format_primer_pair_v2: left allele-specific (left.end == snp) ───


def test_format_pair_left_allele_specific(designer: KASPDesigner) -> None:
    """When left primer 3' end == snp_position, left is allele-specific."""
    snp_pos = 30
    pp = _make_pair(left_start=10, left_end=snp_pos)
    rows = designer._format_primer_pair_v2(
        pp, snp_position=snp_pos, snp_alleles=("A", "G"), variant_sites=[], key="31-0"
    )
    assert len(rows) == 3
    assert "Allele-A" in rows[0].index
    assert "Allele-G" in rows[1].index
    assert "Common" in rows[2].index

    # Allele-specific rows should be LEFT direction
    assert rows[0].direction == "LEFT"
    assert rows[1].direction == "LEFT"
    # Common should be RIGHT
    assert rows[2].direction == "RIGHT"

    # Allele-specific sequences should end with the respective allele
    assert rows[0].primer_seq.endswith("A")
    assert rows[1].primer_seq.endswith("G")


def test_format_pair_right_allele_specific(designer: KASPDesigner) -> None:
    """When left primer 3' end != snp_position, right is allele-specific."""
    snp_pos = 50
    pp = _make_pair(left_start=10, left_end=30, right_start=150, right_end=170)
    rows = designer._format_primer_pair_v2(
        pp, snp_position=snp_pos, snp_alleles=("A", "G"), variant_sites=[], key="31-0"
    )
    # Allele-specific should be RIGHT direction
    assert rows[0].direction == "RIGHT"
    assert rows[1].direction == "RIGHT"
    # Common is LEFT
    assert rows[2].direction == "LEFT"


# ── coordinate format ────────────────────────────────────────────────


def test_format_pair_left_coords_1based(designer: KASPDesigner) -> None:
    """LEFT primer: start_out = start+1, end_out = end+1 (1-based)."""
    pp = _make_pair(left_start=10, left_end=30, right_start=150, right_end=170)
    rows = designer._format_primer_pair_v2(
        pp, snp_position=30, snp_alleles=("A", "G"), variant_sites=[], key="31-0"
    )
    # The Common row is RIGHT (since left.end == snp)
    # Allele rows are LEFT
    allele_row = rows[0]
    assert allele_row.direction == "LEFT"
    assert allele_row.start == 11  # 10 + 1
    assert allele_row.end == 31  # 30 + 1


def test_format_pair_right_coords_swapped(designer: KASPDesigner) -> None:
    """RIGHT primer: start_out = end+1, end_out = start+1 (V2 swap)."""
    pp = _make_pair(left_start=10, left_end=30, right_start=150, right_end=170)
    rows = designer._format_primer_pair_v2(
        pp, snp_position=30, snp_alleles=("A", "G"), variant_sites=[], key="31-0"
    )
    common_row = rows[2]  # Common is RIGHT
    assert common_row.direction == "RIGHT"
    assert common_row.start == 171  # end+1 (5' end, higher)
    assert common_row.end == 151  # start+1 (3' end, lower)


# ── diff_three_all / diff_num ────────────────────────────────────────


def test_format_pair_diff_three_all_yes(designer: KASPDesigner) -> None:
    """var_site parsed from key is in variant_sites -> diff_three_all = YES."""
    pp = _make_pair()
    # key="31-0" → var_site_from_key = 30 (0-based)
    rows = designer._format_primer_pair_v2(
        pp, snp_position=50, snp_alleles=("A", "G"), variant_sites=[30], key="31-0"
    )
    assert rows[0].diff_three_all == "YES"


def test_format_pair_diff_three_all_no(designer: KASPDesigner) -> None:
    """var_site parsed from key is NOT in variant_sites -> diff_three_all = NO."""
    pp = _make_pair()
    rows = designer._format_primer_pair_v2(
        pp, snp_position=50, snp_alleles=("A", "G"), variant_sites=[99], key="31-0"
    )
    assert rows[0].diff_three_all == "NO"


def test_format_pair_invalid_key_diff_three_all_no(designer: KASPDesigner) -> None:
    """Malformed key → var_site_from_key = -1 → diff_three_all = NO."""
    pp = _make_pair()
    rows = designer._format_primer_pair_v2(
        pp, snp_position=50, snp_alleles=("A", "G"), variant_sites=[30], key="bad-key"
    )
    assert rows[0].diff_three_all == "NO"


# ── genomic coordinates ──────────────────────────────────────────────


def test_format_pair_no_genomic_info(designer: KASPDesigner) -> None:
    """When _genomic_start is None, genomic coords are None."""
    pp = _make_pair()
    rows = designer._format_primer_pair_v2(
        pp, snp_position=50, snp_alleles=("A", "G"), variant_sites=[], key="31-0"
    )
    assert rows[0].genomic_start is None
    assert rows[0].genomic_end is None


def test_format_pair_genomic_plus_strand(designer: KASPDesigner) -> None:
    """Plus strand: genomic = genomic_start + template position."""
    designer._genomic_start = 1000
    designer._genomic_strand = Strand.PLUS
    pp = _make_pair(left_start=10, left_end=30)
    rows = designer._format_primer_pair_v2(
        pp, snp_position=30, snp_alleles=("A", "G"), variant_sites=[], key="31-0"
    )
    # Allele row (LEFT direction)
    assert rows[0].genomic_start == 1010  # 1000 + 10
    assert rows[0].genomic_end == 1030  # 1000 + 30


def test_format_pair_genomic_minus_strand(designer: KASPDesigner) -> None:
    """Minus strand: coordinates reversed relative to template length."""
    designer._genomic_start = 2000
    designer._genomic_strand = Strand.MINUS
    designer._template_sequence = "a" * 200
    pp = _make_pair(left_start=10, left_end=30)
    rows = designer._format_primer_pair_v2(
        pp, snp_position=30, snp_alleles=("A", "G"), variant_sites=[], key="31-0"
    )
    # LEFT allele row: minus strand
    # genomic_start = 2000 + (200 - 30 - 1) = 2169
    # genomic_end = 2000 + (200 - 10 - 1) = 2189
    assert rows[0].genomic_start == 2169
    assert rows[0].genomic_end == 2189


# ── output fields completeness ───────────────────────────────────────


def test_format_pair_all_expected_fields_present(designer: KASPDesigner) -> None:
    """All V2 output fields are present in each row."""
    pp = _make_pair()
    rows = designer._format_primer_pair_v2(
        pp, snp_position=50, snp_alleles=("A", "G"), variant_sites=[], key="31-0"
    )
    expected_fields = {
        "index",
        "product_size",
        "direction",
        "start",
        "end",
        "genomic_start",
        "genomic_end",
        "diff_num",
        "diff_three_all",
        "length",
        "tm",
        "gc_percent",
        "self_any",
        "self_end",
        "end_stability",
        "hairpin",
        "primer_seq",
        "reverse_complement",
        "penalty",
        "compl_any",
        "compl_end",
        "score",
        "design_quality",
        "snp_name",
    }
    for row in rows:
        assert set(asdict(row).keys()) == expected_fields


def test_format_pair_penalty_propagated(designer: KASPDesigner) -> None:
    """Pair-level penalty/compl values propagate to all 3 rows."""
    pp = _make_pair()
    rows = designer._format_primer_pair_v2(
        pp, snp_position=50, snp_alleles=("A", "G"), variant_sites=[], key="31-0"
    )
    for row in rows:
        assert row.penalty == 1.5
        assert row.compl_any == 0.1
        assert row.compl_end == 0.2
