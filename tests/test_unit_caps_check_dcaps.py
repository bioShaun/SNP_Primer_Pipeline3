"""Unit tests for CAPSDesigner._check_dcaps — full branch coverage (no Primer3)."""

from __future__ import annotations

from pathlib import Path

import pytest

from snp_primer_pipeline.models import RestrictionEnzyme
from snp_primer_pipeline.primers.caps import CAPSDesigner

pytestmark = pytest.mark.unit


@pytest.fixture
def designer(tmp_path: Path) -> CAPSDesigner:
    fake = tmp_path / "primer3_core"
    fake.write_text("#!/bin/sh\necho\n")
    fake.chmod(0o755)
    return CAPSDesigner(primer3_path=fake, enzyme_file=None)


# ── Branch 1: dCAPS found on forward strand ──────────────────────────


def test_dcaps_found_forward_strand(designer: CAPSDesigner) -> None:
    """SNP falls in enzyme match region, change_pos >1 away from SNP,
    mutation prevents cutting → dcaps=True on first (forward) pass."""
    # Enzyme: GAATTC (6 bases)
    # We mask position 0 → pattern = NAATTC → [ACGT]AATTC
    # Wild:  ...GAATTC... where SNP is within the match region
    # Mut:   ...GCATTC... — mutant changes a base so the full enzyme no longer matches
    #
    # Build sequences where SNP is at position 5 (inside GAATTC at 0..5).
    # Enzyme position i=0 → change_pos=0, abs(5-0)=5 > 1 ✓
    # Mut region "gcattc" does NOT match [acgt]aattc → dcaps found.
    wild = "gaattc" + "a" * 10
    mut = "gcattc" + "a" * 10
    snp_pos = 5  # inside the match (0..5)

    enzyme = RestrictionEnzyme(name="EcoRI", sequence="GAATTC", price=50)
    result = designer._check_dcaps(enzyme, wild, mut, snp_pos)

    assert result.dcaps is True
    assert result.change_pos is not None
    assert result.template_seq != ""
    assert len(result.primer_end_positions) > 0


# ── Branch 2: change_pos too close to SNP (abs <= 1) → skipped ───────


def test_dcaps_skips_when_change_pos_adjacent_to_snp(designer: CAPSDesigner) -> None:
    """If the only viable variable position in the enzyme is ≤1 base from SNP, skip it."""
    # Enzyme: GC (2 bases, palindromic)
    # Place enzyme at 0..1, SNP at position 0.
    # i=0: change_pos=0, abs(0-0)=0 ≤ 1 → skip.
    # i=1: change_pos=1, abs(0-1)=1 ≤ 1 → skip.
    # GC is palindromic → no recursion.
    wild = "gc" + "a" * 10
    mut = "tc" + "a" * 10
    snp_pos = 0

    enzyme = RestrictionEnzyme(name="Short", sequence="GC", price=10)
    result = designer._check_dcaps(enzyme, wild, mut, snp_pos)
    assert result.dcaps is False


# ── Branch 3: SNP outside match region → match skipped ───────────────


def test_dcaps_skips_when_snp_outside_match(designer: CAPSDesigner) -> None:
    """Enzyme match exists but SNP is far away → no dCAPS."""
    # Enzyme at positions 0..5, SNP at position 20
    wild = "gaattc" + "a" * 20
    mut = "gaattc" + "a" * 19 + "t"
    snp_pos = 20

    # RC of GAATTC is GAATTC (palindromic) → no recursion
    enzyme = RestrictionEnzyme(name="EcoRI", sequence="GAATTC", price=50)
    result = designer._check_dcaps(enzyme, wild, mut, snp_pos)
    assert result.dcaps is False


# ── Branch 4: mutation doesn't prevent cutting → skipped ─────────────


def test_dcaps_skips_when_mutation_still_matches(designer: CAPSDesigner) -> None:
    """Even though SNP is inside match & change_pos is far enough,
    if the mutant region still matches the degenerate pattern → no dCAPS."""
    # Enzyme: RGATCY (R=[AG], Y=[CT]) - 6 bases, pattern [ag]gatc[ct]
    # Wild: agatcc at 0..5, SNP at pos 3 ('t')
    # Mut:  agaccc - SNP at pos 3: t->c
    # For every i where abs(snp-change_pos)>1 and SNP in match:
    #   i=0: pattern [acgt]gatc[ct], matches wild "agatcc"; change_pos=0, abs(3-0)=3>1
    #         mut_region "agaccc" vs [acgt]gatc[ct] → 'c' at pos 3 ≠ 't' → no match → would find dcaps
    # We need mut to still match. Let's use SNP at pos 5 (the Y position) with
    # mut value that still matches Y=[CT]:
    # Wild: agatcc (pos 5 = 'c'), Mut: agatct (pos 5 = 't'), both match Y.
    # i=0: pattern [acgt]gatc[ct], change_pos=0, abs(5-0)=5>1
    #   mut_region "agatct" vs [acgt]gatc[ct] → matches! → skip.
    # All other i values where SNP is in match and abs>1 similarly match → no dcaps.
    wild = "agatcc" + "a" * 10
    mut = "agatct" + "a" * 10
    snp_pos = 5

    enzyme = RestrictionEnzyme(name="BamHI_var", sequence="RGATCY", price=50)
    result = designer._check_dcaps(enzyme, wild, mut, snp_pos)
    assert result.dcaps is False


# ── Branch 5: forward finds nothing → recurses with RC ───────────────


def test_dcaps_found_on_reverse_complement_strand(designer: CAPSDesigner) -> None:
    """Forward enzyme doesn't match, but its RC does → dCAPS via recursion."""
    # Use a non-palindromic enzyme so RC ≠ forward.
    # Enzyme: GGATC (5 bases), RC = GATCC
    # Place RC match in wild: gatcc at 0..4, SNP at pos 4
    # Forward GGATC won't match "gatcc" → falls through to RC branch.
    # RC "gatcc" with i=0 → N=NATCC pattern, change_pos=0, abs(4-0)=4>1
    # mut changes pos 4: gatct → mut_region "gatct" vs [acgt]atcc → no match → dcaps!
    wild = "gatcc" + "a" * 10
    mut = "gatct" + "a" * 10
    snp_pos = 4

    enzyme = RestrictionEnzyme(name="AsymEnz", sequence="GGATC", price=30)
    result = designer._check_dcaps(enzyme, wild, mut, snp_pos)
    assert result.dcaps is True


# ── Branch 6: palindromic enzyme, no match → no recursion, returns ───


def test_dcaps_palindromic_no_match_no_recursion(designer: CAPSDesigner) -> None:
    """Palindromic enzyme (RC==self), no match anywhere → returns immediately."""
    wild = "t" * 20
    mut = "t" * 20
    snp_pos = 10

    # GAATTC is palindromic
    enzyme = RestrictionEnzyme(name="EcoRI", sequence="GAATTC", price=50)
    result = designer._check_dcaps(enzyme, wild, mut, snp_pos)
    assert result.dcaps is False
    assert result.template_seq == ""


# ── Branch 7: primer_end_positions direction (change_pos > snp) ──────


def test_dcaps_primer_end_positions_change_after_snp(designer: CAPSDesigner) -> None:
    """When change_pos > snp_position, primer_end_positions = range(snp+1, change_pos)."""
    # Enzyme: GAATTC (6 bases)
    # Place match at 0..5, SNP at position 0.
    # Iteration order: i=0 first. i=0 → change_pos=0, abs(0-0)=0 ≤ 1 → skip.
    # i=1 → change_pos=1, abs(0-1)=1 ≤ 1 → skip.
    # i=2 → pattern GA[ACGT]TTC → change_pos=2, abs(0-2)=2>1 → check mut.
    #   mut = "taattc" → mut_region "taattc" vs ga[acgt]ttc → 't' ≠ 'g' at pos 0 → no match → dcaps!
    #   change_pos(2) > snp_pos(0) → primer_end_positions = range(1, 2) = [1]
    wild = "gaattc" + "a" * 10
    mut = "taattc" + "a" * 10  # SNP at pos 0: g→t
    snp_pos = 0

    enzyme = RestrictionEnzyme(name="EcoRI", sequence="GAATTC", price=50)
    result = designer._check_dcaps(enzyme, wild, mut, snp_pos)

    assert result.dcaps is True
    assert result.change_pos == 3  # 1-based: change_pos=2 → stored as 2+1=3
    assert result.primer_end_positions == [1]  # range(0+1, 2) = [1]


# ── Integration: _test_enzyme delegates to _check_dcaps correctly ────


def test_test_enzyme_delegates_to_check_dcaps(designer: CAPSDesigner) -> None:
    """_test_enzyme falls through to _check_dcaps when cut counts are equal."""
    # Use sequences where enzyme has 0 exact cuts in both wild and mut,
    # so _test_enzyme delegates to _check_dcaps.
    # GAATTC doesn't match "aaattc" exactly → 0 cuts each → equal → _check_dcaps.
    # _check_dcaps: i=0 → pattern [acgt]aattc matches "aaattc", SNP at 1 inside 0..5,
    # change_pos=0, abs(1-0)=1 <= 1 → skip.
    # Other i values: patterns don't match "aaattc" (need 'g' at pos 0 etc).
    # GAATTC is palindromic → no RC recursion.
    wild = "aaattc" + "a" * 10
    mut = "acattc" + "a" * 10
    snp_pos = 1

    enzyme = RestrictionEnzyme(name="EcoRI", sequence="GAATTC", price=50)
    result = designer._test_enzyme(enzyme, wild, mut, snp_pos)
    assert result.dcaps is False
    assert result.caps is False
