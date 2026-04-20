"""Unit tests for KASPDesigner._filter_primers_v2 (no Primer3)."""

from __future__ import annotations

from pathlib import Path

import pytest

from snp_primer_pipeline.models import Primer, PrimerPair
from snp_primer_pipeline.primers.kasp import KASPDesigner

pytestmark = pytest.mark.unit


@pytest.fixture
def designer(tmp_path: Path) -> KASPDesigner:
    """Avoid SoftwarePaths.auto_detect in environments without primer3."""
    fake = tmp_path / "primer3_core"
    fake.write_text("#!/bin/sh\necho\n")
    fake.chmod(0o755)
    return KASPDesigner(primer3_path=fake, config_path=None)


def test_filter_primers_v2_empty_diffarray_keeps_all(designer: KASPDesigner) -> None:
    left = Primer(
        name="L",
        sequence="A" * 10,
        start=0,
        end=35,
        length=36,
        direction="LEFT",
        diff_num=2,
    )
    right = Primer(
        name="R",
        sequence="T" * 10,
        start=150,
        end=200,
        length=51,
        direction="RIGHT",
        diff_num=2,
    )
    pp = PrimerPair(
        left=left,
        right=right,
        product_size=120,
        penalty=1.0,
        compl_any=0.0,
        compl_end=0.0,
        score=0.0,
    )
    raw = {"31-0": (pp, 30)}
    designer.diffarray = {}
    out = designer._filter_primers_v2(raw, snp_position=100, variant_sites=[30], template_length=250)
    assert "31-0" in out


def test_filter_primers_v2_diffarray_filters_when_min_sums_zero(designer: KASPDesigner) -> None:
    left = Primer(
        name="L",
        sequence="A" * 10,
        start=0,
        end=35,
        length=36,
        direction="LEFT",
        diff_num=2,
    )
    right = Primer(
        name="R",
        sequence="T" * 10,
        start=150,
        end=200,
        length=51,
        direction="RIGHT",
        diff_num=2,
    )
    pp = PrimerPair(
        left=left,
        right=right,
        product_size=120,
        penalty=1.0,
        compl_any=0.0,
        compl_end=0.0,
        score=0.0,
    )
    raw = {"31-0": (pp, 30)}
    # Window around left 3' (end=35): rr uses range(max(35-9, 20), 36) -> 26..35
    designer.diffarray = {k: [0] for k in range(26, 36)}
    out = designer._filter_primers_v2(raw, snp_position=100, variant_sites=[30], template_length=250)
    assert "31-0" not in out


def test_filter_primers_v2_diffarray_keeps_when_all_homeologs_differ(designer: KASPDesigner) -> None:
    left = Primer(
        name="L",
        sequence="A" * 10,
        start=0,
        end=35,
        length=36,
        direction="LEFT",
        diff_num=2,
    )
    right = Primer(
        name="R",
        sequence="T" * 10,
        start=150,
        end=200,
        length=51,
        direction="RIGHT",
        diff_num=2,
    )
    pp = PrimerPair(
        left=left,
        right=right,
        product_size=120,
        penalty=1.0,
        compl_any=0.0,
        compl_end=0.0,
        score=0.0,
    )
    raw = {"31-0": (pp, 30)}
    designer.diffarray = {k: [1] for k in range(26, 36)}
    out = designer._filter_primers_v2(raw, snp_position=100, variant_sites=[30], template_length=250)
    assert "31-0" in out
