"""Unit tests for CAPSDesigner enzyme screening (no Primer3)."""

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


def test_reverse_complement_caps_lowercase() -> None:
    from snp_primer_pipeline.models import reverse_complement

    assert reverse_complement("gatc") == "gatc"


def test_test_enzyme_caps_different_cut_counts(designer: CAPSDesigner) -> None:
    """Wild has two sites, mutant has one -> CAPS."""
    enzyme = RestrictionEnzyme(name="EcoRI", sequence="GAATTC", price=50)
    wild = "gaattc" + "n" * 6 + "gaattc"
    mut = "gaattc" + "n" * 6 + "gaattt"
    out = designer._test_enzyme(enzyme, wild, mut, snp_position=6)
    assert out.caps is True
