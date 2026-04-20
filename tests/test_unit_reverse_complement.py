"""Tests for three reverse-complement implementations (documented behavioral differences)."""

from __future__ import annotations

from pathlib import Path

import pytest

from snp_primer_pipeline.models import Primer
from snp_primer_pipeline.primers.caps import CAPSDesigner
from snp_primer_pipeline.primers.kasp import KASPDesigner

pytestmark = pytest.mark.unit


@pytest.fixture
def kasp_designer(tmp_path: Path) -> KASPDesigner:
    fake = tmp_path / "primer3_core"
    fake.write_text("#!/bin/sh\necho\n")
    fake.chmod(0o755)
    return KASPDesigner(primer3_path=fake, config_path=None)


@pytest.fixture
def caps_designer(tmp_path: Path) -> CAPSDesigner:
    fake = tmp_path / "primer3_core"
    fake.write_text("#!/bin/sh\necho\n")
    fake.chmod(0o755)
    return CAPSDesigner(primer3_path=fake, enzyme_file=None)


def test_primer_reverse_complement_atgc() -> None:
    p = Primer(sequence="ATCG")
    assert p.reverse_complement() == "CGAT"


def test_kasp_reverse_complement_passes_unknown_bases(kasp_designer: KASPDesigner) -> None:
    assert kasp_designer._reverse_complement("Nn") == "nN"


def test_caps_vs_kasp_iupac_mismatch_documented(
    kasp_designer: KASPDesigner, caps_designer: CAPSDesigner
) -> None:
    """KASP maps only ATGC; CAPS maps extended IUPAC (merge should unify these)."""
    seq = "R"
    assert kasp_designer._reverse_complement(seq) == "R"
    assert caps_designer._reverse_complement(seq) == "Y"
