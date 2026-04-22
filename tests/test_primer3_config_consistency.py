"""Verify KASP and CAPS designers use the same Primer3 configuration strategy."""

from __future__ import annotations

from pathlib import Path

import pytest

pytestmark = pytest.mark.unit


def test_both_designers_pass_no_settings_file(tmp_path: Path) -> None:
    """Both designers should create Primer3Runner with settings_file=None."""
    fake_p3 = tmp_path / "primer3_core"
    fake_p3.write_text("#!/bin/sh\necho\n")
    fake_p3.chmod(0o755)

    from snp_primer_pipeline.primers.caps import CAPSDesigner
    from snp_primer_pipeline.primers.kasp import KASPDesigner

    kasp = KASPDesigner(primer3_path=fake_p3)
    caps = CAPSDesigner(primer3_path=fake_p3, enzyme_file=None)

    assert kasp.primer3_runner.settings_file is None
    assert caps.primer3_runner.settings_file is None


def test_both_designers_auto_detect_primer3_path(tmp_path: Path) -> None:
    """Both designers should auto-detect primer3 if not explicitly provided."""
    from snp_primer_pipeline.primers.caps import CAPSDesigner
    from snp_primer_pipeline.primers.kasp import KASPDesigner

    kasp = KASPDesigner()
    caps = CAPSDesigner(enzyme_file=None)

    # Both should resolve to an actual path (not the default Path("primer3_core"))
    assert kasp.primer3_runner.primer3_path != Path("primer3_core")
    assert caps.primer3_runner.primer3_path != Path("primer3_core")
