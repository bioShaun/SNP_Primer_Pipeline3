"""Unit tests for KASPDesigner._generate_primer3_inputs_v2.

These assertions pin the V2-compatible Primer3 input generation logic:
- variant site at the SNP is skipped,
- variant sites too far from the SNP (product too large) are skipped,
- force_left_end / force_right_end are emitted in 1-based coordinates with
  the left/right orientation that matches whether the variant site is
  upstream or downstream of the SNP,
- sequence_id uses 1-based var site index.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from snp_primer_pipeline.primers.kasp import KASPDesigner

pytestmark = pytest.mark.unit


@pytest.fixture
def designer(tmp_path: Path) -> KASPDesigner:
    """KASPDesigner with stubbed primer3 executable (not actually run here)."""
    fake = tmp_path / "primer3_core"
    fake.write_text("#!/bin/sh\necho\n")
    fake.chmod(0o755)
    return KASPDesigner(primer3_path=fake, config_path=None)


def _parse_settings(input_str: str) -> dict[str, list[str]]:
    """Parse a Primer3 input block into a dict of key -> list of values."""
    settings: dict[str, list[str]] = {}
    for line in input_str.splitlines():
        if not line or line == "=" or "=" not in line:
            continue
        key, _, value = line.partition("=")
        settings.setdefault(key, []).append(value)
    return settings


def test_generate_inputs_skips_snp_position(designer: KASPDesigner) -> None:
    template = "A" * 400
    snp = 200
    variant_sites = [snp]  # only the SNP itself -> skipped

    inputs = designer._generate_primer3_inputs_v2(
        template=template,
        snp_position=snp,
        variant_sites=variant_sites,
        product_size_range=(50, 250),
    )

    assert inputs == []


def test_generate_inputs_skips_sites_exceeding_product_max(designer: KASPDesigner) -> None:
    template = "A" * 1000
    snp = 500
    # product_max - 35 = 215 -> anything with |var - snp| > 215 must be skipped
    close_site = snp - 100
    far_site = snp + 300  # 300 > 215 -> skipped

    inputs = designer._generate_primer3_inputs_v2(
        template=template,
        snp_position=snp,
        variant_sites=[close_site, far_site],
        product_size_range=(50, 250),
    )

    assert len(inputs) == 1
    assert inputs[0][0] == close_site


def test_generate_inputs_var_before_snp_forces_left_var_right_snp(
    designer: KASPDesigner,
) -> None:
    template = "A" * 400
    snp = 200
    var_site = 150  # upstream

    inputs = designer._generate_primer3_inputs_v2(
        template=template,
        snp_position=snp,
        variant_sites=[var_site],
        product_size_range=(50, 250),
    )

    assert len(inputs) == 1
    var_returned, input_str = inputs[0]
    assert var_returned == var_site

    settings = _parse_settings(input_str)
    assert settings["SEQUENCE_ID"] == [f"var_{var_site + 1}"]
    assert settings["SEQUENCE_FORCE_LEFT_END"] == [str(var_site + 1)]
    assert settings["SEQUENCE_FORCE_RIGHT_END"] == [str(snp + 1)]
    assert settings["SEQUENCE_TEMPLATE"] == [template]


def test_generate_inputs_var_after_snp_forces_left_snp_right_var(
    designer: KASPDesigner,
) -> None:
    template = "A" * 400
    snp = 100
    var_site = 180  # downstream

    inputs = designer._generate_primer3_inputs_v2(
        template=template,
        snp_position=snp,
        variant_sites=[var_site],
        product_size_range=(50, 250),
    )

    assert len(inputs) == 1
    _, input_str = inputs[0]

    settings = _parse_settings(input_str)
    assert settings["SEQUENCE_FORCE_LEFT_END"] == [str(snp + 1)]
    assert settings["SEQUENCE_FORCE_RIGHT_END"] == [str(var_site + 1)]
    assert settings["SEQUENCE_ID"] == [f"var_{var_site + 1}"]


def test_generate_inputs_core_settings_are_v2_defaults(designer: KASPDesigner) -> None:
    inputs = designer._generate_primer3_inputs_v2(
        template="A" * 400,
        snp_position=200,
        variant_sites=[150],
        product_size_range=(50, 250),
    )
    assert len(inputs) == 1
    settings = _parse_settings(inputs[0][1])

    assert settings["PRIMER_MAX_SIZE"] == [str(designer.max_size)]
    assert settings["PRIMER_MIN_TM"] == ["57.0"]
    assert settings["PRIMER_OPT_TM"] == ["60.0"]
    assert settings["PRIMER_MAX_TM"] == [str(designer.max_tm)]
    assert settings["PRIMER_PAIR_MAX_DIFF_TM"] == ["6.0"]
    assert settings["PRIMER_FIRST_BASE_INDEX"] == ["1"]
    assert settings["PRIMER_LIBERAL_BASE"] == ["1"]
    assert settings["PRIMER_NUM_RETURN"] == ["5"]
    assert settings["PRIMER_EXPLAIN_FLAG"] == ["1"]
    # pick_anyway defaults to False -> 0
    assert settings["PRIMER_PICK_ANYWAY"] == ["0"]
    # Product size ranges are the three V2-fixed tiers regardless of input range.
    assert settings["PRIMER_PRODUCT_SIZE_RANGE"] == ["50-100 100-150 150-250"]


def test_generate_inputs_pick_anyway_propagates(tmp_path: Path) -> None:
    fake = tmp_path / "primer3_core"
    fake.write_text("#!/bin/sh\necho\n")
    fake.chmod(0o755)
    designer = KASPDesigner(primer3_path=fake, config_path=None, pick_anyway=True)

    inputs = designer._generate_primer3_inputs_v2(
        template="A" * 400,
        snp_position=200,
        variant_sites=[150],
        product_size_range=(50, 250),
    )
    settings = _parse_settings(inputs[0][1])
    assert settings["PRIMER_PICK_ANYWAY"] == ["1"]
