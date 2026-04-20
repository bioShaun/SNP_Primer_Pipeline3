"""Unit tests for KASPDesigner.format_output / format_simple_output.

These lock down the exact tab-separated layout of KASP output files so that
future refactors cannot silently change column order, formatting (e.g. Tm
precision), or specificity/best-primer integration without failing a test.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pytest

from snp_primer_pipeline.core.specificity import SpecificityResult, SpecificityStatus
from snp_primer_pipeline.primers.kasp import KASPDesigner

if TYPE_CHECKING:
    from typing import Any

pytestmark = pytest.mark.unit


FULL_HEADER = [
    "index",
    "product_size",
    "type",
    "start",
    "end",
    "genomic_start",
    "genomic_end",
    "variation number",
    "3'diffall",
    "length",
    "Tm",
    "GCcontent",
    "any",
    "3'",
    "end_stability",
    "hairpin",
    "primer_seq",
    "ReverseComplement",
    "penalty",
    "compl_any",
    "compl_end",
    "score",
    "specificity",
    "design_quality",
]

SIMPLE_HEADER = [
    "Index",
    "Allele_A",
    "Tm_A",
    "GC_A",
    "Allele_B",
    "Tm_B",
    "GC_B",
    "Common",
    "Tm_C",
    "GC_C",
    "Product_Size",
    "Genomic_Range",
    "Score",
    "Specificity",
    "Design_Quality",
    "Best_Primer",
]


@pytest.fixture
def designer(tmp_path: Path) -> KASPDesigner:
    """KASPDesigner that does not require a real primer3 executable."""
    fake = tmp_path / "primer3_core"
    fake.write_text("#!/bin/sh\necho\n")
    fake.chmod(0o755)
    return KASPDesigner(primer3_path=fake, config_path=None)


def _make_result(
    *,
    index: str,
    direction: str,
    primer_seq: str = "ACGTACGTac",
    reverse_complement: str = "gtACGTACGT",
    diff_num: int = 2,
    diff_three_all: str = "YES",
    tm: float = 60.123,
    gc_percent: float = 45.678,
    start: int = 10,
    end: int = 20,
    genomic_start: int | None = 1000,
    genomic_end: int | None = 1010,
    product_size: int = 150,
    design_quality: str | None = "STRICT",
    score: float = 1.2345,
) -> dict[str, Any]:
    """Build a kasp_result dict matching the shape emitted by design_primers."""
    result: dict[str, Any] = {
        "index": index,
        "product_size": product_size,
        "direction": direction,
        "start": start,
        "end": end,
        "genomic_start": genomic_start,
        "genomic_end": genomic_end,
        "diff_num": diff_num,
        "diff_three_all": diff_three_all,
        "length": len(primer_seq),
        "tm": tm,
        "gc_percent": gc_percent,
        "self_any": 2.3,
        "self_end": 1.0,
        "end_stability": 4.56,
        "hairpin": 0.0,
        "primer_seq": primer_seq,
        "reverse_complement": reverse_complement,
        "penalty": 0.75,
        "compl_any": 3.25,
        "compl_end": 1.5,
        "score": score,
    }
    if design_quality is not None:
        result["design_quality"] = design_quality
    return result


# ---------------------------------------------------------------------------
# format_output
# ---------------------------------------------------------------------------


def test_format_output_writes_header_and_single_row(designer: KASPDesigner, tmp_path: Path) -> None:
    out = tmp_path / "kasp.tsv"
    result = _make_result(index="1-Allele-A", direction="LEFT")

    designer.format_output([result], snp_name="SNP1", output_file=out)

    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == FULL_HEADER
    assert len(lines) == 2

    row = lines[1].split("\t")
    assert len(row) == len(FULL_HEADER)
    assert row[FULL_HEADER.index("index")] == "SNP1-1-Allele-A"
    assert row[FULL_HEADER.index("product_size")] == "150"
    assert row[FULL_HEADER.index("type")] == "LEFT"
    assert row[FULL_HEADER.index("genomic_start")] == "1000"
    assert row[FULL_HEADER.index("genomic_end")] == "1010"
    assert row[FULL_HEADER.index("variation number")] == "2"
    assert row[FULL_HEADER.index("3'diffall")] == "YES"
    assert row[FULL_HEADER.index("Tm")] == "60.12"
    assert row[FULL_HEADER.index("GCcontent")] == "45.68"
    assert row[FULL_HEADER.index("primer_seq")] == "ACGTACGTac"
    assert row[FULL_HEADER.index("ReverseComplement")] == "gtACGTACGT"
    assert row[FULL_HEADER.index("penalty")] == "0.75"
    assert row[FULL_HEADER.index("score")] == "1.23"
    assert row[FULL_HEADER.index("specificity")] == "NA"
    assert row[FULL_HEADER.index("design_quality")] == "STRICT"


def test_format_output_missing_genomic_coords_yields_na(
    designer: KASPDesigner, tmp_path: Path
) -> None:
    out = tmp_path / "kasp.tsv"
    result = _make_result(
        index="1-Allele-A",
        direction="LEFT",
        genomic_start=None,
        genomic_end=None,
    )

    designer.format_output([result], snp_name="SNP1", output_file=out)

    row = out.read_text().splitlines()[1].split("\t")
    assert row[FULL_HEADER.index("genomic_start")] == "NA"
    assert row[FULL_HEADER.index("genomic_end")] == "NA"


def test_format_output_missing_design_quality_yields_na(
    designer: KASPDesigner, tmp_path: Path
) -> None:
    out = tmp_path / "kasp.tsv"
    result = _make_result(index="1-Allele-A", direction="LEFT", design_quality=None)

    designer.format_output([result], snp_name="SNP1", output_file=out)

    row = out.read_text().splitlines()[1].split("\t")
    assert row[FULL_HEADER.index("design_quality")] == "NA"


def test_format_output_uses_specificity_by_base_key(designer: KASPDesigner, tmp_path: Path) -> None:
    out = tmp_path / "kasp.tsv"
    # Triplet per primer set: Allele-A / Allele-B / Common
    a = _make_result(index="1-Allele-A", direction="LEFT")
    b = _make_result(index="1-Allele-B", direction="LEFT")
    c = _make_result(index="1-Common", direction="RIGHT")

    spec_results = {
        "1": SpecificityResult(status=SpecificityStatus.WARNING, reason="weak off-target"),
    }

    designer.format_output(
        [a, b, c],
        snp_name="SNP1",
        output_file=out,
        specificity_results=spec_results,
    )

    lines = out.read_text().splitlines()
    assert len(lines) == 4  # header + 3 rows
    for row_line in lines[1:]:
        row = row_line.split("\t")
        assert row[FULL_HEADER.index("specificity")] == "WARNING"


def test_format_output_variant_sites_block_emitted_only_when_requested(
    designer: KASPDesigner, tmp_path: Path
) -> None:
    out = tmp_path / "kasp.tsv"
    result = _make_result(index="1-Allele-A", direction="LEFT")

    # variant_sites but show_variant_sites=False => no trailing block
    designer.format_output(
        [result],
        snp_name="SNP1",
        output_file=out,
        variant_sites=[4, 9, 14],
        show_variant_sites=False,
    )
    content_hidden = out.read_text()
    assert "Sites that can differ all" not in content_hidden

    # show_variant_sites=True => block appears with 1-based indices
    designer.format_output(
        [result],
        snp_name="SNP1",
        output_file=out,
        variant_sites=[4, 9, 14],
        show_variant_sites=True,
    )
    content_shown = out.read_text()
    assert "Sites that can differ all for SNP1" in content_shown
    assert "5, 10, 15" in content_shown


# ---------------------------------------------------------------------------
# format_simple_output
# ---------------------------------------------------------------------------


def _triplet(base_index: str) -> list[dict[str, Any]]:
    """Build a complete (Allele-A, Allele-B, Common) triplet for a set."""
    a = _make_result(
        index=f"{base_index}-Allele-A",
        direction="LEFT",
        primer_seq="AAAA",
        tm=60.1,
        gc_percent=50.0,
        genomic_start=1000,
        genomic_end=1010,
        score=2.34,
        design_quality="STRICT",
    )
    b = _make_result(
        index=f"{base_index}-Allele-B",
        direction="LEFT",
        primer_seq="TTTT",
        tm=60.5,
        gc_percent=50.1,
        genomic_start=1000,
        genomic_end=1010,
    )
    c = _make_result(
        index=f"{base_index}-Common",
        direction="RIGHT",
        primer_seq="GGGG",
        tm=59.9,
        gc_percent=55.0,
        genomic_start=1100,
        genomic_end=1150,
    )
    return [a, b, c]


def test_format_simple_output_single_set_layout(designer: KASPDesigner, tmp_path: Path) -> None:
    out = tmp_path / "simple.tsv"

    designer.format_simple_output(_triplet("1"), snp_name="SNP1", output_file=out)

    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == SIMPLE_HEADER
    assert len(lines) == 2

    row = lines[1].split("\t")
    assert len(row) == len(SIMPLE_HEADER)
    assert row[SIMPLE_HEADER.index("Index")] == "SNP1-1"
    assert row[SIMPLE_HEADER.index("Allele_A")] == "AAAA"
    assert row[SIMPLE_HEADER.index("Tm_A")] == "60.1"
    assert row[SIMPLE_HEADER.index("GC_A")] == "50.0"
    assert row[SIMPLE_HEADER.index("Allele_B")] == "TTTT"
    assert row[SIMPLE_HEADER.index("Common")] == "GGGG"
    # Genomic range is built from allele_a.genomic_start and common.genomic_end.
    assert row[SIMPLE_HEADER.index("Genomic_Range")] == "1000-1150"
    assert row[SIMPLE_HEADER.index("Score")] == "2.34"
    assert row[SIMPLE_HEADER.index("Specificity")] == "NA"
    assert row[SIMPLE_HEADER.index("Design_Quality")] == "STRICT"
    assert row[SIMPLE_HEADER.index("Best_Primer")] == "NO"


def test_format_simple_output_missing_genomic_range_yields_na(
    designer: KASPDesigner, tmp_path: Path
) -> None:
    out = tmp_path / "simple.tsv"
    triplet = _triplet("1")
    triplet[0]["genomic_start"] = None  # Allele-A has no genomic coords

    designer.format_simple_output(triplet, snp_name="SNP1", output_file=out)

    row = out.read_text().splitlines()[1].split("\t")
    assert row[SIMPLE_HEADER.index("Genomic_Range")] == "NA"


def test_format_simple_output_incomplete_last_set_is_skipped(
    designer: KASPDesigner, tmp_path: Path
) -> None:
    out = tmp_path / "simple.tsv"
    # 1 complete set + 2 dangling results (not a triplet) -> only 1 data row.
    results = _triplet("1") + _triplet("2")[:2]

    designer.format_simple_output(results, snp_name="SNP1", output_file=out)

    lines = out.read_text().splitlines()
    assert len(lines) == 2
    assert lines[1].split("\t")[SIMPLE_HEADER.index("Index")] == "SNP1-1"


def test_format_simple_output_best_primer_key_string(
    designer: KASPDesigner, tmp_path: Path
) -> None:
    out = tmp_path / "simple.tsv"
    results = _triplet("1") + _triplet("2")

    designer.format_simple_output(
        results,
        snp_name="SNP1",
        output_file=out,
        best_primer_key="2",
    )

    rows = [line.split("\t") for line in out.read_text().splitlines()[1:]]
    assert rows[0][SIMPLE_HEADER.index("Best_Primer")] == "NO"
    assert rows[1][SIMPLE_HEADER.index("Best_Primer")] == "YES"


def test_format_simple_output_best_primer_key_dict(designer: KASPDesigner, tmp_path: Path) -> None:
    out = tmp_path / "simple.tsv"
    results = _triplet("1") + _triplet("2")

    designer.format_simple_output(
        results,
        snp_name="SNP1",
        output_file=out,
        best_primer_key={"SNP1": "1"},
    )

    rows = [line.split("\t") for line in out.read_text().splitlines()[1:]]
    assert rows[0][SIMPLE_HEADER.index("Best_Primer")] == "YES"
    assert rows[1][SIMPLE_HEADER.index("Best_Primer")] == "NO"


def test_format_simple_output_specificity_and_quality(
    designer: KASPDesigner, tmp_path: Path
) -> None:
    out = tmp_path / "simple.tsv"
    results = _triplet("1")

    spec_results = {
        "1": SpecificityResult(status=SpecificityStatus.FAIL, reason="off-target"),
    }

    designer.format_simple_output(
        results,
        snp_name="SNP1",
        output_file=out,
        specificity_results=spec_results,
    )

    row = out.read_text().splitlines()[1].split("\t")
    assert row[SIMPLE_HEADER.index("Specificity")] == "FAIL"
    assert row[SIMPLE_HEADER.index("Design_Quality")] == "STRICT"
