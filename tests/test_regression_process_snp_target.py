"""Regression / xfail tests for main.process_snp target-sequence selection."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from snp_primer_pipeline.config import PipelineConfig
from snp_primer_pipeline.main import process_snp
from snp_primer_pipeline.models import SNP, FlankingRegion, Strand

pytestmark = pytest.mark.unit


@pytest.mark.xfail(
    reason="Issue #1: target sequence must match SNP chromosome, not first hit order",
    strict=False,
)
def test_process_snp_prefers_sequence_on_snp_chromosome(tmp_path: Path) -> None:
    inp = tmp_path / "in.csv"
    inp.write_text("x,chrT,ATCG[A/G]TCG\n")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chrT\nATCG\n")
    out = tmp_path / "out"
    out.mkdir()

    config = PipelineConfig(
        input_file=inp,
        reference_file=ref,
        output_dir=out,
        design_kasp=True,
        design_caps=False,
        run_specificity=False,
    )

    snp = SNP(
        name="s1",
        chromosome="chrTarget",
        flanking_sequence="A" * 80,
        snp_position=24,
        allele_a="A",
        allele_b="G",
        iupac_code="R",
    )

    r_wrong = FlankingRegion(
        snp_name="s1",
        chromosome="chrOther",
        start=1,
        end=100,
        strand=Strand.PLUS,
        snp_position_in_region=25,
        allele="A",
    )
    r_target = FlankingRegion(
        snp_name="s1",
        chromosome="chrTarget",
        start=1,
        end=100,
        strand=Strand.PLUS,
        snp_position_in_region=30,
        allele="A",
    )
    seq_wrong = "C" * 100
    seq_target = "G" * 100

    def sid(r: FlankingRegion) -> str:
        return f"{r.snp_name}_{r.chromosome}_{r.allele}_{r.snp_position_in_region}"

    extracted = {sid(r_wrong): seq_wrong, sid(r_target): seq_target}
    captured: list[str] = []

    def capture_design(*args: object, **_kwargs: object) -> list:
        captured.append(args[0])
        return []

    with patch("snp_primer_pipeline.main.KASPDesigner") as kd:
        kd.return_value.design_primers.side_effect = capture_design
        process_snp(snp, [r_wrong, r_target], extracted, config)

    assert captured
    assert captured[0] == seq_target
