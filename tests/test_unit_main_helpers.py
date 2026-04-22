"""Unit tests for main.py helper functions."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from snp_primer_pipeline.config import PipelineConfig
from snp_primer_pipeline.main import (
    _design_caps_primers,
    _design_kasp_primers,
    _parse_input_file,
    _prepare_output_dir,
    process_snp,
)
from snp_primer_pipeline.models import SNP, FlankingRegion, Strand

pytestmark = pytest.mark.unit


def _make_config(tmp_path: Path, **kwargs: object) -> PipelineConfig:
    """Create a minimal valid PipelineConfig for tests."""
    inp = tmp_path / "dummy.csv"
    inp.write_text("s1,chr1,ATCG[A/G]TCG\n")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nATCG\n")
    out = tmp_path / "out"
    out.mkdir()
    defaults = {
        "input_file": inp,
        "reference_file": ref,
        "output_dir": out,
    }
    defaults.update(kwargs)  # type: ignore[typeddict-item]
    return PipelineConfig(**defaults)  # type: ignore[arg-type]


def test_prepare_output_dir_creates_directory(tmp_path: Path) -> None:
    output_dir = tmp_path / "output"
    config = _make_config(tmp_path, output_dir=output_dir)
    _prepare_output_dir(config)
    assert output_dir.exists()
    assert output_dir.is_dir()


def test_parse_input_file_coordinates_format(tmp_path: Path) -> None:
    inp = tmp_path / "in.txt"
    inp.write_text("chr1\t100\tA\tG\n")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nATCG\n")
    out = tmp_path / "out"
    out.mkdir()

    config = PipelineConfig(
        input_file=inp,
        reference_file=ref,
        output_dir=out,
    )

    mock_snp = SNP(
        name="s1",
        chromosome="chr1",
        flanking_sequence="ATCG",
        snp_position=0,
        allele_a="A",
        allele_b="G",
        iupac_code="R",
    )

    with patch("snp_primer_pipeline.main.PolymarkerParser") as MockParser:
        MockParser.detect_format.return_value = "coordinates"
        mock_parser = MockParser.return_value
        mock_parser.parse_coordinates.return_value = [mock_snp]
        mock_parser.fetch_target_flanking.return_value = ({}, {})

        snps, input_format, target_regions_direct, target_sequences_direct, exclude_chromosomes = _parse_input_file(config)

        MockParser.assert_called_once_with(inp)
        mock_parser.parse_coordinates.assert_called_once_with(ref)
        assert len(snps) == 1
        assert input_format == "coordinates"
        assert exclude_chromosomes == {"s1": "chr1"}
        assert target_regions_direct == {}
        assert target_sequences_direct == {}


def test_parse_input_file_standard_format(tmp_path: Path) -> None:
    inp = tmp_path / "in.csv"
    inp.write_text("s1,chr1,ATCG[A/G]TCG\n")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nATCG\n")
    out = tmp_path / "out"
    out.mkdir()

    config = PipelineConfig(
        input_file=inp,
        reference_file=ref,
        output_dir=out,
    )

    mock_snp = SNP(
        name="s1",
        chromosome="chr1",
        flanking_sequence="ATCG[A/G]TCG",
        snp_position=4,
        allele_a="A",
        allele_b="G",
        iupac_code="R",
    )

    with patch("snp_primer_pipeline.main.PolymarkerParser") as MockParser:
        MockParser.detect_format.return_value = "polymarker"
        mock_parser = MockParser.return_value
        mock_parser.parse.return_value = [mock_snp]

        snps, input_format, target_regions_direct, target_sequences_direct, exclude_chromosomes = _parse_input_file(config)

        mock_parser.parse.assert_called_once()
        assert len(snps) == 1
        assert input_format == "polymarker"
        assert exclude_chromosomes is None
        assert target_sequences_direct == {}


def test_design_kasp_primers_runs_designer(tmp_path: Path) -> None:
    config = _make_config(tmp_path, run_specificity=False)
    snp = SNP(
        name="s1",
        chromosome="chr1",
        flanking_sequence="A" * 80,
        snp_position=24,
        allele_a="A",
        allele_b="G",
        iupac_code="R",
    )

    with patch("snp_primer_pipeline.main.KASPDesigner") as MockDesigner:
        mock_designer = MockDesigner.return_value
        mock_designer.design_primers.return_value = [
            {"primer": "left"},
            {"primer": "right"},
            {"primer": "common"},
        ]

        result = _design_kasp_primers(
            snp=snp,
            target_sequence="A" * 100,
            target_snp_position=24,
            genomic_start=1,
            genomic_strand="+",
            config=config,
            snp_output_dir=tmp_path / "out",
            alignment=None,
            sites_diff_all=[],
            diffarray={},
            target_name=None,
            flanking_regions=[],
            has_multiple_sequences=False,
        )

        MockDesigner.assert_called_once_with(
            max_tm=config.max_tm,
            max_size=config.max_primer_size,
            pick_anyway=config.pick_anyway,
        )
        mock_designer.design_primers.assert_called_once()
        assert result[0] is not None
        assert result[1] is None
        assert result[2] is None


def test_design_kasp_primers_with_specificity(tmp_path: Path) -> None:
    config = _make_config(tmp_path, run_specificity=True)
    snp = SNP(
        name="s1",
        chromosome="chr1",
        flanking_sequence="A" * 80,
        snp_position=24,
        allele_a="A",
        allele_b="G",
        iupac_code="R",
    )

    mock_region = FlankingRegion(
        snp_name="s1",
        chromosome="chr1",
        start=1,
        end=100,
        strand=Strand.PLUS,
        snp_position_in_region=25,
        allele="A",
    )

    with (
        patch("snp_primer_pipeline.main.KASPDesigner") as MockDesigner,
        patch("snp_primer_pipeline.main.SpecificityBlastRunner") as MockSpecRunner,
        patch("snp_primer_pipeline.main.SpecificityAssessor") as MockAssessor,
    ):
        mock_designer = MockDesigner.return_value
        mock_designer.design_primers.return_value = [
            {"primer": "left"},
            {"primer": "right"},
            {"primer": "common"},
        ]

        mock_spec_runner = MockSpecRunner.return_value
        mock_spec_runner.prepare_primer_fasta.return_value = tmp_path / "primers.fa"

        mock_assessor = MockAssessor.return_value
        mock_assessor.parse_blast_output.return_value = []
        mock_assessor.assess.return_value = {}
        mock_assessor.select_best.return_value = "set_1"

        result = _design_kasp_primers(
            snp=snp,
            target_sequence="A" * 100,
            target_snp_position=24,
            genomic_start=1,
            genomic_strand="+",
            config=config,
            snp_output_dir=tmp_path / "out",
            alignment=None,
            sites_diff_all=[],
            diffarray={},
            target_name=None,
            flanking_regions=[mock_region],
            has_multiple_sequences=False,
        )

        assert result[0] is not None
        MockSpecRunner.assert_called_once()
        mock_assessor.assess.assert_called_once()
        assert result[2] == "set_1"


def test_design_kasp_primers_without_specificity(tmp_path: Path) -> None:
    config = _make_config(tmp_path, run_specificity=False)
    snp = SNP(
        name="s1",
        chromosome="chr1",
        flanking_sequence="A" * 80,
        snp_position=24,
        allele_a="A",
        allele_b="G",
        iupac_code="R",
    )

    with (
        patch("snp_primer_pipeline.main.KASPDesigner") as MockDesigner,
        patch("snp_primer_pipeline.main.SpecificityBlastRunner") as MockSpecRunner,
    ):
        mock_designer = MockDesigner.return_value
        mock_designer.design_primers.return_value = [
            {"primer": "left"},
            {"primer": "right"},
            {"primer": "common"},
        ]

        result = _design_kasp_primers(
            snp=snp,
            target_sequence="A" * 100,
            target_snp_position=24,
            genomic_start=1,
            genomic_strand="+",
            config=config,
            snp_output_dir=tmp_path / "out",
            alignment=None,
            sites_diff_all=[],
            diffarray={},
            target_name=None,
            flanking_regions=[],
            has_multiple_sequences=False,
        )

        assert result[0] is not None
        MockSpecRunner.assert_not_called()


def test_design_caps_primers_runs_designer(tmp_path: Path) -> None:
    config = _make_config(tmp_path, design_caps=True)
    snp = SNP(
        name="s1",
        chromosome="chr1",
        flanking_sequence="A" * 80,
        snp_position=24,
        allele_a="A",
        allele_b="G",
        iupac_code="R",
    )

    with patch("snp_primer_pipeline.main.CAPSDesigner") as MockDesigner:
        mock_designer = MockDesigner.return_value
        mock_designer.find_usable_enzymes.return_value = ([], [])
        mock_designer.format_output.return_value = None

        result = _design_caps_primers(
            snp=snp,
            target_sequence="A" * 100,
            target_snp_position=24,
            config=config,
            snp_output_dir=tmp_path / "out",
            sites_diff_all=[],
        )

        MockDesigner.assert_called_once()
        mock_designer.find_usable_enzymes.assert_called_once()
        assert result == []


def test_process_snp_skips_design_when_disabled(tmp_path: Path) -> None:
    config = _make_config(tmp_path, design_kasp=False, design_caps=False, run_specificity=False)

    snp = SNP(
        name="s1",
        chromosome="chr1",
        flanking_sequence="A" * 80,
        snp_position=24,
        allele_a="A",
        allele_b="G",
        iupac_code="R",
    )

    region = FlankingRegion(
        snp_name="s1",
        chromosome="chr1",
        start=1,
        end=100,
        strand=Strand.PLUS,
        snp_position_in_region=25,
        allele="A",
    )
    seq = "A" * 100
    extracted = {
        f"{region.snp_name}_{region.chromosome}_{region.allele}_{region.snp_position_in_region}": seq
    }

    with (
        patch("snp_primer_pipeline.main.KASPDesigner") as MockKasp,
        patch("snp_primer_pipeline.main.CAPSDesigner") as MockCaps,
    ):
        result = process_snp(snp, [region], extracted, config)
        MockKasp.assert_not_called()
        MockCaps.assert_not_called()
        assert result["kasp"] == []
        assert result["caps"] == []


def test_process_snp_runs_caps_when_enabled(tmp_path: Path) -> None:
    config = _make_config(tmp_path, design_kasp=False, design_caps=True, run_specificity=False)

    snp = SNP(
        name="s1",
        chromosome="chr1",
        flanking_sequence="A" * 80,
        snp_position=24,
        allele_a="A",
        allele_b="G",
        iupac_code="R",
    )

    region = FlankingRegion(
        snp_name="s1",
        chromosome="chr1",
        start=1,
        end=100,
        strand=Strand.PLUS,
        snp_position_in_region=25,
        allele="A",
    )
    seq = "A" * 100
    extracted = {
        f"{region.snp_name}_{region.chromosome}_{region.allele}_{region.snp_position_in_region}": seq
    }

    with patch("snp_primer_pipeline.main.CAPSDesigner") as MockCaps:
        mock_designer = MockCaps.return_value
        mock_designer.find_usable_enzymes.return_value = ([], [])
        mock_designer.format_output.return_value = None

        result = process_snp(snp, [region], extracted, config)

        MockCaps.assert_called_once()
        assert result["caps"] == []
