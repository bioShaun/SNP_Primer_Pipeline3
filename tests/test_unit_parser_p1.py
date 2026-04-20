"""P1 unit tests for PolymarkerParser (format detection, coordinates path mocked)."""

from __future__ import annotations

from pathlib import Path

import pytest

from snp_primer_pipeline.core.parser import PolymarkerParser
from snp_primer_pipeline.exceptions import ParseError

pytestmark = pytest.mark.unit


def test_detect_format_polymarker(tmp_path: Path) -> None:
    p = tmp_path / "p.csv"
    p.write_text("SNP1,chr1A,ATCGATCG[A/G]TCGATCGA\n")
    assert PolymarkerParser.detect_format(p) == "polymarker"


def test_detect_format_coordinates(tmp_path: Path) -> None:
    p = tmp_path / "c.txt"
    p.write_text("chr1A\t100\tA\tG\n")
    assert PolymarkerParser.detect_format(p) == "coordinates"


def test_detect_format_unknown_raises(tmp_path: Path) -> None:
    p = tmp_path / "x.txt"
    p.write_text("only_one_field\n")
    with pytest.raises(ParseError, match="Unknown input format"):
        PolymarkerParser.detect_format(p)


def test_parse_coordinates_with_mocked_fetch(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    inp = tmp_path / "coords.csv"
    inp.write_text("chr1,100,A,G\n")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\n" + "A" * 500 + "\n")

    def fake_fetch(
        _self: PolymarkerParser,
        coordinates: list,
        _reference: Path,
        flank: int = 50,
    ) -> dict[str, str]:
        del flank  # signature matches; return is fixed-length for this test
        return {c["name"]: "A" * 101 for c in coordinates}

    monkeypatch.setattr(PolymarkerParser, "_fetch_sequences", fake_fetch)

    parser = PolymarkerParser(inp)
    snps = parser.parse_coordinates(ref)
    assert len(snps) == 1
    assert snps[0].name == "chr1-100"
    assert snps[0].chromosome == "chr1"
    assert snps[0].snp_position == 50
