"""Unit tests for the shared FASTA parser."""

from __future__ import annotations

from pathlib import Path

import pytest

from snp_primer_pipeline.core.fasta import parse_fasta, parse_fasta_to_dict

pytestmark = pytest.mark.unit


class TestParseFasta:
    def test_single_sequence(self, tmp_path: Path) -> None:
        fasta = tmp_path / "single.fa"
        fasta.write_text(
            ">seq1\nATCG\nATCG\n",
            encoding="utf-8",
        )
        result = list(parse_fasta(fasta))
        assert result == [("seq1", "ATCGATCG")]

    def test_multiple_sequences(self, tmp_path: Path) -> None:
        fasta = tmp_path / "multi.fa"
        fasta.write_text(
            ">seq1\nAAAA\n>seq2\nCCCC\n>seq3\nGGGG\n",
            encoding="utf-8",
        )
        result = list(parse_fasta(fasta))
        assert result == [
            ("seq1", "AAAA"),
            ("seq2", "CCCC"),
            ("seq3", "GGGG"),
        ]

    def test_empty_lines_ignored(self, tmp_path: Path) -> None:
        fasta = tmp_path / "empty_lines.fa"
        fasta.write_text(
            ">seq1\n\nATCG\n\n\n>seq2\n\nGGGG\n\n",
            encoding="utf-8",
        )
        result = list(parse_fasta(fasta))
        assert result == [
            ("seq1", "ATCG"),
            ("seq2", "GGGG"),
        ]

    def test_empty_file(self, tmp_path: Path) -> None:
        fasta = tmp_path / "empty.fa"
        fasta.write_text("", encoding="utf-8")
        result = list(parse_fasta(fasta))
        assert result == []

    def test_no_sequences(self, tmp_path: Path) -> None:
        fasta = tmp_path / "no_seqs.fa"
        fasta.write_text("\n\n\n", encoding="utf-8")
        result = list(parse_fasta(fasta))
        assert result == []

    def test_id_with_whitespace_preserved(self, tmp_path: Path) -> None:
        fasta = tmp_path / "whitespace_id.fa"
        fasta.write_text(
            ">seq1 extra info\nATCG\n",
            encoding="utf-8",
        )
        result = list(parse_fasta(fasta))
        assert result == [("seq1 extra info", "ATCG")]

    def test_file_not_found_raises_oserror(self, tmp_path: Path) -> None:
        missing = tmp_path / "nonexistent.fa"
        with pytest.raises(OSError):
            list(parse_fasta(missing))


class TestParseFastaToDict:
    def test_multiple_sequences(self, tmp_path: Path) -> None:
        fasta = tmp_path / "multi.fa"
        fasta.write_text(
            ">seq1\nAAAA\n>seq2\nCCCC\n",
            encoding="utf-8",
        )
        result = parse_fasta_to_dict(fasta)
        assert result == {"seq1": "AAAA", "seq2": "CCCC"}

    def test_empty_file(self, tmp_path: Path) -> None:
        fasta = tmp_path / "empty.fa"
        fasta.write_text("", encoding="utf-8")
        result = parse_fasta_to_dict(fasta)
        assert result == {}

    def test_file_not_found_raises_oserror(self, tmp_path: Path) -> None:
        missing = tmp_path / "nonexistent.fa"
        with pytest.raises(OSError):
            parse_fasta_to_dict(missing)

    def test_duplicate_ids_overwritten(self, tmp_path: Path) -> None:
        fasta = tmp_path / "dup.fa"
        fasta.write_text(
            ">seq1\nAAAA\n>seq1\nCCCC\n",
            encoding="utf-8",
        )
        result = parse_fasta_to_dict(fasta)
        assert result == {"seq1": "CCCC"}
