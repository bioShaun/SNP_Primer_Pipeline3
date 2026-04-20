"""Unit tests for FlankingExtractor BLAST helpers (no external BLAST)."""

from __future__ import annotations

from pathlib import Path

import pytest

from snp_primer_pipeline.core.blast import FlankingExtractor
from snp_primer_pipeline.models import BlastHit, Strand

pytestmark = pytest.mark.unit


def _hit(**overrides: object) -> BlastHit:
    defaults: dict = {
        "query_id": "SNP1_chr1A_A",
        "subject_id": "chr1A",
        "identity": 99.0,
        "alignment_length": 50,
        "mismatches": 0,
        "gap_opens": 0,
        "query_start": 1,
        "query_end": 50,
        "subject_start": 1000,
        "subject_end": 1049,
        "evalue": 1e-50,
        "bit_score": 100.0,
        "query_seq": "A" * 50,
        "subject_seq": "A" * 50,
        "subject_length": 10000,
    }
    defaults.update(overrides)
    return BlastHit(**defaults)


@pytest.fixture
def extractor(tmp_path: Path) -> FlankingExtractor:
    return FlankingExtractor(tmp_path / "dummy_ref")


def test_calculate_subject_snp_position_plus_no_gaps(extractor: FlankingExtractor) -> None:
    hit = _hit(
        query_start=1,
        query_end=10,
        subject_start=100,
        subject_end=109,
        query_seq="ATCGATCGAT",
        subject_seq="ATCGATCGAT",
        alignment_length=10,
    )
    pos = extractor._calculate_subject_snp_position(hit, snp_pos=10)
    assert pos == 109


def test_calculate_subject_snp_position_minus_no_gaps(extractor: FlankingExtractor) -> None:
    hit = _hit(
        query_start=1,
        query_end=10,
        subject_start=200,
        subject_end=191,
        query_seq="ATCGATCGAT",
        subject_seq="ATCGATCGAT",
        alignment_length=10,
    )
    assert hit.strand == Strand.MINUS
    pos = extractor._calculate_subject_snp_position(hit, snp_pos=10)
    assert pos == 191


def test_calculate_subject_snp_position_query_gap(extractor: FlankingExtractor) -> None:
    hit = _hit(
        query_start=1,
        query_end=11,
        subject_start=100,
        subject_end=110,
        query_seq="ATCG-ATCGAT",
        subject_seq="ATCGATCGATC",
        alignment_length=11,
    )
    pos = extractor._calculate_subject_snp_position(hit, snp_pos=10)
    assert pos is not None


def test_create_flanking_region_plus_strand(extractor: FlankingExtractor) -> None:
    hit = _hit(subject_id="chr1A", subject_length=5000)
    region = extractor._create_flanking_region(
        "SNP1", hit, subject_snp_pos=1000, flanking_size=10, allele="A"
    )
    assert region.chromosome == "chr1A"
    assert region.start == 990
    assert region.end == 1010
    assert region.snp_position_in_region == 11


def test_parse_query_id_three_parts(extractor: FlankingExtractor) -> None:
    snp, chrom, allele = extractor._parse_query_id("my_snp_chr7A_B_allele", None)
    assert snp == "my_snp_chr7A"
    assert chrom == "B"
    assert allele == "allele"


def test_parse_query_id_with_known_snps_underscore_chrom(extractor: FlankingExtractor) -> None:
    known = ["long_snp_name"]
    snp, chrom, allele = extractor._parse_query_id("long_snp_name_chr3B_alt", known)
    assert snp == "long_snp_name"
    assert chrom == "chr3B"
    assert allele == "alt"


def test_parse_query_id_malformed_fallback(extractor: FlankingExtractor) -> None:
    snp, chrom, allele = extractor._parse_query_id("bad", None)
    assert snp == "bad"
    assert chrom == "unknown"
    assert allele == "unknown"
