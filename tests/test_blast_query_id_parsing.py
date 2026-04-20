from __future__ import annotations

from pathlib import Path

import pytest

from snp_primer_pipeline.core.blast import FlankingExtractor
from snp_primer_pipeline.models import BlastHit

pytestmark = pytest.mark.unit


def test_extract_flanking_regions_with_underscores_in_coordinate_query_id():
    """Coordinate-mode query IDs can contain underscores in both SNP and chromosome names."""
    extractor = FlankingExtractor(Path("dummy_reference"))
    query_id = "AA_Chr01-16099_AA_Chr01_Y"

    hit = BlastHit(
        query_id=query_id,
        subject_id="AA_Chr01",
        identity=99.010,
        alignment_length=101,
        mismatches=1,
        gap_opens=0,
        query_start=1,
        query_end=101,
        subject_start=16049,
        subject_end=16149,
        evalue=4.88e-44,
        bit_score=180.0,
        query_seq="A" * 101,
        subject_seq="A" * 101,
        subject_length=28865182,
    )

    regions = extractor.extract_flanking_regions(
        blast_hits={query_id: [hit]},
        snp_positions={"AA_Chr01-16099": 50},
        flanking_size=500,
        max_hits=6,
    )

    assert len(regions) == 1
    assert regions[0].snp_name == "AA_Chr01-16099"
    assert regions[0].chromosome == "AA_Chr01"
    assert regions[0].allele == "Y"
