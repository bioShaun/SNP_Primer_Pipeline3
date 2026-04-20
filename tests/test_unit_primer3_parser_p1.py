"""P1 unit tests for Primer3OutputParser.parse_string."""

from __future__ import annotations

import pytest

from snp_primer_pipeline.core.primer3_parser import Primer3OutputParser

pytestmark = pytest.mark.unit


def _minimal_primer3_record() -> str:
    return """SEQUENCE_ID=test_seq
PRIMER_PAIR_NUM_RETURNED=1
PRIMER_LEFT_0_SEQUENCE=ACGTACGTACGT
PRIMER_LEFT_0=5,12
PRIMER_LEFT_0_TM=60.0
PRIMER_LEFT_0_GC_PERCENT=50.0
PRIMER_LEFT_0_SELF_ANY_TH=0.0
PRIMER_LEFT_0_SELF_END_TH=0.0
PRIMER_LEFT_0_HAIRPIN_TH=0.0
PRIMER_LEFT_0_END_STABILITY=0.0
PRIMER_RIGHT_0_SEQUENCE=TTTTCCCCAAAA
PRIMER_RIGHT_0=100,12
PRIMER_RIGHT_0_TM=60.0
PRIMER_RIGHT_0_GC_PERCENT=50.0
PRIMER_RIGHT_0_SELF_ANY_TH=0.0
PRIMER_RIGHT_0_SELF_END_TH=0.0
PRIMER_RIGHT_0_HAIRPIN_TH=0.0
PRIMER_RIGHT_0_END_STABILITY=0.0
PRIMER_PAIR_0_PRODUCT_SIZE=80
PRIMER_PAIR_0_PENALTY=0.5
PRIMER_PAIR_0_COMPL_ANY_TH=0.0
PRIMER_PAIR_0_COMPL_END_TH=0.0
"""


def test_parse_string_one_pair() -> None:
    parser = Primer3OutputParser()
    out = parser.parse_string(_minimal_primer3_record() + "=\n")
    assert "test_seq" in out
    pairs = out["test_seq"]
    assert len(pairs) == 1
    assert pairs[0].left.sequence == "ACGTACGTACGT"
    assert pairs[0].product_size == 80


def test_parse_string_zero_returned() -> None:
    parser = Primer3OutputParser()
    block = """SEQUENCE_ID=empty
PRIMER_PAIR_NUM_RETURNED=0
=
"""
    out = parser.parse_string(block)
    assert "empty" not in out or out.get("empty") == []
