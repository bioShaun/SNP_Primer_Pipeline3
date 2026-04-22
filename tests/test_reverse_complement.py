"""Unit tests for reverse_complement — verifies the single canonical implementation."""

from __future__ import annotations

import pytest

from snp_primer_pipeline.models import reverse_complement


class TestReverseComplement:
    """Tests for the canonical reverse_complement function."""

    # --- basic ATGC ---

    def test_basic_atgc(self):
        assert reverse_complement("ATGC") == "GCAT"

    def test_lowercase_atgc(self):
        assert reverse_complement("atgc") == "gcat"

    def test_mixed_case(self):
        assert reverse_complement("AaTt") == "aAtT"

    def test_empty_string(self):
        assert reverse_complement("") == ""

    def test_single_base(self):
        assert reverse_complement("A") == "T"
        assert reverse_complement("C") == "G"

    # --- IUPAC ambiguity codes ---

    def test_iupac_ry(self):
        """R (A/G) ↔ Y (T/C)."""
        assert reverse_complement("R") == "Y"
        assert reverse_complement("Y") == "R"

    def test_iupac_sw(self):
        """S (G/C) is self-complementary; W (A/T) is self-complementary."""
        assert reverse_complement("S") == "S"
        assert reverse_complement("W") == "W"

    def test_iupac_km(self):
        """K (G/T) ↔ M (A/C)."""
        assert reverse_complement("K") == "M"
        assert reverse_complement("M") == "K"

    def test_iupac_bvdh(self):
        """B (not A) ↔ V (not T); D (not C) ↔ H (not G)."""
        assert reverse_complement("B") == "V"
        assert reverse_complement("V") == "B"
        assert reverse_complement("D") == "H"
        assert reverse_complement("H") == "D"

    def test_iupac_n(self):
        """N is self-complementary."""
        assert reverse_complement("N") == "N"

    def test_iupac_lowercase(self):
        assert reverse_complement("r") == "y"
        assert reverse_complement("y") == "r"
        assert reverse_complement("b") == "v"
        assert reverse_complement("d") == "h"

    # --- realistic sequences ---

    def test_primer_like_sequence(self):
        seq = "GAAGGTGACCAAGTTCATGCT"
        rc = reverse_complement(seq)
        assert len(rc) == len(seq)
        # double RC should return the original
        assert reverse_complement(rc) == seq

    def test_unknown_char_passthrough(self):
        """Characters not in the complement table are kept as-is."""
        assert reverse_complement("AXA") == "TXT"
