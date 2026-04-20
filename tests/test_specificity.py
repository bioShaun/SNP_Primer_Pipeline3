#!/usr/bin/env python3
"""
Unit tests for KASP primer specificity assessment.

Tests BTOP parsing, Fail/Warning/PASS logic, target exclusion,
and best primer selection using synthetic BLAST output.
"""

import pytest
from pathlib import Path
import tempfile

from snp_primer_pipeline.core.specificity import (
    parse_btop,
    check_three_prime_match,
    OfftargetHit,
    SpecificityAssessor,
    SpecificityBlastRunner,
    SpecificityResult,
    SpecificityStatus,
)

pytestmark = pytest.mark.unit


class TestBtopParsing:
    """Test BTOP string parsing."""

    def test_all_matches(self):
        """BTOP with only matches."""
        ops = parse_btop("20")
        assert ops == [20]

    def test_single_mismatch(self):
        """BTOP with a single mismatch."""
        ops = parse_btop("10AG5")
        assert ops == [10, ("A", "G"), 5]

    def test_multiple_mismatches(self):
        """BTOP with multiple mismatches."""
        ops = parse_btop("3AG2TC1")
        assert ops == [3, ("A", "G"), 2, ("T", "C"), 1]

    def test_gap_in_query(self):
        """BTOP with gap in query (insertion in subject)."""
        ops = parse_btop("5-C3")
        assert ops == [5, ("-", "C"), 3]

    def test_gap_in_subject(self):
        """BTOP with gap in subject."""
        ops = parse_btop("5A-3")
        assert ops == [5, ("A", "-"), 3]

    def test_complex_btop(self):
        """BTOP with mixed operations."""
        ops = parse_btop("6AG3-C2TA")
        assert ops == [6, ("A", "G"), 3, ("-", "C"), 2, ("T", "A")]

    def test_empty_btop(self):
        ops = parse_btop("")
        assert ops == []


class TestThreePrimeMatch:
    """Test 3' end mismatch detection from BTOP."""

    def test_perfect_match_all(self):
        """All bases match — 3' end should be clean."""
        assert check_three_prime_match("20", 1, 20, 20, n_bases=2) is True

    def test_mismatch_at_3prime_end(self):
        """Mismatch at the very last base."""
        # 18 matches then AG mismatch = mismatch at pos 19 (in 20bp primer pos 1-20)
        assert check_three_prime_match("18AG", 1, 20, 20, n_bases=2) is False

    def test_mismatch_at_penultimate(self):
        """Mismatch at penultimate position."""
        # 17 matches, mismatch, 1 match = mismatch at pos 18
        assert check_three_prime_match("17AG1", 1, 20, 20, n_bases=2) is False

    def test_mismatch_far_from_3prime(self):
        """Mismatch far from 3' end — should pass."""
        # Mismatch at pos 3, rest matches
        assert check_three_prime_match("2AG17", 1, 20, 20, n_bases=2) is True

    def test_alignment_not_reaching_3prime(self):
        """Alignment doesn't reach 3' end of primer."""
        # Primer is 25bp but alignment only covers 1-20
        assert check_three_prime_match("20", 1, 20, 25, n_bases=2) is False

    def test_full_match_short_primer(self):
        """Short primer fully matched."""
        assert check_three_prime_match("18", 1, 18, 18, n_bases=2) is True


class TestSpecificityAssessor:
    """Test the SpecificityAssessor logic."""

    def _make_hit(self, **kwargs):
        """Helper to create OfftargetHit with defaults."""
        defaults = {
            "query_id": "SNP1_1-0_F1",
            "subject_id": "chr2",
            "pident": 100.0,
            "length": 20,
            "mismatch": 0,
            "gapopen": 0,
            "qstart": 1,
            "qend": 20,
            "sstart": 1000,
            "send": 1019,
            "qlen": 20,
            "btop": "20",
        }
        defaults.update(kwargs)
        return OfftargetHit(**defaults)

    def _make_kasp_results(self, base_key="1-0"):
        """Create minimal KASP result dicts (group of 3)."""
        return [
            {"index": f"{base_key}-Allele-A", "primer_seq": "ACGTACGTACGTACGTACGT", "score": 5.0, "product_size": 100},
            {"index": f"{base_key}-Allele-B", "primer_seq": "ACGTACGTACGTACGTACGA", "score": 5.0, "product_size": 100},
            {"index": f"{base_key}-Common", "primer_seq": "TGCATGCATGCATGCATGCA", "score": 5.0, "product_size": 100},
        ]

    def test_pass_no_hits(self):
        """No BLAST hits at all → PASS."""
        assessor = SpecificityAssessor()
        results = assessor.assess(
            self._make_kasp_results(), {}, "SNP1", "chr1", 50000
        )
        assert "1-0" in results
        assert results["1-0"].status == SpecificityStatus.PASS

    def test_pass_hits_on_target(self):
        """Hits within target window → excluded → PASS."""
        assessor = SpecificityAssessor(target_window=2000)
        blast_hits = {
            "SNP1_1-0_F1": [self._make_hit(subject_id="chr1", sstart=50010, send=50029)],
            "SNP1_1-0_R": [self._make_hit(query_id="SNP1_1-0_R", subject_id="chr1", sstart=50100, send=50081)],
        }
        results = assessor.assess(
            self._make_kasp_results(), blast_hits, "SNP1", "chr1", 50000
        )
        assert results["1-0"].status == SpecificityStatus.PASS

    def test_fail_off_target_pair(self):
        """F+R pair on off-target contig with all criteria met → FAIL."""
        assessor = SpecificityAssessor()
        # F1 on chr2 forward strand, R on chr2 reverse strand, 200bp apart
        blast_hits = {
            "SNP1_1-0_F1": [
                self._make_hit(
                    subject_id="chr2", sstart=5000, send=5019,
                    pident=100.0, length=20, mismatch=0, btop="20"
                )
            ],
            "SNP1_1-0_R": [
                self._make_hit(
                    query_id="SNP1_1-0_R", subject_id="chr2",
                    sstart=5219, send=5200,  # reverse strand
                    pident=100.0, length=20, mismatch=0, btop="20"
                )
            ],
        }
        results = assessor.assess(
            self._make_kasp_results(), blast_hits, "SNP1", "chr1", 50000
        )
        assert results["1-0"].status == SpecificityStatus.FAIL
        assert len(results["1-0"].fail_details) > 0

    def test_fail_blocked_by_3prime_mismatch(self):
        """Pair with 3' end mismatch → should NOT fail (param A blocks)."""
        assessor = SpecificityAssessor()
        blast_hits = {
            "SNP1_1-0_F1": [
                self._make_hit(
                    subject_id="chr2", sstart=5000, send=5019,
                    pident=90.0, length=20, mismatch=2,
                    btop="18AG"  # mismatch at 3' end
                )
            ],
            "SNP1_1-0_R": [
                self._make_hit(
                    query_id="SNP1_1-0_R", subject_id="chr2",
                    sstart=5219, send=5200,
                    pident=100.0, length=20, mismatch=0, btop="20"
                )
            ],
        }
        results = assessor.assess(
            self._make_kasp_results(), blast_hits, "SNP1", "chr1", 50000
        )
        assert results["1-0"].status == SpecificityStatus.PASS

    def test_fail_blocked_by_distance(self):
        """Pair too far apart → should NOT fail (param C blocks)."""
        assessor = SpecificityAssessor()
        blast_hits = {
            "SNP1_1-0_F1": [
                self._make_hit(subject_id="chr2", sstart=5000, send=5019, btop="20")
            ],
            "SNP1_1-0_R": [
                self._make_hit(
                    query_id="SNP1_1-0_R", subject_id="chr2",
                    sstart=8000, send=7981,  # >1000bp away
                    btop="20"
                )
            ],
        }
        results = assessor.assess(
            self._make_kasp_results(), blast_hits, "SNP1", "chr1", 50000
        )
        assert results["1-0"].status == SpecificityStatus.PASS

    def test_fail_blocked_by_same_strand(self):
        """Both on same strand → should NOT fail (param D blocks)."""
        assessor = SpecificityAssessor()
        blast_hits = {
            "SNP1_1-0_F1": [
                self._make_hit(subject_id="chr2", sstart=5000, send=5019, btop="20")
            ],
            "SNP1_1-0_R": [
                self._make_hit(
                    query_id="SNP1_1-0_R", subject_id="chr2",
                    sstart=5200, send=5219,  # same strand!
                    btop="20"
                )
            ],
        }
        results = assessor.assess(
            self._make_kasp_results(), blast_hits, "SNP1", "chr1", 50000
        )
        assert results["1-0"].status == SpecificityStatus.PASS

    def test_warning_high_hit_count(self):
        """Single primer with many off-target hits → WARNING."""
        assessor = SpecificityAssessor(warning_threshold=5)
        # Create 10 hits on different contigs (all identity >= 85%, coverage >= 80%)
        hits = []
        for j in range(10):
            hits.append(
                self._make_hit(
                    query_id="SNP1_1-0_F1",
                    subject_id=f"scaffold_{j}",
                    pident=90.0,
                    length=18,
                    mismatch=2,
                    qlen=20,
                    btop="18AG",
                )
            )
        blast_hits = {"SNP1_1-0_F1": hits}
        results = assessor.assess(
            self._make_kasp_results(), blast_hits, "SNP1", "chr1", 50000
        )
        assert results["1-0"].status == SpecificityStatus.WARNING

    def test_warning_escalation_to_fail(self):
        """Extreme repeat count → escalated to FAIL."""
        assessor = SpecificityAssessor(warning_threshold=5)
        # 26 hits → >5*5=25 → escalate to FAIL
        hits = []
        for j in range(26):
            hits.append(
                self._make_hit(
                    query_id="SNP1_1-0_F1",
                    subject_id=f"scaffold_{j}",
                    pident=90.0,
                    length=18,
                    mismatch=2,
                    qlen=20,
                    btop="18AG",
                )
            )
        blast_hits = {"SNP1_1-0_F1": hits}
        results = assessor.assess(
            self._make_kasp_results(), blast_hits, "SNP1", "chr1", 50000
        )
        assert results["1-0"].status == SpecificityStatus.FAIL

    def test_target_exclusion(self):
        """Hits within target window should be excluded from warning count."""
        assessor = SpecificityAssessor(warning_threshold=5, target_window=2000)
        # 10 hits, but all on target chromosome within window
        hits = []
        for j in range(10):
            hits.append(
                self._make_hit(
                    query_id="SNP1_1-0_F1",
                    subject_id="chr1",
                    sstart=50000 + j * 10,
                    send=50000 + j * 10 + 19,
                    pident=95.0,
                    length=20,
                    qlen=20,
                    btop="20",
                )
            )
        blast_hits = {"SNP1_1-0_F1": hits}
        results = assessor.assess(
            self._make_kasp_results(), blast_hits, "SNP1", "chr1", 50000
        )
        assert results["1-0"].status == SpecificityStatus.PASS


class TestBestPrimerSelection:
    """Test best primer selection logic."""

    def _make_kasp_results(self, base_key, score):
        return [
            {"index": f"{base_key}-Allele-A", "primer_seq": "A" * 20, "score": score, "product_size": 100},
            {"index": f"{base_key}-Allele-B", "primer_seq": "A" * 20, "score": score, "product_size": 100},
            {"index": f"{base_key}-Common", "primer_seq": "A" * 20, "score": score, "product_size": 100},
        ]

    def test_pass_beats_warning(self):
        """PASS with lower score should beat WARNING with higher score."""
        kasp = self._make_kasp_results("1-0", 3.0) + self._make_kasp_results("2-0", 8.0)
        spec = {
            "1-0": SpecificityResult(status=SpecificityStatus.PASS),
            "2-0": SpecificityResult(status=SpecificityStatus.WARNING),
        }
        best = SpecificityAssessor.select_best(kasp, spec)
        assert best == "1-0"

    def test_higher_score_wins_same_status(self):
        """Among PASS primers, higher score wins."""
        kasp = self._make_kasp_results("1-0", 3.0) + self._make_kasp_results("2-0", 8.0)
        spec = {
            "1-0": SpecificityResult(status=SpecificityStatus.PASS),
            "2-0": SpecificityResult(status=SpecificityStatus.PASS),
        }
        best = SpecificityAssessor.select_best(kasp, spec)
        assert best == "2-0"

    def test_warning_beats_fail(self):
        """WARNING should beat FAIL."""
        kasp = self._make_kasp_results("1-0", 5.0) + self._make_kasp_results("2-0", 5.0)
        spec = {
            "1-0": SpecificityResult(status=SpecificityStatus.WARNING),
            "2-0": SpecificityResult(status=SpecificityStatus.FAIL),
        }
        best = SpecificityAssessor.select_best(kasp, spec)
        assert best == "1-0"

    def test_no_results(self):
        """Empty results → None."""
        best = SpecificityAssessor.select_best([], {})
        assert best is None


class TestSpecificityBlastRunner:
    """Test primer FASTA preparation."""

    def test_prepare_primer_fasta(self):
        """Test tail stripping and FASTA writing."""
        kasp = [
            {"index": "1-0-Allele-A", "primer_seq": "GAAGGTGACCAAGTTCATGCTACGTACGTACGT", "score": 5.0, "product_size": 100},
            {"index": "1-0-Allele-B", "primer_seq": "GAAGGTCGGAGTCAACGGATTACGTACGTACGA", "score": 5.0, "product_size": 100},
            {"index": "1-0-Common", "primer_seq": "TGCATGCATGCATGCA", "score": 5.0, "product_size": 100},
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            runner = SpecificityBlastRunner(Path("/dummy/ref"), threads=1)
            fasta = runner.prepare_primer_fasta(kasp, "SNP1", Path(tmpdir))
            assert fasta is not None
            content = fasta.read_text()
            # FAM tail should be stripped from Allele-A
            assert "GAAGGTGACCAAGTTCATGCT" not in content
            # HEX tail should be stripped from Allele-B
            assert "GAAGGTCGGAGTCAACGGATT" not in content
            # Stripped sequences should be present
            assert "ACGTACGTACGT" in content
            # Common primer should be present as-is
            assert "TGCATGCATGCATGCA" in content
            # Should have 3 entries
            assert content.count(">") == 3

    def test_no_primers(self):
        """Empty primer list → None."""
        with tempfile.TemporaryDirectory() as tmpdir:
            runner = SpecificityBlastRunner(Path("/dummy/ref"), threads=1)
            fasta = runner.prepare_primer_fasta([], "SNP1", Path(tmpdir))
            assert fasta is None


class TestBlastOutputParsing:
    """Test parsing of synthetic BLAST output."""

    def test_parse_blast_tsv(self):
        """Parse a synthetic blastn-short output file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            # Write a synthetic BLAST hit line
            f.write("SNP1_1-0_F1\tchr2\t100.000\t20\t0\t0\t1\t20\t5000\t5019\t20\t20\n")
            f.write("SNP1_1-0_R\tchr2\t95.000\t20\t1\t0\t1\t20\t5200\t5181\t20\t18AG\n")
            tsv_path = Path(f.name)

        assessor = SpecificityAssessor()
        hits = assessor.parse_blast_output(tsv_path)

        assert "SNP1_1-0_F1" in hits
        assert len(hits["SNP1_1-0_F1"]) == 1
        assert hits["SNP1_1-0_F1"][0].pident == 100.0
        assert hits["SNP1_1-0_F1"][0].btop == "20"

        assert "SNP1_1-0_R" in hits
        assert hits["SNP1_1-0_R"][0].strand == "-"

        tsv_path.unlink()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
