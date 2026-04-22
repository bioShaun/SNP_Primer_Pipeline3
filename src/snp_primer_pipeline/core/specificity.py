#!/usr/bin/env python3
"""
KASP primer specificity assessment module.

This module evaluates the off-target binding risk of designed KASP primers
by running blastn-short against the reference genome and applying Fail/Warning
criteria based on 3' end anchoring, identity/coverage, amplicon distance,
and primer orientation.
"""

from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING

from ..exceptions import BlastError

if TYPE_CHECKING:
    from ..primers.kasp import KASPResult

logger = logging.getLogger(__name__)


# KASP platform tail sequences (5' end of allele-specific primers)
FAM_TAIL = "GAAGGTGACCAAGTTCATGCT"
HEX_TAIL = "GAAGGTCGGAGTCAACGGATT"


class SpecificityStatus(Enum):
    """Specificity assessment result."""

    PASS = "PASS"
    WARNING = "WARNING"
    FAIL = "FAIL"


@dataclass
class SpecificityResult:
    """Result of specificity assessment for a primer set."""

    status: SpecificityStatus = SpecificityStatus.PASS
    reason: str = ""
    fail_details: list[str] = field(default_factory=list)
    warning_details: list[str] = field(default_factory=list)


@dataclass
class OfftargetHit:
    """Parsed BLAST hit for specificity analysis."""

    query_id: str
    subject_id: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    qlen: int
    btop: str

    @property
    def coverage(self) -> float:
        """Alignment coverage as fraction of query length."""
        return self.length / self.qlen if self.qlen > 0 else 0.0

    @property
    def strand(self) -> str:
        """'+' if forward, '-' if reverse."""
        return "+" if self.sstart < self.send else "-"

    @property
    def subject_start_min(self) -> int:
        """Leftmost position on subject."""
        return min(self.sstart, self.send)

    @property
    def subject_end_max(self) -> int:
        """Rightmost position on subject."""
        return max(self.sstart, self.send)


def parse_btop(btop: str) -> list:
    """
    Parse BTOP (Blast Traceback Operations) string into a list of operations.

    Each element is either:
    - int: number of matching bases
    - (query_base, subject_base): a mismatch pair
    - ('-', subject_base): insertion in subject (gap in query)
    - (query_base, '-'): deletion in subject (gap in subject)

    Example BTOP: "6AG3-C2TA" -> [6, ('A','G'), 3, ('-','C'), 2, ('T','A')]
    """
    ops = []
    i = 0
    while i < len(btop):
        # Check for number (matching bases)
        if btop[i].isdigit():
            num_str = ""
            while i < len(btop) and btop[i].isdigit():
                num_str += btop[i]
                i += 1
            ops.append(int(num_str))
        # Gap in query (insertion in subject): -X
        elif btop[i] == "-" and i + 1 < len(btop):
            ops.append(("-", btop[i + 1]))
            i += 2
        # Gap in subject (deletion from subject): X-
        elif i + 1 < len(btop) and btop[i + 1] == "-":
            ops.append((btop[i], "-"))
            i += 2
        # Mismatch: two letters
        elif btop[i].isalpha() and i + 1 < len(btop) and btop[i + 1].isalpha():
            ops.append((btop[i], btop[i + 1]))
            i += 2
        else:
            # Skip unexpected characters
            i += 1
    return ops


def check_three_prime_match(
    btop: str, _query_start: int, query_end: int, qlen: int, n_bases: int = 2
) -> bool:
    """
    Check if the 3' end of the primer has 0 mismatches.

    For a forward-strand hit (qstart < qend), 3' is at the end of BTOP.
    For a reverse-strand hit, we look at the alignment direction.
    In blastn output, qstart < qend always for the query, so 3' end
    of the original primer maps to the end of the alignment.

    Args:
        btop: BTOP string from BLAST
        query_start: qstart from BLAST (1-based)
        query_end: qend from BLAST (1-based)
        qlen: query length
        n_bases: number of 3' bases to check

    Returns:
        True if the last n_bases of the primer alignment are all matches
    """
    ops = parse_btop(btop)
    if not ops:
        return False

    # Build per-position alignment from BTOP
    # Each position in the query is either a match or mismatch/gap
    positions = []  # list of 'match', 'mismatch', or 'gap'
    for op in ops:
        if isinstance(op, int):
            positions.extend(["match"] * op)
        elif isinstance(op, tuple):
            if op[0] == "-":
                # Gap in query — does not consume a query position
                continue
            if op[1] == "-":
                # Gap in subject — consumes a query position
                positions.append("gap")
            else:
                # Mismatch
                positions.append("mismatch")

    if len(positions) == 0:
        return False

    # The alignment covers query positions qstart..qend (1-based)
    # 3' end of the original primer is at position qlen
    # We need to check if the last n_bases of the query that are
    # within the alignment are all matches.

    # If alignment doesn't reach the 3' end, it's not a match
    if query_end < qlen:
        return False

    # Check the last n_bases positions in the alignment
    if len(positions) < n_bases:
        return False

    return all(p == "match" for p in positions[-n_bases:])


class SpecificityBlastRunner:
    """Run blastn-short on KASP primers for off-target assessment."""

    BLAST_FIELD_COUNT = 12

    def __init__(self, reference: Path, threads: int = 1):
        self.reference = Path(reference)
        self.threads = threads

    def prepare_primer_fasta(
        self,
        kasp_results: list[KASPResult],
        snp_name: str,
        output_dir: Path,
    ) -> Path | None:
        """
        Strip FAM/HEX tails from KASP primers and write to FASTA.

        Args:
            kasp_results: KASP primer result dicts (groups of 3: alleleA, alleleB, common)
            snp_name: SNP identifier
            output_dir: directory to write FASTA into

        Returns:
            Path to primers_no_tail.fasta, or None if no primers
        """
        if not kasp_results:
            return None

        fasta_path = output_dir / f"primers_no_tail_{snp_name}.fasta"
        written = set()

        with open(fasta_path, "w") as f:
            for i in range(0, len(kasp_results), 3):
                if i + 2 >= len(kasp_results):
                    break

                allele_a = kasp_results[i]
                allele_b = kasp_results[i + 1]
                common = kasp_results[i + 2]

                # Base key for this primer set
                base_key = allele_a.index.rsplit("-", 2)[0]

                # Allele-specific primers: strip tails
                for label, result in [("F1", allele_a), ("F2", allele_b)]:
                    seq = result.primer_seq.upper()
                    # Strip FAM or HEX tail
                    for tail in [FAM_TAIL, HEX_TAIL]:
                        if seq.startswith(tail):
                            seq = seq[len(tail) :]
                            break
                    primer_id = f"{snp_name}_{base_key}_{label}"
                    if primer_id not in written:
                        f.write(f">{primer_id}\n{seq}\n")
                        written.add(primer_id)

                # Common primer: no tail
                seq = common.primer_seq.upper()
                primer_id = f"{snp_name}_{base_key}_R"
                if primer_id not in written:
                    f.write(f">{primer_id}\n{seq}\n")
                    written.add(primer_id)

        if not written:
            return None

        return fasta_path

    def run_blast(self, query_fasta: Path, output_file: Path) -> Path:
        """
        Run blastn-short for off-target detection.

        Returns:
            Path to BLAST output TSV
        """
        cmd = [
            "blastn",
            "-task",
            "blastn-short",
            "-query",
            str(query_fasta),
            "-db",
            str(self.reference),
            "-evalue",
            "1000",
            "-word_size",
            "7",
            "-dust",
            "no",
            "-soft_masking",
            "false",
            "-max_target_seqs",
            "10000",
            "-max_hsps",
            "1",
            "-outfmt",
            "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen btop",
            "-num_threads",
            str(self.threads),
            "-out",
            str(output_file),
        ]

        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            raise BlastError(f"Specificity BLAST failed: {e.stderr}") from e
        except FileNotFoundError as e:
            raise BlastError(f"BLAST not found in PATH: {e}") from e
        else:
            return output_file


class SpecificityAssessor:
    """
    Assess KASP primer specificity from blastn-short results.

    Applies Fail (off-target amplification) and Warning (primer sponge)
    criteria per the assessment standard document.
    """

    BLAST_FIELD_COUNT = 12

    def __init__(
        self,
        target_window: int = 2000,
        warning_threshold: int = 50,
        fail_identity: float = 85.0,
        fail_coverage: float = 0.85,
        fail_max_mismatch: int = 3,
        fail_distance_range: tuple[int, int] = (50, 1000),
        fail_three_prime_bases: int = 2,
        warning_identity: float = 85.0,
        warning_coverage: float = 0.80,
    ):
        self.target_window = target_window
        self.warning_threshold = warning_threshold
        self.fail_identity = fail_identity
        self.fail_coverage = fail_coverage
        self.fail_max_mismatch = fail_max_mismatch
        self.fail_distance_range = fail_distance_range
        self.fail_three_prime_bases = fail_three_prime_bases
        self.warning_identity = warning_identity
        self.warning_coverage = warning_coverage

    def parse_blast_output(self, blast_file: Path) -> dict[str, list[OfftargetHit]]:
        """
        Parse blastn-short output into per-primer hit lists.

        Returns:
            Dict mapping primer_id -> list of OfftargetHit
        """
        hits: dict[str, list[OfftargetHit]] = {}

        if not blast_file.exists():
            return hits

        with open(blast_file) as f:
            for raw_line in f:
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < self.BLAST_FIELD_COUNT:
                    continue
                try:
                    hit = OfftargetHit(
                        query_id=fields[0],
                        subject_id=fields[1],
                        pident=float(fields[2]),
                        length=int(fields[3]),
                        mismatch=int(fields[4]),
                        gapopen=int(fields[5]),
                        qstart=int(fields[6]),
                        qend=int(fields[7]),
                        sstart=int(fields[8]),
                        send=int(fields[9]),
                        qlen=int(fields[10]),
                        btop=fields[11],
                    )
                    hits.setdefault(hit.query_id, []).append(hit)
                except (ValueError, IndexError):
                    continue

        return hits

    def _is_on_target(
        self, hit: OfftargetHit, target_chrom: str, target_snp_pos: int | None
    ) -> bool:
        """Check if a hit falls within the target region."""
        if hit.subject_id != target_chrom:
            return False
        if target_snp_pos is None:
            return True  # Same chromosome, no position info → assume on-target
        hit_mid = (hit.subject_start_min + hit.subject_end_max) // 2
        return abs(hit_mid - target_snp_pos) <= self.target_window

    def _check_fail_single_hit(self, hit: OfftargetHit) -> bool:
        """Check if a single hit passes Fail criteria A + B (individual primer)."""
        # Parameter B: identity & coverage & mismatch cap
        if hit.pident < self.fail_identity:
            return False
        if hit.coverage < self.fail_coverage:
            return False
        if hit.mismatch > self.fail_max_mismatch:
            return False

        # Parameter A: 3' end 2bp perfect match
        return bool(
            check_three_prime_match(
                hit.btop, hit.qstart, hit.qend, hit.qlen, self.fail_three_prime_bases
            )
        )

    def _check_fail_pair(self, f_hit: OfftargetHit, r_hit: OfftargetHit) -> bool:
        """
        Check Fail criteria C + D for a pair of hits on the same contig.

        Args:
            f_hit: hit from an F primer (allele-specific)
            r_hit: hit from the R primer (common)

        Returns:
            True if this pair constitutes a FAIL
        """
        # Must be on the same contig
        if f_hit.subject_id != r_hit.subject_id:
            return False

        # Parameter D: converging orientation
        # One must be on + strand and the other on - strand
        if f_hit.strand == r_hit.strand:
            return False

        # Determine which is forward and which is reverse
        if f_hit.strand == "+":
            fwd_hit, rev_hit = f_hit, r_hit
        else:
            fwd_hit, rev_hit = r_hit, f_hit

        # 3' end of forward hit is at its rightmost position
        fwd_3prime = fwd_hit.subject_end_max
        # 3' end of reverse hit is at its leftmost position
        rev_3prime = rev_hit.subject_start_min

        # For converging amplification, fwd 3' should be left of rev 3'
        if fwd_3prime >= rev_3prime:
            return False

        # Parameter C: distance between 3' ends
        distance = rev_3prime - fwd_3prime
        return bool(self.fail_distance_range[0] <= distance <= self.fail_distance_range[1])

    def _count_effective_hits(
        self,
        hits: list[OfftargetHit],
        target_chrom: str,
        target_snp_pos: int | None,
    ) -> int:
        """Count effective hits for Warning assessment."""
        count = 0
        for hit in hits:
            if self._is_on_target(hit, target_chrom, target_snp_pos):
                continue
            if hit.pident >= self.warning_identity and hit.coverage >= self.warning_coverage:
                count += 1
        return count

    def assess(
        self,
        kasp_results: list[KASPResult],
        blast_hits: dict[str, list[OfftargetHit]],
        snp_name: str,
        target_chrom: str,
        target_snp_pos: int | None = None,
    ) -> dict[str, SpecificityResult]:
        """
        Assess specificity for all primer sets of a SNP.

        Args:
            kasp_results: list of KASP result dicts (groups of 3)
            blast_hits: parsed BLAST hits keyed by primer_id
            snp_name: SNP name
            target_chrom: target chromosome ID
            target_snp_pos: genomic position of SNP on target chromosome

        Returns:
            Dict mapping primer_set_base_key -> SpecificityResult
        """
        results: dict[str, SpecificityResult] = {}

        for i in range(0, len(kasp_results), 3):
            if i + 2 >= len(kasp_results):
                break

            allele_a = kasp_results[i]
            base_key = allele_a.index.rsplit("-", 2)[0]

            f1_id = f"{snp_name}_{base_key}_F1"
            f2_id = f"{snp_name}_{base_key}_F2"
            r_id = f"{snp_name}_{base_key}_R"

            result = SpecificityResult()

            # --- Fail check ---
            # Get off-target hits that pass individual criteria A+B
            f_primer_ids = [f1_id, f2_id]
            r_primer_id = r_id

            r_hits_raw = blast_hits.get(r_primer_id, [])
            r_offtarget = [
                h for h in r_hits_raw if not self._is_on_target(h, target_chrom, target_snp_pos)
            ]
            r_qualified = [h for h in r_offtarget if self._check_fail_single_hit(h)]

            for f_id in f_primer_ids:
                f_hits_raw = blast_hits.get(f_id, [])
                f_offtarget = [
                    h for h in f_hits_raw if not self._is_on_target(h, target_chrom, target_snp_pos)
                ]
                f_qualified = [h for h in f_offtarget if self._check_fail_single_hit(h)]

                for fh in f_qualified:
                    for rh in r_qualified:
                        if self._check_fail_pair(fh, rh):
                            detail = (
                                f"{f_id} + {r_primer_id} off-target on "
                                f"{fh.subject_id}:{fh.subject_start_min}-{rh.subject_end_max}"
                            )
                            result.fail_details.append(detail)

            if result.fail_details:
                result.status = SpecificityStatus.FAIL
                result.reason = (
                    f"Off-target amplification detected ({len(result.fail_details)} site(s))"
                )
                results[base_key] = result
                continue

            # --- Warning check ---
            for pid in [f1_id, f2_id, r_id]:
                hits = blast_hits.get(pid, [])
                count = self._count_effective_hits(hits, target_chrom, target_snp_pos)
                if count >= self.warning_threshold * 5:
                    # Escalate to FAIL
                    result.fail_details.append(
                        f"{pid}: {count} effective hits (>{self.warning_threshold * 5})"
                    )
                    result.status = SpecificityStatus.FAIL
                    result.reason = f"Extreme repeat: {pid} has {count} genome-wide hits"
                elif count >= self.warning_threshold:
                    result.warning_details.append(
                        f"{pid}: {count} effective hits (>={self.warning_threshold})"
                    )

            if result.status == SpecificityStatus.FAIL:
                results[base_key] = result
                continue

            if result.warning_details:
                result.status = SpecificityStatus.WARNING
                result.reason = f"Primer sponge risk ({len(result.warning_details)} primer(s))"

            results[base_key] = result

        return results

    @staticmethod
    def select_best(
        kasp_results: list[KASPResult],
        specificity_results: dict[str, SpecificityResult],
    ) -> str | None:
        """
        Select the best primer set based on specificity + score.

        Priority: PASS > WARNING > FAIL, then highest score wins.

        Returns:
            base_key of the best primer set, or None
        """
        candidates: list[tuple[int, float, str]] = []

        # Status priority: PASS=0, WARNING=1, FAIL=2
        status_priority = {
            SpecificityStatus.PASS: 0,
            SpecificityStatus.WARNING: 1,
            SpecificityStatus.FAIL: 2,
        }

        for i in range(0, len(kasp_results), 3):
            if i + 2 >= len(kasp_results):
                break

            allele_a = kasp_results[i]
            base_key = allele_a.index.rsplit("-", 2)[0]
            score = allele_a.score
            if isinstance(score, str):
                try:
                    score = float(score)
                except ValueError:
                    score = 0.0

            spec = specificity_results.get(base_key)
            priority = status_priority.get(spec.status if spec else SpecificityStatus.PASS, 0)

            # Lower priority number is better, higher score is better
            candidates.append((priority, -score, base_key))

        if not candidates:
            return None

        candidates.sort()
        return candidates[0][2]
