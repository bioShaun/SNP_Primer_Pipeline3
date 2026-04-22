#!/usr/bin/env python3
"""
Main pipeline module for SNP Primer Pipeline.

This module provides the main entry point and orchestrates the entire pipeline.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import TYPE_CHECKING

from .config import PipelineConfig
from .core.alignment import MultipleSequenceAligner
from .core.blast import BlastParser, BlastRunner, FlankingExtractor
from .core.fasta import parse_fasta
from .core.parser import PolymarkerParser
from .core.specificity import SpecificityAssessor, SpecificityBlastRunner
from .exceptions import AlignmentError, BlastError, PipelineError, PrimerDesignError
from .primers.caps import CAPSDesigner
from .primers.kasp import KASPDesigner, KASPResult

if TYPE_CHECKING:
    from typing import Any

MIN_SEQUENCES_FOR_ALIGNMENT = 2


def setup_logging(log_level: str = "INFO") -> None:
    """Setup logging configuration."""
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
        ],
    )


def _prepare_output_dir(config: PipelineConfig) -> None:
    """Create output directory and subdirectories.

    Args:
        config: Pipeline configuration.
    """
    config.output_dir.mkdir(parents=True, exist_ok=True)


def _parse_input_file(
    config: PipelineConfig,
) -> tuple[list[Any], str, dict[str, Any], dict[str, str], dict[str, str] | None]:
    """Parse input file and return SNPs with coordinate-format metadata.

    Args:
        config: Pipeline configuration.

    Returns:
        Tuple of (
            snps, input_format, target_regions_direct,
            target_sequences_direct, exclude_chromosomes
        ).
    """
    logger = logging.getLogger(__name__)
    logger.info("Step 1: Parsing input file...")

    input_format = PolymarkerParser.detect_format(config.input_file)
    logger.info(f"Detected input format: {input_format}")

    parser = PolymarkerParser(config.input_file)

    if input_format == "coordinates":
        snps = parser.parse_coordinates(config.reference_file)
    else:
        snps = parser.parse()

    logger.info(f"Parsed {len(snps)} SNPs")

    target_sequences_direct: dict[str, str] = {}
    target_regions_direct: dict[str, Any] = {}
    exclude_chromosomes: dict[str, str] | None = None

    if input_format == "coordinates":
        logger.info("Step 2b: Directly extracting target flanking from coordinates...")
        target_regions_direct, target_sequences_direct = parser.fetch_target_flanking(
            config.reference_file, config.flanking_size
        )
        exclude_chromosomes = {snp.name: snp.chromosome for snp in snps}

    return snps, input_format, target_regions_direct, target_sequences_direct, exclude_chromosomes


def _run_blast_and_extract(
    config: PipelineConfig,
    snps: list[Any],
    target_regions_direct: dict[str, Any],
    target_sequences_direct: dict[str, str],
    exclude_chromosomes: dict[str, str] | None,
) -> tuple[list[Any], dict[str, str]]:
    """Run BLAST, extract flanking regions, and load sequences.

    Args:
        config: Pipeline configuration.
        snps: List of parsed SNPs.
        target_regions_direct: Direct target regions from coordinate input.
        target_sequences_direct: Direct target sequences from coordinate input.
        exclude_chromosomes: Chromosomes to exclude from BLAST hits.

    Returns:
        Tuple of (flanking_regions, extracted_sequences).
    """
    logger = logging.getLogger(__name__)

    # Step 2: Convert to FASTA for BLAST
    logger.info("Step 2: Converting to FASTA format...")
    fasta_file = config.output_dir / "for_blast.fa"
    parser = PolymarkerParser(config.input_file)
    parser.snps = snps
    parser.to_fasta(fasta_file)

    # Step 3: Run BLAST
    logger.info("Step 3: Running BLAST search...")
    blast_runner = BlastRunner(config.reference_file, config.threads)
    blast_output = config.output_dir / "blast_out.txt"
    blast_runner.run(fasta_file, blast_output)

    # Step 4: Parse BLAST results
    logger.info("Step 4: Parsing BLAST results...")
    blast_parser = BlastParser(blast_output)
    blast_hits = blast_parser.parse()

    # Step 5: Extract flanking regions
    logger.info("Step 5: Extracting flanking regions...")
    flanking_extractor = FlankingExtractor(config.reference_file)

    snp_positions = {snp.name: snp.snp_position for snp in snps}

    flanking_regions = flanking_extractor.extract_flanking_regions(
        blast_hits,
        snp_positions,
        config.flanking_size,
        config.max_hits,
        exclude_chromosomes=exclude_chromosomes,
    )

    logger.info(f"Extracted {len(flanking_regions)} flanking regions from BLAST")

    # Step 6: Extract sequences and merge with direct target sequences
    logger.info("Step 6: Extracting flanking sequences...")
    flanking_fasta = config.output_dir / "flanking_sequences.fa"
    flanking_extractor.extract_sequences(flanking_regions, flanking_fasta)

    extracted_sequences: dict[str, str] = {}
    if flanking_fasta.exists():
        for seq_id, sequence in parse_fasta(flanking_fasta):
            extracted_sequences[seq_id.split()[0]] = sequence

    if target_sequences_direct:
        extracted_sequences.update(target_sequences_direct)
        for _snp_name, region in target_regions_direct.items():
            flanking_regions.insert(0, region)
        logger.info(f"Merged {len(target_sequences_direct)} direct target sequences")

    return flanking_regions, extracted_sequences


def _merge_and_write_kasp_results(
    config: PipelineConfig,
    all_kasp_results: list[KASPResult],
    all_specificity_results: dict[str, Any],
    all_best_primer_keys: dict[str, str],
) -> None:
    """Merge and write KASP results for all SNPs.

    Args:
        config: Pipeline configuration.
        all_kasp_results: Combined KASP primer results.
        all_specificity_results: Combined specificity assessment results.
        all_best_primer_keys: Mapping of SNP name to best primer key.
    """
    logger = logging.getLogger(__name__)
    logger.info("Step 8: Writing merged KASP results...")

    kasp_designer = KASPDesigner(
        max_tm=config.max_tm,
        max_size=config.max_primer_size,
        pick_anyway=config.pick_anyway,
    )

    merged_complete = config.output_dir / "all_KASP_primers.txt"
    merged_simple = config.output_dir / "all_KASP_primers_summary.txt"

    kasp_designer.format_output(
        all_kasp_results,
        "",
        merged_complete,
        None,
        show_variant_sites=False,
        specificity_results=all_specificity_results if all_specificity_results else None,
    )
    kasp_designer.format_simple_output(
        all_kasp_results,
        "",
        merged_simple,
        specificity_results=all_specificity_results if all_specificity_results else None,
        best_primer_key=all_best_primer_keys if all_best_primer_keys else None,
    )


def run_pipeline(config: PipelineConfig) -> None:
    """
    Run the complete SNP primer design pipeline.

    Args:
        config: Pipeline configuration
    """
    logger = logging.getLogger(__name__)

    _prepare_output_dir(config)

    logger.info("Starting SNP Primer Pipeline")
    logger.info(f"Input file: {config.input_file}")
    logger.info(f"Reference: {config.reference_file}")
    logger.info(f"Output directory: {config.output_dir}")

    try:
        snps, _input_format, target_regions_direct, target_sequences_direct, exclude_chromosomes = (
            _parse_input_file(config)
        )

        if not snps:
            logger.warning("No valid SNPs found in input file")
            return

        flanking_regions, extracted_sequences = _run_blast_and_extract(
            config, snps, target_regions_direct, target_sequences_direct, exclude_chromosomes
        )

        logger.info("Step 7: Processing SNPs for primer design...")
        snp_regions: dict[str, list] = {}
        for region in flanking_regions:
            if region.snp_name not in snp_regions:
                snp_regions[region.snp_name] = []
            snp_regions[region.snp_name].append(region)

        all_kasp_results: list[KASPResult] = []
        all_specificity_results: dict[str, Any] = {}
        all_best_primer_keys: dict[str, str] = {}
        failed_snps: list[str] = []
        skipped_snps: list[str] = []
        for snp in snps:
            if snp.name not in snp_regions:
                logger.warning(f"No flanking regions found for SNP {snp.name}")
                skipped_snps.append(snp.name)
                continue

            logger.info(f"Processing SNP: {snp.name}")

            try:
                results = process_snp(snp, snp_regions[snp.name], extracted_sequences, config)
                if results.get("kasp"):
                    for r in results["kasp"]:
                        r.snp_name = snp.name
                    all_kasp_results.extend(results["kasp"])
                    if results.get("specificity"):
                        all_specificity_results.update(results["specificity"])
                    if results.get("best_primer_key"):
                        all_best_primer_keys[snp.name] = results["best_primer_key"]
            except (PipelineError, PrimerDesignError, BlastError, AlignmentError, OSError):
                logger.exception("Failed to process SNP %s", snp.name)
                failed_snps.append(snp.name)
                continue

        if all_kasp_results:
            _merge_and_write_kasp_results(
                config, all_kasp_results, all_specificity_results, all_best_primer_keys
            )

        # Pipeline summary
        processed = len(snps) - len(failed_snps) - len(skipped_snps)
        logger.info(
            f"Pipeline completed: {processed}/{len(snps)} SNPs processed, "
            f"{len(failed_snps)} failed, {len(skipped_snps)} skipped"
        )
        if failed_snps:
            logger.warning(f"Failed SNPs: {', '.join(failed_snps)}")
        if skipped_snps:
            logger.warning(f"Skipped SNPs (no flanking regions): {', '.join(skipped_snps)}")

    # Catches all unhandled exceptions to provide context logging before re-raising
    except Exception:
        logger.exception("Pipeline failed")
        raise


def _prepare_snp_sequences(
    snp: Any,
    flanking_regions: list[Any],
    extracted_sequences: dict[str, str],
    config: PipelineConfig,
) -> tuple[
    str | None,
    int | None,
    int | None,
    str | None,
    dict[str, str],
    Path,
    Any | None,
    list[int],
    list[int],
    dict[int, list[int]],
    str | None,
]:
    """Prepare sequences for a SNP, select target, write FASTA, and run MSA.

    Args:
        snp: SNP object.
        flanking_regions: List of flanking regions for this SNP.
        extracted_sequences: Dictionary of extracted sequences by ID.
        config: Pipeline configuration.

    Returns:
        Tuple of (
            target_sequence, target_snp_position, genomic_start, genomic_strand,
            sequences, snp_output_dir, alignment, sites_diff_all, sites_diff_any,
            diffarray, target_name
        ).
    """
    logger = logging.getLogger(__name__)

    snp_output_dir = config.output_dir / f"SNP_{snp.name}"
    snp_output_dir.mkdir(exist_ok=True)

    snp_fasta = snp_output_dir / f"flanking_{snp.name}.fa"

    sequences: dict[str, str] = {}
    target_sequence: str | None = None
    target_snp_position: int | None = None
    genomic_start: int | None = None
    genomic_strand: str | None = None
    target_name: str | None = None

    with open(snp_fasta, "w") as f:
        for region in flanking_regions:
            seq_id = f"{region.snp_name}_{region.chromosome}_{region.allele}_{region.snp_position_in_region}"

            if seq_id not in extracted_sequences:
                logger.warning(f"Sequence not found for {seq_id}, skipping")
                continue

            sequence = extracted_sequences[seq_id]
            sequences[seq_id] = sequence
            f.write(f">{seq_id}\n{sequence}\n")

            if target_sequence is None and region.chromosome == snp.chromosome:
                target_sequence = sequence
                target_snp_position = region.snp_position_in_region - 1
                genomic_start = region.start
                genomic_strand = region.strand
                target_name = seq_id

    # Fallback: if no match found, use first available region with warning
    if target_sequence is None and flanking_regions:
        for region in flanking_regions:
            seq_id = f"{region.snp_name}_{region.chromosome}_{region.allele}_{region.snp_position_in_region}"
            if seq_id in extracted_sequences:
                target_sequence = extracted_sequences[seq_id]
                target_snp_position = region.snp_position_in_region - 1
                genomic_start = region.start
                genomic_strand = region.strand
                target_name = seq_id
                logger.warning(
                    f"No flanking region on target chromosome {snp.chromosome} for SNP {snp.name}. Using first hit."
                )
                break

    alignment: Any | None = None
    sites_diff_all: list[int] = []
    sites_diff_any: list[int] = []
    diffarray: dict[int, list[int]] = {}

    if target_sequence is None or target_snp_position is None:
        return (
            None,
            None,
            None,
            None,
            sequences,
            snp_output_dir,
            alignment,
            sites_diff_all,
            sites_diff_any,
            diffarray,
            target_name,
        )

    if len(sequences) < MIN_SEQUENCES_FOR_ALIGNMENT:
        logger.info(f"Single sequence for SNP {snp.name}, skipping alignment (no homeologs)")
    else:
        try:
            aligner = MultipleSequenceAligner()
            alignment = aligner.align_file(snp_fasta, snp_output_dir / f"alignment_{snp.name}.fa")

            if target_name is None:
                target_name = next(iter(sequences))
            alignment.set_target_sequence(target_name)

            sites_diff_all, sites_diff_any, diffarray = alignment.find_variant_sites_v2(
                target_name, target_snp_position
            )

            logger.info(f"Found {len(sites_diff_all)} sites that differ from all homeologs")

        except AlignmentError as e:
            logger.warning(f"Alignment failed for SNP {snp.name}: {e}")

    return (
        target_sequence,
        target_snp_position,
        genomic_start,
        genomic_strand,
        sequences,
        snp_output_dir,
        alignment,
        sites_diff_all,
        sites_diff_any,
        diffarray,
        target_name,
    )


def _design_kasp_primers(
    snp: Any,
    target_sequence: str,
    target_snp_position: int,
    genomic_start: int | None,
    genomic_strand: str | None,
    config: PipelineConfig,
    snp_output_dir: Path,
    alignment: Any | None,
    sites_diff_all: list[int],
    diffarray: dict[int, list[int]],
    target_name: str | None,
    flanking_regions: list[Any],
    has_multiple_sequences: bool,
) -> tuple[list[KASPResult] | None, dict[str, Any] | None, str | None]:
    """Design KASP primers and optionally assess specificity.

    Args:
        snp: SNP object.
        target_sequence: Target template sequence.
        target_snp_position: 0-based SNP position in target sequence.
        genomic_start: Genomic start coordinate of target region.
        genomic_strand: Strand of target region.
        config: Pipeline configuration.
        snp_output_dir: SNP-specific output directory.
        alignment: Multiple sequence alignment result.
        sites_diff_all: Sites that differ from all homeologs.
        diffarray: Difference array for variant filtering.
        target_name: Name of target sequence in alignment.
        flanking_regions: List of flanking regions for this SNP.
        has_multiple_sequences: Whether there are multiple sequences for this SNP.

    Returns:
        Tuple of (kasp_primers, specificity_results, best_primer_key).
    """
    logger = logging.getLogger(__name__)
    specificity_results: dict[str, Any] | None = None
    best_primer_key: str | None = None

    try:
        logger.info(f"Designing KASP primers for {snp.name}")

        kasp_designer = KASPDesigner(
            max_tm=config.max_tm,
            max_size=config.max_primer_size,
            pick_anyway=config.pick_anyway,
        )

        kasp_primers = kasp_designer.design_primers(
            target_sequence,
            target_snp_position,
            (snp.allele_a, snp.allele_b),
            config.primer_product_size_range,
            sites_diff_all,
            alignment if has_multiple_sequences else None,
            target_name if has_multiple_sequences else None,
            snp_output_dir,
            diffarray,
            genomic_start=genomic_start,
            genomic_strand=genomic_strand,
        )

        # Step 7.5: Specificity assessment
        if config.run_specificity and kasp_primers:
            try:
                logger.info(f"Assessing primer specificity for {snp.name}")
                spec_runner = SpecificityBlastRunner(config.reference_file, config.threads)
                primer_fasta = spec_runner.prepare_primer_fasta(
                    kasp_primers, snp.name, snp_output_dir
                )
                if primer_fasta:
                    spec_blast_out = snp_output_dir / f"specificity_blast_{snp.name}.tsv"
                    spec_runner.run_blast(primer_fasta, spec_blast_out)

                    assessor = SpecificityAssessor(
                        target_window=config.specificity_target_window,
                        warning_threshold=config.specificity_warning_threshold,
                    )
                    blast_hits = assessor.parse_blast_output(spec_blast_out)

                    target_chrom = flanking_regions[0].chromosome if flanking_regions else ""
                    target_snp_genomic = (
                        genomic_start + target_snp_position if genomic_start else None
                    )

                    specificity_results = assessor.assess(
                        kasp_primers, blast_hits, snp.name, target_chrom, target_snp_genomic
                    )
                    best_primer_key = assessor.select_best(kasp_primers, specificity_results)

                    for key, res in specificity_results.items():
                        logger.info(
                            f"  Primer set {key}: {res.status.value}"
                            + (f" - {res.reason}" if res.reason else "")
                        )
                    if best_primer_key:
                        logger.info(f"  Best primer set: {best_primer_key}")

            except BlastError as e:
                logger.warning(f"Specificity assessment failed for {snp.name}: {e}")

        logger.info(f"Designed {len(kasp_primers) // 3} KASP primer pairs for {snp.name}")

    except PrimerDesignError as e:
        logger.warning(f"KASP design failed for SNP {snp.name}: {e}")
        return None, None, None

    return kasp_primers, specificity_results, best_primer_key


def _write_snp_kasp_results(
    kasp_results: list[KASPResult],
    snp_name: str,
    snp_output_dir: Path,
    specificity_results: dict[str, Any] | None,
    best_primer_key: str | None,
    sites_diff_all: list[int],
    show_variant_sites: bool,
) -> None:
    """Write KASP primer results for a single SNP.

    Args:
        kasp_results: List of KASP primer designs.
        snp_name: Name of the SNP.
        snp_output_dir: SNP-specific output directory.
        specificity_results: Specificity assessment results.
        best_primer_key: Key of the best primer set.
        sites_diff_all: Sites that differ from all homeologs.
        show_variant_sites: Whether to show variant sites in output.
    """
    kasp_designer = KASPDesigner(
        max_tm=63.0,
        max_size=25,
        pick_anyway=False,
    )

    kasp_output = snp_output_dir / f"KASP_primers_{snp_name}.txt"
    kasp_designer.format_output(
        kasp_results,
        snp_name,
        kasp_output,
        sites_diff_all,
        show_variant_sites=show_variant_sites,
        specificity_results=specificity_results,
    )

    kasp_summary = snp_output_dir / f"KASP_primers_{snp_name}_summary.txt"
    kasp_designer.format_simple_output(
        kasp_results,
        snp_name,
        kasp_summary,
        specificity_results=specificity_results,
        best_primer_key=best_primer_key,
    )


def _design_caps_primers(
    snp: Any,
    target_sequence: str,
    target_snp_position: int,
    config: PipelineConfig,
    snp_output_dir: Path,
    sites_diff_all: list[int],
) -> list[dict[str, Any]] | None:
    """Design CAPS primers for a SNP.

    Args:
        snp: SNP object.
        target_sequence: Target template sequence.
        target_snp_position: 0-based SNP position in target sequence.
        config: Pipeline configuration.
        snp_output_dir: SNP-specific output directory.
        sites_diff_all: Sites that differ from all homeologs.

    Returns:
        List of CAPS primer results, or None if design failed.
    """
    logger = logging.getLogger(__name__)

    try:
        logger.info(f"Designing CAPS primers for {snp.name}")

        enzyme_file = Path(__file__).parent.parent.parent / "resources" / "NEB_parsed_REs.txt"

        caps_designer = CAPSDesigner(
            enzyme_file=enzyme_file,
            max_tm=config.max_tm,
            max_size=config.max_primer_size,
            pick_anyway=config.pick_anyway,
        )

        caps_enzymes, dcaps_enzymes = caps_designer.find_usable_enzymes(
            target_sequence, target_snp_position, (snp.allele_a, snp.allele_b), config.max_price
        )

        logger.info(f"Found {len(caps_enzymes)} CAPS and {len(dcaps_enzymes)} dCAPS enzymes")

        all_caps_primers: list[Any] = []

        for enzyme in caps_enzymes + dcaps_enzymes:
            try:
                primers = caps_designer.design_caps_primers(
                    target_sequence,
                    target_snp_position,
                    (snp.allele_a, snp.allele_b),
                    enzyme,
                    (300, 900),
                    sites_diff_all,
                    snp_output_dir,
                )
                all_caps_primers.extend(primers)
            except PrimerDesignError:
                continue

        caps_output = snp_output_dir / f"CAPS_primers_{snp.name}.txt"
        caps_designer.format_output(
            all_caps_primers, snp.name, caps_output, caps_enzymes, dcaps_enzymes, sites_diff_all
        )

        logger.info(f"Designed {len(all_caps_primers)} CAPS primer pairs for {snp.name}")

    except PrimerDesignError as e:
        logger.warning(f"CAPS design failed for SNP {snp.name}: {e}")
        return None

    return all_caps_primers


def process_snp(
    snp: Any,
    flanking_regions: list[Any],
    extracted_sequences: dict[str, str],
    config: PipelineConfig,
) -> dict[str, Any]:
    """
    Process a single SNP for primer design.

    Args:
        snp: SNP object
        flanking_regions: List of flanking regions for this SNP
        extracted_sequences: Dictionary of extracted sequences by ID
        config: Pipeline configuration

    Returns:
        Dictionary mapping design type ('kasp', 'caps') to lists of results
    """
    logger = logging.getLogger(__name__)
    results: dict[str, Any] = {"kasp": [], "caps": []}

    (
        target_sequence,
        target_snp_position,
        genomic_start,
        genomic_strand,
        sequences,
        snp_output_dir,
        alignment,
        sites_diff_all,
        _sites_diff_any,
        diffarray,
        target_name,
    ) = _prepare_snp_sequences(snp, flanking_regions, extracted_sequences, config)

    if target_sequence is None or target_snp_position is None:
        logger.warning(f"No valid target sequence found for SNP {snp.name}")
        return results

    # Design KASP primers
    if config.design_kasp:
        kasp_primers, specificity_results, best_primer_key = _design_kasp_primers(
            snp,
            target_sequence,
            target_snp_position,
            genomic_start,
            genomic_strand,
            config,
            snp_output_dir,
            alignment,
            sites_diff_all,
            diffarray,
            target_name,
            flanking_regions,
            len(sequences) > 1,
        )
        if kasp_primers is not None:
            _write_snp_kasp_results(
                kasp_primers,
                snp.name,
                snp_output_dir,
                specificity_results,
                best_primer_key,
                sites_diff_all,
                config.show_variant_sites,
            )
            results["kasp"] = kasp_primers
            results["specificity"] = specificity_results
            results["best_primer_key"] = best_primer_key

    # Design CAPS primers
    if config.design_caps:
        caps_primers = _design_caps_primers(
            snp, target_sequence, target_snp_position, config, snp_output_dir, sites_diff_all
        )
        if caps_primers is not None:
            results["caps"] = caps_primers

    return results


def main() -> None:
    """Main entry point for command line interface."""
    parser = argparse.ArgumentParser(
        description="SNP Primer Pipeline - Design KASP and CAPS primers for SNPs"
    )

    parser.add_argument("input_file", type=Path, help="Input polymarker CSV file")

    parser.add_argument("reference", type=Path, help="Reference genome BLAST database")

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("output"),
        help="Output directory (default: output)",
    )

    parser.add_argument("--no-kasp", action="store_true", help="Skip KASP primer design")

    parser.add_argument("--no-caps", action="store_true", help="Skip CAPS primer design")

    parser.add_argument(
        "--max-price", type=int, default=200, help="Maximum enzyme price for CAPS (default: 200)"
    )

    parser.add_argument(
        "--max-tm", type=float, default=63.0, help="Maximum primer Tm (default: 63.0)"
    )

    parser.add_argument(
        "--max-size", type=int, default=25, help="Maximum primer size (default: 25)"
    )

    parser.add_argument(
        "--pick-anyway", action="store_true", help="Pick primers even if constraints are violated"
    )

    parser.add_argument(
        "--threads", type=int, default=1, help="Number of threads for BLAST (default: 1)"
    )

    parser.add_argument(
        "--show-variant-sites",
        action="store_true",
        help="Show variant sites in KASP output (default: hidden for cleaner output)",
    )

    parser.add_argument(
        "--no-specificity",
        action="store_true",
        help="Skip off-target specificity assessment",
    )

    parser.add_argument(
        "--specificity-warning-threshold",
        type=int,
        default=50,
        help="Warning threshold for primer repeat hit count (default: 50)",
    )

    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level (default: INFO)",
    )

    args = parser.parse_args()

    # Setup logging
    setup_logging(args.log_level)

    # Create configuration
    config = PipelineConfig(
        input_file=args.input_file,
        reference_file=args.reference,
        output_dir=args.output,
        design_kasp=not args.no_kasp,
        design_caps=not args.no_caps,
        max_price=args.max_price,
        max_tm=args.max_tm,
        max_primer_size=args.max_size,
        pick_anyway=args.pick_anyway,
        threads=args.threads,
        show_variant_sites=args.show_variant_sites,
        run_specificity=not args.no_specificity,
        specificity_warning_threshold=args.specificity_warning_threshold,
    )

    try:
        run_pipeline(config)
    except PipelineError:
        logging.exception("Pipeline failed")
        sys.exit(1)
    except KeyboardInterrupt:
        logging.info("Pipeline interrupted by user")
        sys.exit(1)


if __name__ == "__main__":
    main()
