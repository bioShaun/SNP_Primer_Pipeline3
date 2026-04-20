"""SNP Primer Design Pipeline 3.

A modern, modular pipeline for designing KASP and CAPS/dCAPS primers
for SNP genotyping in any species. This is a refactored version of
SNP_Primer_Pipeline2 with improved architecture and testability.
"""

__version__ = "3.0.0"
__author__ = "Junli Zhang, Claude Assistant"

from .config import PipelineConfig, SoftwarePaths
from .core import (
    AlignedSequence,
    BlastParser,
    BlastRunner,
    FlankingExtractor,
    MultipleSequenceAligner,
    MultipleSequenceAlignment,
    PolymarkerParser,
    Primer3Input,
    Primer3OutputParser,
    Primer3Runner,
    SpecificityAssessor,
    SpecificityBlastRunner,
)
from .main import process_snp, run_pipeline
from .models import SNP, BlastHit, FlankingRegion, Primer, PrimerPair, RestrictionEnzyme
from .primers import CAPSDesigner, KASPDesigner

__all__ = [
    "SNP",
    "AlignedSequence",
    "BlastHit",
    "BlastParser",
    "BlastRunner",
    "CAPSDesigner",
    "FlankingExtractor",
    "FlankingRegion",
    "KASPDesigner",
    "MultipleSequenceAligner",
    "MultipleSequenceAlignment",
    "PipelineConfig",
    "PolymarkerParser",
    "Primer",
    "Primer3Input",
    "Primer3OutputParser",
    "Primer3Runner",
    "PrimerPair",
    "RestrictionEnzyme",
    "SoftwarePaths",
    "SpecificityAssessor",
    "SpecificityBlastRunner",
    "__author__",
    "__version__",
    "process_snp",
    "run_pipeline",
]
