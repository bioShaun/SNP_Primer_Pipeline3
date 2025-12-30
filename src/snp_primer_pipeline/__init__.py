"""SNP Primer Design Pipeline 3 Claude.

A modern, modular pipeline for designing KASP and CAPS/dCAPS primers
for SNP genotyping in any species. This is a refactored version of
SNP_Primer_Pipeline2 with improved architecture and testability.
"""

__version__ = "3.0.0"
__author__ = "Junli Zhang, Claude Assistant"

from .config import PipelineConfig, SoftwarePaths
from .models import SNP, BlastHit, FlankingRegion, Primer, PrimerPair, RestrictionEnzyme
from .core import (
    PolymarkerParser,
    BlastRunner, BlastParser, FlankingExtractor,
    AlignedSequence, MultipleSequenceAlignment, MultipleSequenceAligner,
    Primer3Input, Primer3Runner, Primer3OutputParser
)
from .primers import KASPDesigner, CAPSDesigner
from .main import run_pipeline, process_snp

__all__ = [
    "__version__",
    "__author__",
    "PipelineConfig",
    "SoftwarePaths",
    "SNP",
    "BlastHit", 
    "FlankingRegion",
    "Primer",
    "PrimerPair",
    "RestrictionEnzyme",
    "PolymarkerParser",
    "BlastRunner",
    "BlastParser", 
    "FlankingExtractor",
    "AlignedSequence",
    "MultipleSequenceAlignment",
    "MultipleSequenceAligner",
    "Primer3Input",
    "Primer3Runner",
    "Primer3OutputParser",
    "KASPDesigner",
    "CAPSDesigner",
    "run_pipeline",
    "process_snp"
]