"""Data models for SNP primer pipeline."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Tuple


class Strand(Enum):
    """DNA strand orientation."""
    PLUS = "plus"
    MINUS = "minus"


@dataclass
class SNP:
    """SNP data model."""
    
    name: str
    chromosome: str
    flanking_sequence: str
    snp_position: int  # 0-based position in flanking sequence
    allele_a: str
    allele_b: str
    iupac_code: str
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: Dict) -> "SNP":
        """Create from dictionary."""
        return cls(**data)
    
    @property
    def alleles(self) -> Tuple[str, str]:
        """Get alleles as tuple."""
        return (self.allele_a, self.allele_b)


@dataclass
class BlastHit:
    """BLAST alignment result."""
    
    query_id: str
    subject_id: str
    identity: float
    alignment_length: int
    mismatches: int
    gap_opens: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    evalue: float
    bit_score: float
    query_seq: str
    subject_seq: str
    subject_length: int
    
    @property
    def strand(self) -> Strand:
        """Get strand orientation."""
        return Strand.PLUS if self.subject_start < self.subject_end else Strand.MINUS
    
    @property
    def percent_identity(self) -> float:
        """Calculate percent identity avoiding gaps."""
        return 100 - (self.mismatches + self.gap_opens) / self.alignment_length * 100
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: Dict) -> "BlastHit":
        """Create from dictionary."""
        return cls(**data)


@dataclass
class FlankingRegion:
    """Flanking sequence region for extraction."""
    
    snp_name: str
    chromosome: str
    start: int
    end: int
    strand: Strand
    snp_position_in_region: int
    allele: str
    
    @property
    def length(self) -> int:
        """Get region length."""
        return abs(self.end - self.start) + 1
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        data = asdict(self)
        data['strand'] = self.strand.value  # Convert enum to string
        return data
    
    @classmethod
    def from_dict(cls, data: Dict) -> "FlankingRegion":
        """Create from dictionary."""
        if isinstance(data['strand'], str):
            data['strand'] = Strand(data['strand'])
        return cls(**data)


@dataclass
class Primer:
    """Primer data model."""
    
    name: str = ""
    sequence: str = ""
    start: int = 0
    end: int = 0
    length: int = 0
    tm: float = 0.0
    gc_percent: float = 0.0
    self_any: float = 0.0
    self_end: float = 0.0
    hairpin: float = 0.0
    end_stability: float = 0.0
    direction: str = ""  # LEFT or RIGHT
    diff_three_all: bool = False
    diff_num: int = 0
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: Dict) -> "Primer":
        """Create from dictionary."""
        return cls(**data)
    
    def reverse_complement(self) -> str:
        """Get reverse complement of primer sequence."""
        complement_map = {
            'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
            'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
            'K': 'M', 'M': 'K', 'N': 'N',
            'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
            'k': 'm', 'm': 'k', 'n': 'n',
        }
        return ''.join(complement_map.get(base, base) for base in reversed(self.sequence))


@dataclass
class PrimerPair:
    """Primer pair data model."""
    
    left: Primer = field(default_factory=Primer)
    right: Primer = field(default_factory=Primer)
    product_size: int = 0
    penalty: float = 0.0
    compl_any: float = 0.0
    compl_end: float = 0.0
    score: float = 0.0
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            'left': self.left.to_dict(),
            'right': self.right.to_dict(),
            'product_size': self.product_size,
            'penalty': self.penalty,
            'compl_any': self.compl_any,
            'compl_end': self.compl_end,
            'score': self.score,
        }
    
    @classmethod
    def from_dict(cls, data: Dict) -> "PrimerPair":
        """Create from dictionary."""
        return cls(
            left=Primer.from_dict(data['left']),
            right=Primer.from_dict(data['right']),
            product_size=data['product_size'],
            penalty=data['penalty'],
            compl_any=data['compl_any'],
            compl_end=data['compl_end'],
            score=data['score'],
        )


@dataclass
class RestrictionEnzyme:
    """Restriction enzyme data model."""
    
    name: str
    sequence: str
    price: int
    caps: bool = False
    dcaps: bool = False
    template_seq: str = ""
    change_pos: Optional[int] = None
    primer_end_positions: List[int] = field(default_factory=list)
    all_positions: List[int] = field(default_factory=list)
    
    @property
    def length(self) -> int:
        """Get recognition sequence length."""
        return len(self.sequence)
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: Dict) -> "RestrictionEnzyme":
        """Create from dictionary."""
        return cls(**data)