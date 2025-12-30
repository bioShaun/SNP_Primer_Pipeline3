"""Polymarker input file parser."""

from __future__ import annotations

import re
from pathlib import Path
from typing import List

from loguru import logger

from ..exceptions import ParseError
from ..models import SNP


class PolymarkerParser:
    """Parser for polymarker format input files."""
    
    # IUPAC ambiguity code mapping
    IUPAC_MAP = {
        "[A/G]": "R", "[G/A]": "R",
        "[C/T]": "Y", "[T/C]": "Y", 
        "[G/C]": "S", "[C/G]": "S",
        "[A/T]": "W", "[T/A]": "W",
        "[G/T]": "K", "[T/G]": "K",
        "[A/C]": "M", "[C/A]": "M",
    }
    
    def __init__(self, input_file: Path):
        """Initialize parser with input file path."""
        self.input_file = Path(input_file)
        self.snps: List[SNP] = []
        
        if not self.input_file.exists():
            raise ParseError(f"Input file not found: {self.input_file}")
    
    def parse(self) -> List[SNP]:
        """Parse the input file and return list of SNP objects."""
        self.snps = []
        line_number = 0
        
        logger.info(f"Parsing polymarker input file: {self.input_file}")
        
        try:
            with open(self.input_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line_number += 1
                    line = line.strip()
                    
                    # Skip empty lines
                    if not line:
                        continue
                    
                    # Skip comment lines
                    if line.startswith('#'):
                        continue
                    
                    try:
                        snp = self._parse_line(line, line_number)
                        if snp:
                            self.snps.append(snp)
                    except ParseError as e:
                        logger.warning(f"Skipping invalid line {line_number}: {e}")
                        continue
                    except Exception as e:
                        logger.warning(f"Unexpected error parsing line {line_number}: {e}")
                        continue
        
        except IOError as e:
            raise ParseError(f"Failed to read input file: {e}")
        
        logger.info(f"Successfully parsed {len(self.snps)} SNPs")
        return self.snps
    
    def _parse_line(self, line: str, line_number: int) -> SNP | None:
        """Parse a single line of polymarker input."""
        # Remove any spaces and split by comma
        line = line.replace(" ", "")
        parts = line.split(",")
        
        if len(parts) != 3:
            raise ParseError(
                f"Expected 3 comma-separated fields, got {len(parts)}",
                line_number=line_number,
                line_content=line
            )
        
        snp_name, chromosome, sequence = parts
        
        # Validate SNP name
        if not snp_name:
            raise ParseError("Empty SNP name", line_number=line_number, line_content=line)
        
        # Replace underscores with hyphens in SNP name
        snp_name = snp_name.replace("_", "-")
        
        # Validate chromosome
        if not chromosome:
            raise ParseError("Empty chromosome", line_number=line_number, line_content=line)
        
        # Find and validate SNP position
        snp_match = re.search(r'\[([ATGC])/([ATGC])\]', sequence)
        if not snp_match:
            raise ParseError(
                "No valid SNP found in sequence (expected format: [A/G])",
                line_number=line_number,
                line_content=line
            )
        
        snp_position = snp_match.start()
        allele_a, allele_b = snp_match.groups()
        
        # Convert to IUPAC code
        bracket_code = snp_match.group(0)
        iupac_code = self.convert_iupac(bracket_code)
        
        # Validate flanking sequence length
        flanking_sequence = sequence.replace(bracket_code, iupac_code)
        if len(flanking_sequence) < 20:
            logger.warning(f"Short flanking sequence ({len(flanking_sequence)} bp) for SNP {snp_name}")
        
        return SNP(
            name=snp_name,
            chromosome=chromosome,
            flanking_sequence=flanking_sequence,
            snp_position=snp_position,
            allele_a=allele_a,
            allele_b=allele_b,
            iupac_code=iupac_code
        )
    
    @staticmethod
    def convert_iupac(bracket_code: str) -> str:
        """Convert bracket format IUPAC code to single letter."""
        return PolymarkerParser.IUPAC_MAP.get(bracket_code, "N")
    
    def to_fasta(self, output_file: Path) -> Path:
        """Convert parsed SNPs to FASTA format for BLAST."""
        if not self.snps:
            raise ParseError("No SNPs to convert. Run parse() first.")
        
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Writing FASTA file: {output_file}")
        
        try:
            with open(output_file, 'w') as f:
                for snp in self.snps:
                    # FASTA header format: >snp_name_chromosome_iupac_code
                    header = f">{snp.name}_{snp.chromosome}_{snp.iupac_code}"
                    f.write(f"{header}\n{snp.flanking_sequence}\n")
        
        except IOError as e:
            raise ParseError(f"Failed to write FASTA file: {e}")
        
        logger.info(f"Successfully wrote {len(self.snps)} sequences to FASTA")
        return output_file
    
    def get_snp_by_name(self, name: str) -> SNP | None:
        """Get SNP by name."""
        for snp in self.snps:
            if snp.name == name:
                return snp
        return None
    
    def filter_by_chromosome(self, chromosome: str) -> List[SNP]:
        """Filter SNPs by chromosome."""
        return [snp for snp in self.snps if snp.chromosome == chromosome]
    
    def get_statistics(self) -> dict:
        """Get parsing statistics."""
        if not self.snps:
            return {"total_snps": 0}
        
        chromosomes = set(snp.chromosome for snp in self.snps)
        alleles = set()
        for snp in self.snps:
            alleles.add(snp.allele_a)
            alleles.add(snp.allele_b)
        
        sequence_lengths = [len(snp.flanking_sequence) for snp in self.snps]
        
        return {
            "total_snps": len(self.snps),
            "chromosomes": len(chromosomes),
            "unique_alleles": len(alleles),
            "min_sequence_length": min(sequence_lengths),
            "max_sequence_length": max(sequence_lengths),
            "avg_sequence_length": sum(sequence_lengths) / len(sequence_lengths),
        }