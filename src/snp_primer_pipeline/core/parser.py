"""Polymarker input file parser."""

from __future__ import annotations

import re
import subprocess
import tempfile
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
    
    @staticmethod
    def detect_format(input_file: Path) -> str:
        """
        Detect input file format.
        
        Returns:
            'polymarker' or 'coordinates'
        """
        input_file = Path(input_file)
        if not input_file.exists():
            raise ParseError(f"Input file not found: {input_file}")
        
        try:
            with open(input_file, 'r', encoding='utf-8') as f:
                # Read first non-empty, non-comment line
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    # Dynamic split: if comma is present, split by comma; otherwise split by whitespace
                    parts = [p.strip() for p in line.split(',')] if ',' in line else line.split()
                    
                    # Check for polymarker format (must contain brackets and have 3 fields)
                    if '[' in line and ']' in line and len(parts) == 3:
                        return 'polymarker'
                    
                    # Check for coordinate format (4 fields)
                    if len(parts) == 4:
                        # Validate that position is numeric
                        try:
                            int(parts[1])
                            return 'coordinates'
                        except ValueError:
                            pass
                    
                    # If first line doesn't match either format, it's an error
                    raise ParseError(
                        f"Unknown input format. Expected either:\n"
                        f"  - Polymarker (CSV/Space/Tab): SNP_ID Chromosome Sequence[A/G]Sequence\n"
                        f"  - Coordinates (CSV/Space/Tab): Chromosome Position Ref Alt"
                    )
        except IOError as e:
            raise ParseError(f"Failed to read input file: {e}")
        
        raise ParseError("Empty input file")
    
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
        # Split by comma if present, else by whitespace
        if ',' in line:
            parts = [p.strip() for p in line.split(',')]
        else:
            parts = line.split()
            
        # Remove internal spaces in strings if any (polymarker usually doesn't have internal spaces)
        parts = [p.replace(" ", "") for p in parts]
        
        if len(parts) != 3:
            raise ParseError(
                f"Expected 3 fields, got {len(parts)}",
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
    
    def parse_coordinates(self, reference: Path) -> List[SNP]:
        """
        Parse coordinate-based input file and fetch sequences from reference.
        
        Args:
            reference: Path to BLAST reference database
            
        Returns:
            List of SNP objects with fetched sequences
        """
        logger.info(f"Parsing coordinate input file: {self.input_file}")
        
        # Step 1: Parse coordinates
        self._coordinates = []
        line_number = 0
        
        try:
            with open(self.input_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line_number += 1
                    line = line.strip()
                    
                    if not line or line.startswith('#'):
                        continue
                    
                    if ',' in line:
                        parts = [p.strip() for p in line.split(',')]
                    else:
                        parts = line.split()
                        
                    if len(parts) != 4:
                        logger.warning(f"Skipping line {line_number}: expected 4 fields, got {len(parts)}")
                        continue
                    
                    chrom, pos_str, ref, alt = parts
                    
                    try:
                        pos = int(pos_str)
                    except ValueError:
                        logger.warning(f"Skipping line {line_number}: invalid position '{pos_str}'")
                        continue
                    
                    # Validate alleles
                    if ref.upper() not in 'ATGC' or alt.upper() not in 'ATGC':
                        logger.warning(f"Skipping line {line_number}: invalid alleles {ref}/{alt}")
                        continue
                    
                    self._coordinates.append({
                        'chrom': chrom,
                        'pos': pos,
                        'ref': ref.upper(),
                        'alt': alt.upper(),
                        'name': f"{chrom}-{pos}"
                    })
        except IOError as e:
            raise ParseError(f"Failed to read input file: {e}")
        
        logger.info(f"Parsed {len(self._coordinates)} coordinate records")
        
        if not self._coordinates:
            return []
        
        # Step 2: Fetch short sequences for BLAST (±50bp, used to find homeologs)
        logger.info("Fetching short flanking sequences for homeolog search...")
        sequences = self._fetch_sequences(self._coordinates, reference, flank=50)
        
        # Step 3: Create SNP objects
        self.snps = []
        for coord in self._coordinates:
            seq = sequences.get(coord['name'])
            if not seq:
                logger.warning(f"No sequence retrieved for {coord['name']}, skipping")
                continue
            
            # Insert [Ref/Alt] at position 50 (center of 101bp sequence)
            snp_seq = seq[:50] + f"[{coord['ref']}/{coord['alt']}]" + seq[51:]
            
            # Convert to IUPAC
            bracket_code = f"[{coord['ref']}/{coord['alt']}]"
            iupac_code = self.convert_iupac(bracket_code)
            flanking_seq = seq[:50] + iupac_code + seq[51:]
            
            snp = SNP(
                name=coord['name'],
                chromosome=coord['chrom'],
                flanking_sequence=flanking_seq,
                snp_position=50,  # SNP is at center of fetched sequence
                allele_a=coord['ref'],
                allele_b=coord['alt'],
                iupac_code=iupac_code
            )
            self.snps.append(snp)
        
        logger.info(f"Successfully created {len(self.snps)} SNP objects")
        return self.snps

    def fetch_target_flanking(self, reference: Path, flanking_size: int = 500):
        """
        Directly extract target flanking regions from known coordinates.

        This skips the BLAST→re-extract cycle for the target chromosome.
        Must be called after parse_coordinates().

        Args:
            reference: Path to BLAST database
            flanking_size: bp to extract on each side of the SNP

        Returns:
            Tuple of (target_regions, target_sequences):
                target_regions: dict mapping snp_name -> FlankingRegion
                target_sequences: dict mapping seq_id -> sequence string
        """
        from ..models import FlankingRegion, Strand

        if not hasattr(self, '_coordinates') or not self._coordinates:
            raise ParseError("No coordinates available. Run parse_coordinates() first.")

        logger.info(f"Directly extracting ±{flanking_size}bp target flanking from known coordinates...")

        # Fetch long flanking sequences
        sequences = self._fetch_sequences(self._coordinates, reference, flank=flanking_size)

        target_regions = {}
        target_sequences = {}

        for coord in self._coordinates:
            seq = sequences.get(coord['name'])
            if not seq:
                logger.warning(f"No target flanking for {coord['name']}, skipping")
                continue

            snp_pos_in_region = min(flanking_size, coord['pos'] - 1) + 1  # 1-based
            genomic_start = max(1, coord['pos'] - flanking_size)

            region = FlankingRegion(
                snp_name=coord['name'],
                chromosome=coord['chrom'],
                start=genomic_start,
                end=coord['pos'] + flanking_size,
                strand=Strand.PLUS,
                snp_position_in_region=snp_pos_in_region,
                allele=coord['ref'],
            )

            # Sequence ID matching the format used by FlankingExtractor
            seq_id = f"{region.snp_name}_{region.chromosome}_{region.allele}_{region.snp_position_in_region}"

            target_regions[coord['name']] = region
            target_sequences[seq_id] = seq

        logger.info(f"Extracted {len(target_regions)} target flanking regions")
        return target_regions, target_sequences

    def _fetch_sequences(self, coordinates: list, reference: Path, flank: int = 50) -> dict:
        """
        Fetch sequences from BLAST database using blastdbcmd.
        
        Args:
            coordinates: List of coordinate dictionaries
            reference: Path to BLAST database
            flank: Number of bp to extract on each side of position
            
        Returns:
            Dictionary mapping SNP names to sequences
        """
        sequences = {}
        
        for coord in coordinates:
            try:
                start = max(1, coord['pos'] - flank)
                end = coord['pos'] + flank
                
                cmd = [
                    "blastdbcmd",
                    "-db", str(reference),
                    "-entry", coord['chrom'],
                    "-range", f"{start}-{end}"
                ]
                
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=True
                )
                
                # Parse FASTA output (skip header line)
                lines = result.stdout.strip().split('\n')
                if len(lines) >= 2:
                    seq = ''.join(lines[1:])  # Join all sequence lines
                    sequences[coord['name']] = seq
                else:
                    logger.warning(f"No sequence retrieved for {coord['name']}")
                    
            except subprocess.CalledProcessError as e:
                logger.warning(f"Failed to fetch sequence for {coord['name']}: {e.stderr}")
                continue
            except Exception as e:
                logger.warning(f"Error fetching sequence for {coord['name']}: {e}")
                continue
        
        return sequences
    
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