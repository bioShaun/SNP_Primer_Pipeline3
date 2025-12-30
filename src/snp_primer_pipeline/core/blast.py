#!/usr/bin/env python3
"""
BLAST processing module for SNP Primer Pipeline.

This module handles BLAST execution, result parsing, and flanking region extraction.
"""

import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import Counter
from dataclasses import dataclass

from ..models import BlastHit, FlankingRegion, Strand
from ..exceptions import BlastError


class BlastRunner:
    """BLAST execution wrapper."""
    
    def __init__(self, reference: Path, threads: int = 1):
        """
        Initialize BLAST runner.
        
        Args:
            reference: Path to BLAST database
            threads: Number of threads to use
        """
        self.reference = Path(reference)
        self.threads = threads
    
    def run(self, query_file: Path, output_file: Path) -> Path:
        """
        Execute BLAST search.
        
        Args:
            query_file: Input FASTA file
            output_file: Output file path
            
        Returns:
            Path to output file
            
        Raises:
            BlastError: If BLAST execution fails
        """
        cmd = [
            "blastn",
            "-task", "blastn",  # V2 compatibility
            "-query", str(query_file),
            "-db", str(self.reference),
            "-out", str(output_file),
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq slen",
            "-word_size", "11",  # V2 compatibility - important for finding all homeologs
            "-num_threads", str(self.threads),
            "-max_target_seqs", "20"
        ]
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            return Path(output_file)
        except subprocess.CalledProcessError as e:
            raise BlastError(f"BLAST execution failed: {e.stderr}") from e
        except FileNotFoundError as e:
            raise BlastError(f"BLAST not found in PATH: {e}") from e


class BlastParser:
    """BLAST output parser."""
    
    def __init__(self, blast_file: Path):
        """
        Initialize parser.
        
        Args:
            blast_file: Path to BLAST output file
        """
        self.blast_file = Path(blast_file)
        self.hits: Dict[str, List[BlastHit]] = {}
    
    def parse(self) -> Dict[str, List[BlastHit]]:
        """
        Parse BLAST output file.
        
        Returns:
            Dictionary mapping query IDs to lists of BlastHit objects
            
        Raises:
            BlastError: If parsing fails
        """
        try:
            with open(self.blast_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    fields = line.split('\t')
                    if len(fields) < 15:
                        continue
                    
                    try:
                        hit = BlastHit(
                            query_id=fields[0],
                            subject_id=fields[1],
                            identity=float(fields[2]),
                            alignment_length=int(fields[3]),
                            mismatches=int(fields[4]),
                            gap_opens=int(fields[5]),
                            query_start=int(fields[6]),
                            query_end=int(fields[7]),
                            subject_start=int(fields[8]),
                            subject_end=int(fields[9]),
                            evalue=float(fields[10]),
                            bit_score=float(fields[11]),
                            query_seq=fields[12],
                            subject_seq=fields[13],
                            subject_length=int(fields[14])
                        )
                        
                        if hit.query_id not in self.hits:
                            self.hits[hit.query_id] = []
                        self.hits[hit.query_id].append(hit)
                        
                    except (ValueError, IndexError) as e:
                        # Skip malformed lines
                        continue
            
            return self.hits
            
        except IOError as e:
            raise BlastError(f"Failed to read BLAST file {self.blast_file}: {e}") from e


class FlankingExtractor:
    """Extract flanking regions from BLAST hits."""
    
    def __init__(self, reference: Path):
        """
        Initialize extractor.
        
        Args:
            reference: Path to BLAST database
        """
        self.reference = Path(reference)
    
    def extract_flanking_regions(
        self,
        blast_hits: Dict[str, List[BlastHit]],
        snp_positions: Dict[str, int],
        flanking_size: int = 500,
        max_hits: int = 6,
        min_identity: float = 88.0,
        min_alignment_ratio: float = 0.9
    ) -> List[FlankingRegion]:
        """
        Extract flanking regions from BLAST hits.
        
        Args:
            blast_hits: Dictionary of BLAST hits by query ID
            snp_positions: Dictionary mapping SNP names to positions in query
            flanking_size: Size of flanking regions to extract
            max_hits: Maximum number of hits per SNP
            min_identity: Minimum identity percentage
            min_alignment_ratio: Minimum alignment length ratio
            
        Returns:
            List of FlankingRegion objects
        """
        flanking_regions = []
        snp_hit_counts = Counter()
        top_hits = {}  # Best hit for each query
        
        # First pass: collect valid hits and count
        valid_hits = []
        
        for query_id, hits in blast_hits.items():
            snp_name, qchrom, allele = self._parse_query_id(query_id)
            
            if snp_name not in snp_positions:
                continue
                
            # Convert 0-based SNP position to 1-based for BLAST coordinate comparison
            snp_pos_0based = snp_positions[snp_name]
            snp_pos_1based = snp_pos_0based + 1
            
            for hit in hits:
                # Calculate adjusted identity (accounting for gaps)
                pct_identity = 100 - (hit.mismatches + hit.gap_opens) / hit.alignment_length * 100
                
                # Set minimum alignment length
                min_align_length = max(50, hit.alignment_length * min_alignment_ratio)
                
                # Filter by identity and alignment length
                if pct_identity <= min_identity or hit.alignment_length <= min_align_length:
                    continue
                
                # Check if SNP is within alignment (using 1-based coordinates)
                if snp_pos_1based < hit.query_start or snp_pos_1based > hit.query_end:
                    continue
                
                # Calculate SNP position in subject considering gaps
                subject_snp_pos = self._calculate_subject_snp_position(
                    hit, snp_pos_1based
                )
                
                if subject_snp_pos is None:
                    continue
                
                # Create flanking region
                region = self._create_flanking_region(
                    snp_name, hit, subject_snp_pos, flanking_size, allele
                )
                
                valid_hits.append((query_id, hit, region))
                snp_hit_counts[query_id] += 1
                
                # Track top hit (first hit is best due to BLAST sorting)
                if query_id not in top_hits:
                    top_hits[query_id] = (hit, region)
        
        # Second pass: select regions based on hit count and target chromosome
        target_regions = {}  # query_id -> region for target chromosome hits
        
        for query_id, hit, region in valid_hits:
            if snp_hit_counts[query_id] > max_hits:
                continue
            
            snp_name, qchrom, allele = self._parse_query_id(query_id)
            
            # Prefer hits on target chromosome
            if qchrom == hit.subject_id:
                target_regions[query_id] = region
            
            flanking_regions.append(region)
        
        # Handle SNPs with no target chromosome hits
        for query_id in snp_hit_counts:
            if (snp_hit_counts[query_id] <= max_hits and 
                query_id not in target_regions and 
                query_id in top_hits):
                
                print(f"Warning: no hits on target chromosome for {query_id}. Using best hit.")
                # The region was already added in the first pass
        
        return flanking_regions
    
    def _parse_query_id(self, query_id: str) -> Tuple[str, str, str]:
        """Parse query ID into SNP name, chromosome, and allele."""
        parts = query_id.split("_")
        if len(parts) >= 3:
            snp_name = "_".join(parts[:-2])
            chrom = parts[-2]
            allele = parts[-1]
        else:
            # Fallback for malformed IDs
            snp_name = query_id
            chrom = "unknown"
            allele = "unknown"
        
        return snp_name, chrom, allele
    
    def _calculate_subject_snp_position(self, hit: BlastHit, snp_pos: int) -> Optional[int]:
        """
        Calculate SNP position in subject sequence considering gaps.
        
        Args:
            hit: BLAST hit object
            snp_pos: SNP position in query (1-based)
            
        Returns:
            SNP position in subject sequence, or None if calculation fails
        """
        # Distance from query start to SNP (0-based)
        temp = snp_pos - hit.query_start
        
        # Count gaps before SNP position
        qgaps = [i + 1 for i, char in enumerate(hit.query_seq) if char == "-"]
        sgaps = [i + 1 for i, char in enumerate(hit.subject_seq) if char == "-"]
        
        # Adjust for query gaps
        nqgap = sum(1 for gap_pos in qgaps if gap_pos <= temp)
        temp += nqgap
        
        # Adjust for subject gaps
        nsgap = sum(1 for gap_pos in sgaps if gap_pos <= temp)
        
        # Calculate position in subject
        if hit.strand == Strand.PLUS:
            pos = hit.subject_start + (temp - nsgap)
        else:
            pos = hit.subject_start - (temp - nsgap)
        
        return pos
    
    def _create_flanking_region(
        self,
        snp_name: str,
        hit: BlastHit,
        subject_snp_pos: int,
        flanking_size: int,
        allele: str
    ) -> FlankingRegion:
        """Create FlankingRegion object from hit and SNP position."""
        
        if hit.strand == Strand.PLUS:
            start = max(1, subject_snp_pos - flanking_size)
            end = min(hit.subject_length, subject_snp_pos + flanking_size)
            snp_pos_in_region = subject_snp_pos - start + 1
        else:
            start = max(1, subject_snp_pos - flanking_size)
            end = min(hit.subject_length, subject_snp_pos + flanking_size)
            snp_pos_in_region = end - subject_snp_pos + 1
        
        return FlankingRegion(
            snp_name=snp_name,
            chromosome=hit.subject_id,
            start=start,
            end=end,
            strand=hit.strand,
            snp_position_in_region=snp_pos_in_region,
            allele=allele
        )
    
    def extract_sequences(self, regions: List[FlankingRegion], output_file: Path) -> Path:
        """
        Extract sequences using blastdbcmd.
        
        Args:
            regions: List of flanking regions to extract
            output_file: Output FASTA file path
            
        Returns:
            Path to output file
            
        Raises:
            BlastError: If sequence extraction fails
        """
        try:
            with open(output_file, 'w') as f:
                for region in regions:
                    # Format sequence ID
                    seq_id = f"{region.snp_name}_{region.chromosome}_{region.allele}_{region.snp_position_in_region}"
                    
                    # Extract sequence using blastdbcmd
                    cmd = [
                        "blastdbcmd",
                        "-db", str(self.reference),
                        "-entry", region.chromosome,
                        "-range", f"{region.start}-{region.end}",
                        "-strand", "plus" if region.strand == Strand.PLUS else "minus"
                    ]
                    
                    result = subprocess.run(
                        cmd,
                        capture_output=True,
                        text=True,
                        check=True
                    )
                    
                    # Parse output and write to file
                    lines = result.stdout.strip().split('\n')
                    if len(lines) >= 2:
                        sequence = ''.join(lines[1:])  # Skip header line
                        f.write(f">{seq_id}\n{sequence}\n")
            
            return Path(output_file)
            
        except subprocess.CalledProcessError as e:
            raise BlastError(f"Sequence extraction failed: {e.stderr}") from e
        except IOError as e:
            raise BlastError(f"Failed to write output file {output_file}: {e}") from e