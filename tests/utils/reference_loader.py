#!/usr/bin/env python3
"""
Reference data loader for consistency tests.

Loads reference output data from SNP_Primer_Pipeline2/examples directory.
"""

from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import csv
import re


@dataclass
class KASPPrimerRecord:
    """KASP primer record from reference data."""
    index: str
    product_size: int
    primer_type: str  # LEFT, RIGHT
    start: int
    end: int
    variation_number: int
    diff_three_all: str  # YES/NO
    length: int
    tm: float
    gc_content: float
    self_any: float
    self_three: float
    end_stability: float
    hairpin: float
    primer_seq: str
    reverse_complement: str
    penalty: float
    compl_any: float
    compl_end: float
    score: float


@dataclass
class CAPSPrimerRecord:
    """CAPS primer record from reference data."""
    index: str
    product_size: int
    primer_type: str
    start: int
    end: int
    diff_number: int
    diff_three_all: str
    length: int
    tm: float
    gc_content: float
    self_any: float
    self_three: float
    end_stability: float
    hairpin: float
    primer_seq: str
    reverse_complement: str
    penalty: float
    compl_any: float
    compl_end: float


@dataclass
class BlastHitRecord:
    """BLAST hit record from reference data."""
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
    query_seq: str = ""
    subject_seq: str = ""
    subject_length: int = 0


class ReferenceDataLoader:
    """Loads reference data from SNP_Primer_Pipeline2/examples."""
    
    def __init__(self, examples_dir: Path):
        self.examples_dir = Path(examples_dir)
        if not self.examples_dir.exists():
            raise FileNotFoundError(f"Examples directory not found: {examples_dir}")
    
    def load_kasp_primers(self) -> Dict[str, List[KASPPrimerRecord]]:
        """Load KASP primer reference data."""
        kasp_file = self.examples_dir / "Potential_KASP_primers.tsv"
        if not kasp_file.exists():
            raise FileNotFoundError(f"KASP primers file not found: {kasp_file}")
        
        primers_by_snp = {}
        current_snp = None
        
        with open(kasp_file, 'r') as f:
            lines = f.readlines()
        
        # Skip header and process data
        in_data_section = False
        for line in lines:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith("index\tproduct_size"):
                in_data_section = True
                continue
            
            if line.startswith("Sites that can differ all"):
                in_data_section = False
                continue
                
            if not in_data_section:
                continue
            
            fields = line.split('\t')
            if len(fields) < 20:
                continue
            
            try:
                record = KASPPrimerRecord(
                    index=fields[0],
                    product_size=int(fields[1]),
                    primer_type=fields[2],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    variation_number=int(fields[5]),
                    diff_three_all=fields[6],
                    length=int(fields[7]),
                    tm=float(fields[8]),
                    gc_content=float(fields[9]),
                    self_any=float(fields[10]),
                    self_three=float(fields[11]),
                    end_stability=float(fields[12]),
                    hairpin=float(fields[13]),
                    primer_seq=fields[14],
                    reverse_complement=fields[15],
                    penalty=float(fields[16]),
                    compl_any=float(fields[17]),
                    compl_end=float(fields[18]),
                    score=float(fields[19])
                )
                
                # Extract SNP name from index
                snp_name = record.index.split('-')[0] + '-' + record.index.split('-')[1]
                if snp_name not in primers_by_snp:
                    primers_by_snp[snp_name] = []
                primers_by_snp[snp_name].append(record)
                
            except (ValueError, IndexError) as e:
                print(f"Warning: Failed to parse KASP line: {line[:50]}... Error: {e}")
                continue
        
        return primers_by_snp
    
    def load_caps_primers(self) -> Dict[str, List[CAPSPrimerRecord]]:
        """Load CAPS primer reference data."""
        caps_file = self.examples_dir / "Potential_CAPS_primers.tsv"
        if not caps_file.exists():
            raise FileNotFoundError(f"CAPS primers file not found: {caps_file}")
        
        primers_by_snp = {}
        
        with open(caps_file, 'r') as f:
            lines = f.readlines()
        
        # Skip header and process data
        in_data_section = False
        for line in lines:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith("index\tproduct_size"):
                in_data_section = True
                continue
            
            if line.startswith("Sites that can differ all"):
                in_data_section = False
                continue
                
            if not in_data_section:
                continue
            
            fields = line.split('\t')
            if len(fields) < 19:
                continue
            
            try:
                record = CAPSPrimerRecord(
                    index=fields[0],
                    product_size=int(fields[1]),
                    primer_type=fields[2],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    diff_number=int(fields[5]),
                    diff_three_all=fields[6],
                    length=int(fields[7]),
                    tm=float(fields[8]),
                    gc_content=float(fields[9]),
                    self_any=float(fields[10]),
                    self_three=float(fields[11]),
                    end_stability=float(fields[12]),
                    hairpin=float(fields[13]),
                    primer_seq=fields[14],
                    reverse_complement=fields[15],
                    penalty=float(fields[16]),
                    compl_any=float(fields[17]),
                    compl_end=float(fields[18])
                )
                
                # Extract SNP name from index
                snp_name = record.index.split('-')[0] + '-' + record.index.split('-')[1]
                if snp_name not in primers_by_snp:
                    primers_by_snp[snp_name] = []
                primers_by_snp[snp_name].append(record)
                
            except (ValueError, IndexError) as e:
                print(f"Warning: Failed to parse CAPS line: {line[:50]}... Error: {e}")
                continue
        
        return primers_by_snp
    
    def load_blast_results(self) -> Dict[str, List[BlastHitRecord]]:
        """Load BLAST results reference data."""
        blast_file = self.examples_dir / "blast_out.txt"
        if not blast_file.exists():
            raise FileNotFoundError(f"BLAST results file not found: {blast_file}")
        
        hits_by_query = {}
        
        with open(blast_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                fields = line.split('\t')
                if len(fields) < 12:
                    continue
                
                try:
                    record = BlastHitRecord(
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
                        query_seq=fields[12] if len(fields) > 12 else "",
                        subject_seq=fields[13] if len(fields) > 13 else "",
                        subject_length=int(fields[14]) if len(fields) > 14 else 0
                    )
                    
                    if record.query_id not in hits_by_query:
                        hits_by_query[record.query_id] = []
                    hits_by_query[record.query_id].append(record)
                    
                except (ValueError, IndexError) as e:
                    print(f"Warning: Failed to parse BLAST line: {line[:50]}... Error: {e}")
                    continue
        
        return hits_by_query
    
    def load_alignment(self, snp_name: str) -> Dict[str, str]:
        """Load multiple sequence alignment reference data."""
        alignment_file = self.examples_dir / f"alignment_raw_{snp_name}.fa"
        if not alignment_file.exists():
            return {}
        
        sequences = {}
        current_seq_name = None
        
        with open(alignment_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    current_seq_name = line[1:].split()[0]
                    sequences[current_seq_name] = ""
                elif current_seq_name:
                    sequences[current_seq_name] += line
        
        return sequences
    
    def load_flanking_sequences(self, snp_name: str) -> Dict[str, str]:
        """Load flanking sequences reference data."""
        # Look for flanking files matching the pattern
        flanking_files = list(self.examples_dir.glob(f"flanking_temp_marker_{snp_name}_*.fa"))
        
        sequences = {}
        for flanking_file in flanking_files:
            current_seq_name = None
            with open(flanking_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        current_seq_name = line[1:].split()[0]
                        sequences[current_seq_name] = ""
                    elif current_seq_name:
                        sequences[current_seq_name] += line
        
        return sequences
    
    def load_variation_sites(self, snp_name: str) -> List[int]:
        """Extract variation sites from KASP output file."""
        kasp_file = self.examples_dir / "Potential_KASP_primers.tsv"
        if not kasp_file.exists():
            return []
        
        with open(kasp_file, 'r') as f:
            content = f.read()
        
        # Look for the variation sites section
        pattern = f"Sites that can differ all for {snp_name}\\s*\\n([0-9, ]+)"
        match = re.search(pattern, content)
        
        if match:
            sites_str = match.group(1).strip()
            try:
                # Parse comma-separated numbers
                sites = [int(x.strip()) for x in sites_str.split(',') if x.strip()]
                return sites
            except ValueError:
                return []
        
        return []
    
    def get_snp_names(self) -> List[str]:
        """Get all SNP names from the reference data."""
        snp_names = set()
        
        # Extract from KASP primers
        kasp_primers = self.load_kasp_primers()
        snp_names.update(kasp_primers.keys())
        
        # Extract from CAPS primers
        caps_primers = self.load_caps_primers()
        snp_names.update(caps_primers.keys())
        
        return sorted(list(snp_names))
    
    def load_polymarker_input(self) -> List[Dict[str, str]]:
        """Load polymarker input file."""
        input_file = self.examples_dir / "polymarker_input_example.csv"
        if not input_file.exists():
            raise FileNotFoundError(f"Polymarker input file not found: {input_file}")
        
        snps = []
        with open(input_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split(',')
                if len(parts) >= 3:
                    snps.append({
                        'name': parts[0],
                        'chromosome': parts[1],
                        'sequence': parts[2]
                    })
        
        return snps