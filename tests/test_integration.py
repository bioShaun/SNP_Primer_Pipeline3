#!/usr/bin/env python3
"""
Integration tests for SNP Primer Pipeline.

These tests verify that the main components work together correctly.
"""

import pytest
import tempfile
from pathlib import Path

from snp_primer_pipeline import (
    PolymarkerParser, 
    KASPDesigner, 
    CAPSDesigner,
    SNP
)


def test_polymarker_parser_integration():
    """Test that PolymarkerParser can parse a simple input."""
    # Create a temporary input file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("SNP1,chr1A,ATCGATCG[A/G]TCGATCGA\n")
        f.write("SNP2,chr2B,GCTAGCTA[C/T]AGCTAGCT\n")
        temp_file = Path(f.name)
    
    try:
        parser = PolymarkerParser(temp_file)
        snps = parser.parse()
        
        assert len(snps) == 2
        assert snps[0].name == "SNP1"
        assert snps[0].chromosome == "chr1A"
        assert snps[0].allele_a == "A"
        assert snps[0].allele_b == "G"
        assert snps[0].iupac_code == "R"
        
        assert snps[1].name == "SNP2"
        assert snps[1].chromosome == "chr2B"
        assert snps[1].allele_a == "C"
        assert snps[1].allele_b == "T"
        assert snps[1].iupac_code == "Y"
        
    finally:
        temp_file.unlink()


def test_kasp_designer_initialization():
    """Test that KASPDesigner can be initialized."""
    designer = KASPDesigner()
    assert designer is not None
    assert designer.max_tm == 63.0
    assert designer.max_size == 25


def test_caps_designer_initialization():
    """Test that CAPSDesigner can be initialized."""
    designer = CAPSDesigner()
    assert designer is not None
    assert designer.max_tm == 63.0
    assert designer.max_size == 25


def test_caps_designer_with_enzyme_file():
    """Test that CAPSDesigner can load enzyme database."""
    # Create a temporary enzyme file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("EcoRI,50\tGAATTC\n")
        f.write("BamHI,75\tGGATCC\n")
        temp_file = Path(f.name)
    
    try:
        designer = CAPSDesigner(enzyme_file=temp_file)
        designer.load_enzymes()
        
        assert len(designer.enzymes) == 2
        assert "EcoRI" in designer.enzymes
        assert "BamHI" in designer.enzymes
        assert designer.enzymes["EcoRI"].sequence == "GAATTC"
        assert designer.enzymes["EcoRI"].price == 50
        
    finally:
        temp_file.unlink()


def test_snp_model_serialization():
    """Test SNP model serialization and deserialization."""
    snp = SNP(
        name="test_snp",
        chromosome="chr1",
        flanking_sequence="ATCGATCGATCG[A/G]TCGATCGATCG",
        snp_position=12,
        allele_a="A",
        allele_b="G",
        iupac_code="R"
    )
    
    # Test serialization
    snp_dict = snp.to_dict()
    assert snp_dict["name"] == "test_snp"
    assert snp_dict["chromosome"] == "chr1"
    assert snp_dict["allele_a"] == "A"
    assert snp_dict["allele_b"] == "G"
    
    # Test deserialization
    snp2 = SNP.from_dict(snp_dict)
    assert snp2.name == snp.name
    assert snp2.chromosome == snp.chromosome
    assert snp2.allele_a == snp.allele_a
    assert snp2.allele_b == snp.allele_b


def test_package_imports():
    """Test that all main components can be imported."""
    from snp_primer_pipeline import (
        __version__,
        PipelineConfig,
        SNP,
        PolymarkerParser,
        KASPDesigner,
        CAPSDesigner,
        run_pipeline
    )
    
    assert __version__ == "3.0.0"
    assert PipelineConfig is not None
    assert SNP is not None
    assert PolymarkerParser is not None
    assert KASPDesigner is not None
    assert CAPSDesigner is not None
    assert run_pipeline is not None