"""Basic tests to verify project structure."""

import pytest
from pathlib import Path

def test_package_imports():
    """Test that basic package imports work."""
    from snp_primer_pipeline import __version__, PipelineConfig, SNP
    
    assert __version__ == "3.0.0"
    assert PipelineConfig is not None
    assert SNP is not None

def test_config_creation():
    """Test basic config creation."""
    from snp_primer_pipeline.config import PipelineConfig
    
    # Create a temporary input file for testing
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("test,chr1,ATCG[A/G]ATCG\n")
        temp_input = Path(f.name)
    
    try:
        config = PipelineConfig(
            input_file=temp_input,
            reference_file=Path("dummy.fasta")  # Won't be validated in basic test
        )
        assert config.input_file == temp_input
        assert config.max_tm == 63.0
        assert config.design_kasp is True
    finally:
        temp_input.unlink()

def test_exceptions():
    """Test that custom exceptions work."""
    from snp_primer_pipeline.exceptions import PipelineError, ParseError
    
    with pytest.raises(PipelineError):
        raise PipelineError("Test error")
    
    with pytest.raises(ParseError):
        raise ParseError("Parse error", line_number=1)

def test_models():
    """Test basic model functionality."""
    from snp_primer_pipeline.models import SNP, Primer
    
    snp = SNP(
        name="test_snp",
        chromosome="chr1", 
        flanking_sequence="ATCGATCG",
        snp_position=4,
        allele_a="A",
        allele_b="G",
        iupac_code="R"
    )
    
    assert snp.name == "test_snp"
    assert snp.alleles == ("A", "G")
    
    # Test serialization
    snp_dict = snp.to_dict()
    snp2 = SNP.from_dict(snp_dict)
    assert snp == snp2