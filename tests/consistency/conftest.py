#!/usr/bin/env python3
"""
Shared fixtures for consistency tests.
"""

import pytest
from pathlib import Path
import tempfile
import shutil
from typing import Dict, Any

from ..utils.reference_loader import ReferenceDataLoader
from ..utils.output_comparator import OutputComparator
from ..utils.consistency_reporter import ConsistencyReporter


@pytest.fixture(scope="session")
def examples_dir():
    """SNP_Primer_Pipeline2 examples directory."""
    examples_path = Path(__file__).parent.parent.parent.parent / "SNP_Primer_Pipeline2" / "examples"
    if not examples_path.exists():
        pytest.skip(f"SNP_Primer_Pipeline2 examples directory not found: {examples_path}")
    return examples_path


@pytest.fixture(scope="session")
def reference_loader(examples_dir):
    """Reference data loader."""
    return ReferenceDataLoader(examples_dir)


@pytest.fixture
def comparator():
    """Output comparator."""
    return OutputComparator()


@pytest.fixture
def reporter(tmp_path):
    """Consistency reporter."""
    return ConsistencyReporter(tmp_path)


@pytest.fixture(scope="session")
def v3_pipeline_config():
    """V3 pipeline configuration."""
    return {
        'blast_params': {
            'evalue': 1e-5,
            'word_size': 11,
            'max_target_seqs': 1000
        },
        'primer3_params': {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 5,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_END': 3,
            'PRIMER_PAIR_MAX_COMPL_ANY': 8,
            'PRIMER_PAIR_MAX_DIFF_TM': 4,
            'PRIMER_MAX_END_STABILITY': 9,
            'PRIMER_MAX_SELF_ANY': 8,
            'PRIMER_MAX_SELF_END': 3,
            'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING': 12,
            'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH': 40,
            'PRIMER_MAX_TEMPLATE_MISPRIMING': 12,
            'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': 40,
            'PRIMER_MIN_END_QUALITY': 0,
            'PRIMER_MIN_QUALITY': 0,
            'PRIMER_QUALITY_RANGE_MIN': 0,
            'PRIMER_QUALITY_RANGE_MAX': 100,
            'PRIMER_PRODUCT_SIZE_RANGE': [[75, 1000]],
            'PRIMER_NUM_RETURN': 5
        },
        'alignment_params': {
            'gap_open': 10,
            'gap_extend': 0.5,
            'match': 2,
            'mismatch': -1
        }
    }


@pytest.fixture
def v3_pipeline_runner(v3_pipeline_config, tmp_path):
    """V3 pipeline runner function."""
    def runner(input_file: Path, reference_db: Path, output_dir: Path = None, **kwargs):
        """Run V3 pipeline with given inputs."""
        if output_dir is None:
            output_dir = tmp_path / "v3_output"
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Import V3 pipeline components
        try:
            from snp_primer_pipeline.config import PipelineConfig
            from snp_primer_pipeline.main import run_pipeline
        except ImportError:
            pytest.skip("SNP_Primer_Pipeline3 not properly installed")
        
        # Merge config with any overrides
        config_dict = {**v3_pipeline_config, **kwargs}
        
        # Create config object
        config = PipelineConfig(
            input_file=str(input_file),
            reference_file=str(reference_db),
            output_dir=str(output_dir),
            **config_dict
        )
        
        # Run pipeline
        run_pipeline(config)
        
        return output_dir
    
    return runner


@pytest.fixture(scope="session")
def reference_kasp_data(reference_loader):
    """KASP primer reference data."""
    return reference_loader.load_kasp_primers()


@pytest.fixture(scope="session")
def reference_caps_data(reference_loader):
    """CAPS primer reference data."""
    return reference_loader.load_caps_primers()


@pytest.fixture(scope="session")
def reference_blast_data(reference_loader):
    """BLAST results reference data."""
    return reference_loader.load_blast_results()


@pytest.fixture(scope="session")
def reference_polymarker_input(reference_loader):
    """Polymarker input reference data."""
    return reference_loader.load_polymarker_input()


@pytest.fixture(scope="session")
def snp_names(reference_loader):
    """List of SNP names from reference data."""
    return reference_loader.get_snp_names()


@pytest.fixture
def test_reference_db(examples_dir, tmp_path):
    """Create a test reference database."""
    # Look for reference files in examples
    ref_files = list(examples_dir.glob("*.fa")) + list(examples_dir.glob("*.fasta"))
    
    if not ref_files:
        # Create a minimal test reference if none found
        test_ref = tmp_path / "test_reference.fa"
        with open(test_ref, 'w') as f:
            f.write(">chr7A\n")
            f.write("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
        return test_ref
    
    # Copy the first reference file found
    test_ref = tmp_path / "test_reference.fa"
    shutil.copy2(ref_files[0], test_ref)
    return test_ref


@pytest.fixture
def test_input_file(examples_dir, tmp_path):
    """Create a test input file."""
    input_file = examples_dir / "polymarker_input_example.csv"
    
    if input_file.exists():
        test_input = tmp_path / "test_input.csv"
        shutil.copy2(input_file, test_input)
        return test_input
    
    # Create a minimal test input if none found
    test_input = tmp_path / "test_input.csv"
    with open(test_input, 'w') as f:
        f.write("chr7A-7659,chr7A,ATCGATCGATCGATC[T/C]GATCGATCGATCGATC\n")
        f.write("chr7A-7716,chr7A,GATCGATCGATCGAT[A/G]TCGATCGATCGATCGA\n")
    
    return test_input


@pytest.fixture
def consistency_test_data(examples_dir):
    """Comprehensive test data for consistency tests."""
    return {
        'examples_dir': examples_dir,
        'kasp_file': examples_dir / "Potential_KASP_primers.tsv",
        'caps_file': examples_dir / "Potential_CAPS_primers.tsv",
        'blast_file': examples_dir / "blast_out.txt",
        'input_file': examples_dir / "polymarker_input_example.csv",
        'alignment_files': list(examples_dir.glob("alignment_raw_*.fa")),
        'flanking_files': list(examples_dir.glob("flanking_temp_marker_*.fa"))
    }


@pytest.fixture(autouse=True)
def setup_test_environment(tmp_path, monkeypatch):
    """Set up test environment."""
    # Change to temp directory for test isolation
    monkeypatch.chdir(tmp_path)
    
    # Set environment variables if needed
    monkeypatch.setenv("PRIMER3_CONFIG_DIR", str(tmp_path))
    
    yield
    
    # Cleanup is automatic with tmp_path


@pytest.fixture
def mock_blast_db(tmp_path):
    """Create a mock BLAST database for testing."""
    db_path = tmp_path / "test_db"
    
    # Create minimal BLAST database files
    for ext in ['.nhr', '.nin', '.nsq']:
        (db_path.parent / f"{db_path.name}{ext}").touch()
    
    return db_path