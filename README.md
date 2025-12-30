# SNP Primer Pipeline 3 Claude

A modern, modular Python pipeline for designing KASP and CAPS/dCAPS primers for SNP genotyping in any species. This is a refactored version of SNP_Primer_Pipeline2 with improved architecture, testability, and maintainability.

## Features

- **KASP Primer Design**: Design allele-specific primers for KASP genotyping
- **CAPS/dCAPS Primer Design**: Design primers for restriction enzyme-based genotyping
- **Multiple Species Support**: Works with any species with a reference genome
- **Modular Architecture**: Clean separation of concerns with testable components
- **Modern Python**: Uses Python 3.10+ with type hints and dataclasses
- **Comprehensive Testing**: Property-based testing with hypothesis
- **CLI Interface**: Easy-to-use command-line interface

## Installation

### Prerequisites

- Python 3.10 or higher
- BLAST+ (blastn, blastdbcmd)
- MUSCLE or MAFFT for multiple sequence alignment
- Primer3 for primer design

### Install from source

```bash
git clone <repository-url>
cd SNP_Primer_Pipeline3_claude
pip install -e .
```

## Quick Start

### Basic Usage

```bash
# Design both KASP and CAPS primers
python -m snp_primer_pipeline.main input.csv reference_db -o output_dir

# Design only KASP primers
python -m snp_primer_pipeline.main input.csv reference_db --no-caps

# Design only CAPS primers with custom enzyme price limit
python -m snp_primer_pipeline.main input.csv reference_db --no-kasp --max-price 150
```

### Input Format

The input file should be a CSV file with the following format:

```
SNP_name,chromosome,flanking_sequence
SNP1,chr1A,ATCGATCGATCG[A/G]TCGATCGATCG
SNP2,chr2B,GCTAGCTAGCTA[C/T]AGCTAGCTAGT
```

Where:
- `SNP_name`: Unique identifier for the SNP
- `chromosome`: Target chromosome name
- `flanking_sequence`: DNA sequence with SNP in IUPAC bracket notation

### Python API

```python
from snp_primer_pipeline import (
    PipelineConfig, 
    run_pipeline,
    PolymarkerParser,
    KASPDesigner,
    CAPSDesigner
)

# Create configuration
config = PipelineConfig(
    input_file="input.csv",
    reference_file="reference_db",
    output_dir="output",
    design_kasp=True,
    design_caps=True
)

# Run complete pipeline
run_pipeline(config)

# Or use individual components
parser = PolymarkerParser("input.csv")
snps = parser.parse()

kasp_designer = KASPDesigner()
primers = kasp_designer.design_primers(
    template_sequence="ATCGATCG...",
    snp_position=100,
    snp_alleles=("A", "G")
)
```

## Architecture

The pipeline is organized into several modules:

### Core Modules (`snp_primer_pipeline.core`)

- **`parser.py`**: Parse polymarker input files and convert IUPAC codes
- **`blast.py`**: BLAST execution, result parsing, and flanking region extraction
- **`alignment.py`**: Multiple sequence alignment and variant site identification
- **`primer3_parser.py`**: Primer3 interface for primer design

### Primer Design Modules (`snp_primer_pipeline.primers`)

- **`kasp.py`**: KASP primer design with allele-specific primers
- **`caps.py`**: CAPS/dCAPS primer design with restriction enzyme analysis

### Data Models (`snp_primer_pipeline.models`)

- **`SNP`**: SNP information with alleles and position
- **`BlastHit`**: BLAST alignment results
- **`FlankingRegion`**: Genomic regions around SNPs
- **`Primer`**: Individual primer properties
- **`PrimerPair`**: Primer pair with scoring
- **`RestrictionEnzyme`**: Enzyme properties for CAPS design

### Configuration (`snp_primer_pipeline.config`)

- **`PipelineConfig`**: Main pipeline configuration
- **`SoftwarePaths`**: External software path management

## Output

The pipeline generates the following outputs:

### KASP Primers
- `KASP_primers_{SNP_name}.txt`: Designed KASP primers with scores
- Allele-specific primers for each SNP allele
- Common primers for amplification

### CAPS Primers
- `CAPS_primers_{SNP_name}.txt`: Designed CAPS/dCAPS primers
- Enzyme information and cutting sites
- Primer pairs for restriction enzyme assays

### Intermediate Files
- `for_blast.fa`: FASTA file for BLAST search
- `blast_out.txt`: BLAST results
- `flanking_sequences.fa`: Extracted flanking sequences
- `alignment_{SNP_name}.fa`: Multiple sequence alignments

## Configuration Options

### Command Line Options

```
--max-tm FLOAT          Maximum primer Tm (default: 63.0)
--max-size INT          Maximum primer size (default: 25)
--max-price INT         Maximum enzyme price for CAPS (default: 200)
--pick-anyway           Pick primers even if constraints violated
--threads INT           Number of BLAST threads (default: 1)
--log-level LEVEL       Logging level (DEBUG/INFO/WARNING/ERROR)
```

### Configuration File

Create a `config.yaml` file for advanced configuration:

```yaml
# Pipeline settings
flanking_size: 500
max_hits: 6
primer_product_size_range: [50, 250]

# Primer design settings
max_tm: 63.0
max_primer_size: 25
pick_anyway: false

# CAPS settings
max_price: 200

# Software paths (auto-detected if not specified)
primer3_path: "/usr/local/bin/primer3_core"
muscle_path: "/usr/local/bin/muscle"
```

## Testing

The pipeline includes comprehensive tests:

```bash
# Run all tests
pytest tests/

# Run with coverage
pytest tests/ --cov=snp_primer_pipeline --cov-report=html

# Run property-based tests
pytest tests/ -k "property"
```

## Development

### Project Structure

```
SNP_Primer_Pipeline3_claude/
├── src/snp_primer_pipeline/     # Main package
│   ├── core/                    # Core processing modules
│   ├── primers/                 # Primer design modules
│   ├── config.py               # Configuration management
│   ├── models.py               # Data models
│   ├── exceptions.py           # Custom exceptions
│   └── main.py                 # CLI interface
├── tests/                      # Test suite
├── resources/                  # Resource files
│   └── NEB_parsed_REs.txt     # Enzyme database
├── pyproject.toml             # Project configuration
└── README.md                  # This file
```

### Adding New Features

1. **New Primer Types**: Extend the `primers` module
2. **New Input Formats**: Extend the `parser` module
3. **New Alignment Tools**: Extend the `alignment` module
4. **New Data Models**: Add to `models.py`

### Code Quality

The project follows modern Python best practices:

- Type hints throughout
- Dataclasses for data models
- Property-based testing with hypothesis
- Comprehensive error handling
- Modular, testable architecture

## Comparison with SNP_Primer_Pipeline2

### Improvements

1. **Modern Python**: Uses Python 3.10+ features
2. **Modular Design**: Clean separation of concerns
3. **Type Safety**: Full type hints for better IDE support
4. **Testability**: Comprehensive test suite with property-based testing
5. **Error Handling**: Robust error handling with custom exceptions
6. **Documentation**: Comprehensive documentation and examples
7. **Maintainability**: Clean code structure following best practices

### Migration Guide

For users of SNP_Primer_Pipeline2:

1. **Input Format**: Same polymarker CSV format
2. **Output Format**: Similar but with improved organization
3. **Dependencies**: Requires Python 3.10+ instead of Python 2/3
4. **Installation**: Uses modern pip/pyproject.toml instead of manual setup
5. **Configuration**: YAML configuration instead of hardcoded values

## License

This project is licensed under the GNU General Public License v2.0 - see the original SNP_Primer_Pipeline2 for details.

## Citation

If you use this pipeline in your research, please cite:

> Zhang, J. et al. SNP Primer Pipeline: A modern Python pipeline for designing KASP and CAPS primers for SNP genotyping.

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## Support

For questions and support:

1. Check the documentation
2. Review existing issues
3. Create a new issue with detailed information
4. Include input files and error messages when reporting bugs