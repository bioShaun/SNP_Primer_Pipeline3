# Implementation Plan: Support for SNP Position Input Format

## Objective
Enable the pipeline to accept simple coordinate-based SNP input files (`Chromosome`, `Position`, `Ref`, `Alt`) by automatically fetching flanking sequences from the reference genome.

## Input Format Comparison
| Format | Example | Fields |
|--------|---------|--------|
| **Polymarker CSV** | `chr7A-7659,chr7A,ATCG...[T/C]...GCTA` | SNP_ID, Chromosome, Sequence |
| **Coordinate (NEW)** | `chr7A	7659	T	C` | Chrom, Pos, Ref, Alt (tab/space) |

## Architecture Decision

### Approach: Unified Parser with Format Auto-Detection

Instead of creating a separate `CoordinateParser` class, we will:
1. **Extend `PolymarkerParser`** with auto-detection logic
2. **Add a `from_coordinates()` class method** to construct SNP objects from coordinate data
3. **Reuse existing `blastdbcmd` logic** already in `FlankingExtractor.extract_sequences()`

This reduces code duplication and maintains a single parsing interface.

## Proposed Changes

### 1. Parser Enhancement (`parser.py`)

```python
class PolymarkerParser:
    # Existing code...
    
    @classmethod
    def detect_format(cls, input_file: Path) -> str:
        """Detect if input is 'polymarker' or 'coordinates' format."""
        with open(input_file) as f:
            first_line = f.readline().strip()
            if ',' in first_line and '[' in first_line:
                return 'polymarker'
            elif len(first_line.split()) == 4:
                return 'coordinates'
        raise ParseError("Unknown input format")
    
    def parse_coordinates(self, reference: Path) -> List[SNP]:
        """Parse coordinate file and fetch sequences."""
        # 1. Parse lines into (chrom, pos, ref, alt)
        # 2. Prepare batch file for blastdbcmd
        # 3. Fetch sequences
        # 4. Insert [Ref/Alt] and create SNP objects
```

### 2. Main Pipeline Update (`main.py`)

```python
# Step 1: Detect format and parse
input_format = PolymarkerParser.detect_format(config.input_file)
parser = PolymarkerParser(config.input_file)

if input_format == 'coordinates':
    snps = parser.parse_coordinates(config.reference_file)
else:
    snps = parser.parse()
```

## Task List

- [ ] **Extend `PolymarkerParser`**
    - [ ] Add `detect_format()` static method
    - [ ] Add `parse_coordinates(reference_path)` method
    - [ ] Implement `blastdbcmd` batch sequence fetching (reuse `FlankingExtractor` logic)
    - [ ] Implement sequence reconstruction with `[Ref/Alt]` insertion

- [ ] **Update `main.py`**
    - [ ] Add format detection at Step 1
    - [ ] Branch to `parse_coordinates()` when needed
    - [ ] Pass reference path to parser for coordinate mode

- [ ] **Verification**
    - [ ] Test with `test_data/v2_examples/snp_pos_example.txt`
    - [ ] Compare with existing Polymarker input results
    - [ ] Ensure backward compatibility

## Benefits of This Approach
1. **Code Reuse**: Leverages existing `FlankingExtractor` and `blastdbcmd` integration
2. **Single Interface**: Users only need to know about `PolymarkerParser`
3. **Auto-Detection**: No need for explicit `--input-format` flag
4. **Minimal Changes**: Only modifies `parser.py` and `main.py` (Step 1 only)
