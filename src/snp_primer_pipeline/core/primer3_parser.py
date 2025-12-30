#!/usr/bin/env python3
"""
Primer3 interface module for SNP Primer Pipeline.

This module handles Primer3 input generation, execution, and output parsing.
"""

import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import re

from ..models import Primer, PrimerPair
from ..exceptions import PrimerDesignError


class Primer3Input:
    """Primer3 input generator."""
    
    def __init__(self, config_path: Optional[Path] = None):
        """
        Initialize Primer3 input generator.
        
        Args:
            config_path: Path to Primer3 configuration file
        """
        self.config_path = Path(config_path) if config_path else None
        self.settings: Dict[str, Any] = {}
        self._load_default_settings()
    
    def _load_default_settings(self) -> None:
        """Load default Primer3 settings."""
        # Set thermodynamic parameters path - try multiple locations
        config_dir = None
        possible_paths = [
            Path(__file__).parent.parent.parent / "bin" / "primer3_config",
            Path.cwd() / "bin" / "primer3_config",
            Path(__file__).parent.parent.parent.parent / "bin" / "primer3_config"
        ]
        
        for path in possible_paths:
            if path.exists():
                config_dir = path
                break
        
        self.settings = {
            "PRIMER_TASK": "generic",
            "PRIMER_PICK_LEFT_PRIMER": 1,
            "PRIMER_PICK_INTERNAL_OLIGO": 0,
            "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_OPT_SIZE": 20,
            "PRIMER_MIN_SIZE": 18,
            "PRIMER_MAX_SIZE": 25,
            "PRIMER_OPT_TM": 60.0,
            "PRIMER_MIN_TM": 57.0,
            "PRIMER_MAX_TM": 63.0,
            "PRIMER_MIN_GC": 20.0,
            "PRIMER_MAX_GC": 80.0,
            "PRIMER_MAX_POLY_X": 100,
            "PRIMER_INTERNAL_MAX_POLY_X": 100,
            "PRIMER_SALT_MONOVALENT": 50.0,
            "PRIMER_DNA_CONC": 50.0,
            "PRIMER_MAX_NS_ACCEPTED": 0,
            "PRIMER_MAX_SELF_ANY": 12,
            "PRIMER_MAX_SELF_END": 8,
            "PRIMER_PAIR_MAX_COMPL_ANY": 12,
            "PRIMER_PAIR_MAX_COMPL_END": 8,
            "PRIMER_PRODUCT_SIZE_RANGE": "50-75 75-100 100-125 125-150 150-175 175-200 200-225",
            "PRIMER_NUM_RETURN": 5,
            "PRIMER_FIRST_BASE_INDEX": 1,  # V2 compatibility: use 1-based indexing
            "PRIMER_LIBERAL_BASE": 1,  # V2 compatibility: allow ambiguous bases
        }
        
        # Add thermodynamic parameters path if config directory exists
        if config_dir and config_dir.exists():
            # Use relative path from current working directory
            try:
                relative_path = config_dir.relative_to(Path.cwd())
                self.settings["PRIMER_THERMODYNAMIC_PARAMETERS_PATH"] = str(relative_path) + "/"
            except ValueError:
                # If relative path fails, use absolute path
                self.settings["PRIMER_THERMODYNAMIC_PARAMETERS_PATH"] = str(config_dir) + "/"
    
    def load_config(self, config_path: Path) -> "Primer3Input":
        """
        Load settings from configuration file.
        
        Args:
            config_path: Path to configuration file
            
        Returns:
            Self for method chaining
        """
        try:
            with open(config_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    if '=' in line:
                        key, value = line.split('=', 1)
                        key = key.strip()
                        value = value.strip()
                        
                        # Try to convert to appropriate type
                        if value.replace('.', '').replace('-', '').isdigit():
                            if '.' in value:
                                value = float(value)
                            else:
                                value = int(value)
                        
                        self.settings[key] = value
        except IOError as e:
            raise PrimerDesignError(f"Failed to load config file: {e}") from e
        
        return self
    
    def set_template(self, sequence: str) -> "Primer3Input":
        """Set template sequence."""
        self.settings["SEQUENCE_TEMPLATE"] = sequence
        return self
    
    def set_product_size_range(self, ranges: List[Tuple[int, int]]) -> "Primer3Input":
        """Set product size ranges."""
        range_str = " ".join(f"{start}-{end}" for start, end in ranges)
        self.settings["PRIMER_PRODUCT_SIZE_RANGE"] = range_str
        return self
    
    def set_force_left_end(self, position: int) -> "Primer3Input":
        """Force left primer to end at specific position."""
        self.settings["SEQUENCE_FORCE_LEFT_END"] = position
        return self
    
    def set_force_right_end(self, position: int) -> "Primer3Input":
        """Force right primer to end at specific position."""
        self.settings["SEQUENCE_FORCE_RIGHT_END"] = position
        return self
    
    def set_target(self, start: int, length: int) -> "Primer3Input":
        """Set target region."""
        self.settings["TARGET"] = f"{start},{length}"
        return self
    
    def set_excluded_region(self, start: int, length: int) -> "Primer3Input":
        """Set excluded region."""
        if "PRIMER_EXCLUDED_REGION" not in self.settings:
            self.settings["PRIMER_EXCLUDED_REGION"] = []
        if not isinstance(self.settings["PRIMER_EXCLUDED_REGION"], list):
            self.settings["PRIMER_EXCLUDED_REGION"] = [self.settings["PRIMER_EXCLUDED_REGION"]]
        self.settings["PRIMER_EXCLUDED_REGION"].append(f"{start},{length}")
        return self
    
    def set_setting(self, key: str, value: Any) -> "Primer3Input":
        """Set arbitrary Primer3 setting."""
        self.settings[key] = value
        return self
    
    def generate(self, sequence_id: str) -> str:
        """
        Generate Primer3 input string.
        
        Args:
            sequence_id: Sequence identifier
            
        Returns:
            Primer3 input string
        """
        lines = [f"SEQUENCE_ID={sequence_id}"]
        
        for key, value in self.settings.items():
            if isinstance(value, list):
                for v in value:
                    lines.append(f"{key}={v}")
            else:
                lines.append(f"{key}={value}")
        
        lines.append("=")  # End marker
        return "\n".join(lines)


class Primer3Runner:
    """Primer3 execution wrapper."""
    
    def __init__(self, primer3_path: Optional[Path] = None, settings_file: Optional[Path] = None):
        """
        Initialize Primer3 runner.
        
        Args:
            primer3_path: Path to primer3_core executable
            settings_file: Path to global settings file
        """
        self.primer3_path = primer3_path or Path("primer3_core")
        self.settings_file = settings_file
    
    def run(self, input_file: Path, output_file: Path) -> Path:
        """
        Execute Primer3.
        
        Args:
            input_file: Input file path
            output_file: Output file path
            
        Returns:
            Path to output file
            
        Raises:
            PrimerDesignError: If Primer3 execution fails
        """
        cmd = [str(self.primer3_path)]
        
        if self.settings_file:
            cmd.extend(["-p3_settings_file", str(self.settings_file)])
        
        try:
            with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
                result = subprocess.run(
                    cmd,
                    stdin=infile,
                    stdout=outfile,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True
                )
            
            return Path(output_file)
            
        except subprocess.CalledProcessError as e:
            raise PrimerDesignError(f"Primer3 execution failed: {e.stderr}") from e
        except FileNotFoundError as e:
            raise PrimerDesignError(f"Primer3 not found: {e}") from e
    
    def run_string(self, input_string: str) -> str:
        """
        Execute Primer3 with string input.
        
        Args:
            input_string: Primer3 input string
            
        Returns:
            Primer3 output string
        """
        cmd = [str(self.primer3_path)]
        
        if self.settings_file:
            cmd.extend(["-p3_settings_file", str(self.settings_file)])
        
        try:
            result = subprocess.run(
                cmd,
                input=input_string,
                capture_output=True,
                text=True,
                check=True
            )
            
            return result.stdout
            
        except subprocess.CalledProcessError as e:
            raise PrimerDesignError(f"Primer3 execution failed: {e.stderr}") from e
        except FileNotFoundError as e:
            raise PrimerDesignError(f"Primer3 not found: {e}") from e


class Primer3OutputParser:
    """Primer3 output parser."""
    
    def __init__(self):
        """Initialize parser."""
        pass
    
    def parse(self, output_file: Path, num_pairs: int = 5) -> Dict[str, List[PrimerPair]]:
        """
        Parse Primer3 output file.
        
        Args:
            output_file: Primer3 output file
            num_pairs: Maximum number of primer pairs to parse
            
        Returns:
            Dictionary mapping sequence IDs to lists of PrimerPair objects
        """
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            return self.parse_string(content, num_pairs)
        except IOError as e:
            raise PrimerDesignError(f"Failed to read output file: {e}") from e
    
    def parse_string(self, output_string: str, num_pairs: int = 5) -> Dict[str, List[PrimerPair]]:
        """
        Parse Primer3 output string.
        
        Args:
            output_string: Primer3 output string
            num_pairs: Maximum number of primer pairs to parse
            
        Returns:
            Dictionary mapping sequence IDs to lists of PrimerPair objects
        """
        results = {}
        
        # Split output into records
        records = output_string.strip().split("=\n")
        
        for record in records:
            if not record.strip():
                continue
            
            # Parse record
            data = self._parse_record(record)
            
            if "SEQUENCE_ID" not in data:
                continue
            
            sequence_id = data["SEQUENCE_ID"]
            primer_pairs = self._extract_primer_pairs(data, num_pairs)
            
            if primer_pairs:
                results[sequence_id] = primer_pairs
        
        return results
    
    def _parse_record(self, record: str) -> Dict[str, str]:
        """Parse a single Primer3 output record."""
        data = {}
        
        for line in record.split('\n'):
            line = line.strip()
            if not line or '=' not in line:
                continue
            
            key, value = line.split('=', 1)
            data[key] = value
        
        return data
    
    def _extract_primer_pairs(self, data: Dict[str, str], num_pairs: int) -> List[PrimerPair]:
        """Extract primer pairs from parsed data."""
        primer_pairs = []
        
        # Check if any primers were found
        primers_returned = int(data.get("PRIMER_PAIR_NUM_RETURNED", 0))
        
        if primers_returned == 0:
            return primer_pairs
        
        for i in range(min(primers_returned, num_pairs)):
            try:
                # Extract left primer
                left_seq = data.get(f"PRIMER_LEFT_{i}_SEQUENCE", "")
                left_pos = data.get(f"PRIMER_LEFT_{i}", "")
                left_tm = float(data.get(f"PRIMER_LEFT_{i}_TM", 0))
                left_gc = float(data.get(f"PRIMER_LEFT_{i}_GC_PERCENT", 0))
                left_self_any = float(data.get(f"PRIMER_LEFT_{i}_SELF_ANY_TH", 0))
                left_self_end = float(data.get(f"PRIMER_LEFT_{i}_SELF_END_TH", 0))
                left_hairpin = float(data.get(f"PRIMER_LEFT_{i}_HAIRPIN_TH", 0))
                left_end_stability = float(data.get(f"PRIMER_LEFT_{i}_END_STABILITY", 0))
                
                # Parse position (format: start,length)
                # Note: Assuming 1-based indexing from Primer3 (PRIMER_FIRST_BASE_INDEX=1)
                if left_pos and ',' in left_pos:
                    left_p3_start, left_length = map(int, left_pos.split(','))
                    left_start = left_p3_start - 1  # Convert to 0-based
                    left_end = left_start + left_length - 1
                else:
                    left_start = left_end = left_length = 0
                
                left_primer = Primer(
                    name=f"LEFT_{i}",
                    sequence=left_seq,
                    start=left_start,
                    end=left_end,
                    length=left_length,
                    tm=left_tm,
                    gc_percent=left_gc,
                    self_any=left_self_any,
                    self_end=left_self_end,
                    hairpin=left_hairpin,
                    end_stability=left_end_stability,
                    direction="LEFT"
                )
                
                # Extract right primer
                right_seq = data.get(f"PRIMER_RIGHT_{i}_SEQUENCE", "")
                right_pos = data.get(f"PRIMER_RIGHT_{i}", "")
                right_tm = float(data.get(f"PRIMER_RIGHT_{i}_TM", 0))
                right_gc = float(data.get(f"PRIMER_RIGHT_{i}_GC_PERCENT", 0))
                right_self_any = float(data.get(f"PRIMER_RIGHT_{i}_SELF_ANY_TH", 0))
                right_self_end = float(data.get(f"PRIMER_RIGHT_{i}_SELF_END_TH", 0))
                right_hairpin = float(data.get(f"PRIMER_RIGHT_{i}_HAIRPIN_TH", 0))
                right_end_stability = float(data.get(f"PRIMER_RIGHT_{i}_END_STABILITY", 0))
                
                # Parse position (format: start,length)
                # Note: Assuming 1-based indexing from Primer3
                # For RIGHT primer, pos is the 5' end (highest index)
                if right_pos and ',' in right_pos:
                    right_p3_start, right_length = map(int, right_pos.split(','))
                    right_end = right_p3_start - 1  # Convert to 0-based 5' end
                    right_start = right_end - right_length + 1  # Calculate 3' end
                else:
                    right_start = right_end = right_length = 0
                
                right_primer = Primer(
                    name=f"RIGHT_{i}",
                    sequence=right_seq,
                    start=right_start,
                    end=right_end,
                    length=right_length,
                    tm=right_tm,
                    gc_percent=right_gc,
                    self_any=right_self_any,
                    self_end=right_self_end,
                    hairpin=right_hairpin,
                    end_stability=right_end_stability,
                    direction="RIGHT"
                )
                
                # Extract pair information
                product_size = int(data.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", 0))
                penalty = float(data.get(f"PRIMER_PAIR_{i}_PENALTY", 0))
                compl_any = float(data.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH", 0))
                compl_end = float(data.get(f"PRIMER_PAIR_{i}_COMPL_END_TH", 0))
                
                primer_pair = PrimerPair(
                    left=left_primer,
                    right=right_primer,
                    product_size=product_size,
                    penalty=penalty,
                    compl_any=compl_any,
                    compl_end=compl_end
                )
                
                primer_pairs.append(primer_pair)
                
            except (ValueError, KeyError) as e:
                # Skip malformed primer pairs
                continue
        
        return primer_pairs