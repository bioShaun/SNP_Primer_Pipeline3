"""Configuration management for SNP primer pipeline."""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import yaml

from .exceptions import ConfigurationError


@dataclass
class PipelineConfig:
    """Pipeline configuration settings."""
    
    input_file: Path
    reference_file: Path
    output_dir: Path = Path("output")
    design_kasp: bool = True
    design_caps: bool = True
    max_price: int = 200
    max_tm: float = 63.0
    max_primer_size: int = 25
    pick_anyway: bool = False
    threads: int = 1
    flanking_size: int = 500
    max_hits: int = 6
    primer_product_size_range: Tuple[int, int] = (50, 250)
    log_level: str = "INFO"
    
    def __post_init__(self):
        """Validate configuration after initialization."""
        self.input_file = Path(self.input_file)
        self.reference_file = Path(self.reference_file)
        self.output_dir = Path(self.output_dir)
        
        if not self.input_file.exists():
            raise ConfigurationError(f"Input file not found: {self.input_file}")
        
        if self.max_tm <= 0 or self.max_tm > 100:
            raise ConfigurationError(f"Invalid max_tm: {self.max_tm}")
            
        if self.max_primer_size <= 0 or self.max_primer_size > 100:
            raise ConfigurationError(f"Invalid max_primer_size: {self.max_primer_size}")
            
        if self.threads <= 0:
            raise ConfigurationError(f"Invalid threads: {self.threads}")
    
    @classmethod
    def from_yaml(cls, yaml_file: Path) -> "PipelineConfig":
        """Load configuration from YAML file."""
        yaml_file = Path(yaml_file)
        if not yaml_file.exists():
            raise ConfigurationError(f"Config file not found: {yaml_file}")
        
        try:
            with open(yaml_file, 'r') as f:
                data = yaml.safe_load(f)
            
            # Convert paths to Path objects
            if 'input_file' in data:
                data['input_file'] = Path(data['input_file'])
            if 'reference_file' in data:
                data['reference_file'] = Path(data['reference_file'])
            if 'output_dir' in data:
                data['output_dir'] = Path(data['output_dir'])
                
            return cls(**data)
        except yaml.YAMLError as e:
            raise ConfigurationError(f"Invalid YAML config: {e}")
        except TypeError as e:
            raise ConfigurationError(f"Invalid config parameters: {e}")
    
    @classmethod
    def from_args(cls, args: dict) -> "PipelineConfig":
        """Create configuration from command-line arguments."""
        # Map command-line argument names to config field names
        arg_mapping = {
            'input': 'input_file',
            'reference': 'reference_file',
            'output': 'output_dir',
            'caps': 'design_caps',
            'kasp': 'design_kasp',
            'price': 'max_price',
            'max_tm': 'max_tm',
            'max_size': 'max_primer_size',
            'pick_anyway': 'pick_anyway',
            'threads': 'threads',
            'log_level': 'log_level',
        }
        
        config_args = {}
        for arg_name, config_name in arg_mapping.items():
            if arg_name in args and args[arg_name] is not None:
                config_args[config_name] = args[arg_name]
        
        # Convert boolean flags (0/1 to bool)
        if 'design_caps' in config_args:
            config_args['design_caps'] = bool(config_args['design_caps'])
        if 'design_kasp' in config_args:
            config_args['design_kasp'] = bool(config_args['design_kasp'])
        if 'pick_anyway' in config_args:
            config_args['pick_anyway'] = bool(config_args['pick_anyway'])
        
        return cls(**config_args)


@dataclass
class SoftwarePaths:
    """External software paths."""
    
    primer3_path: Path
    muscle_path: Path
    
    @classmethod
    def auto_detect(cls) -> "SoftwarePaths":
        """Auto-detect software paths based on platform."""
        base_path = Path(__file__).parent.parent.parent / "bin"
        
        if sys.platform.startswith('linux'):
            primer3_path = base_path / "primer3_core"
            # Use V2 muscle path to ensure compatibility
            muscle_path = Path(__file__).parent.parent.parent.parent / "SNP_Primer_Pipeline2" / "bin" / "muscle"
        elif sys.platform == "win32" or sys.platform == "cygwin":
            primer3_path = base_path / "primer3_core.exe"
            muscle_path = base_path / "muscle.exe"
        elif sys.platform == "darwin":  # macOS
            primer3_path = base_path / "primer3_core_darwin64"
            muscle_path = base_path / "muscle3.8.31_i86darwin64"
        else:
            raise ConfigurationError(f"Unsupported platform: {sys.platform}")
        
        # Check if files exist, fall back to system PATH
        if not primer3_path.exists():
            primer3_path = cls._find_in_path("primer3_core")
        if not muscle_path.exists():
            muscle_path = cls._find_in_path("muscle")
        
        return cls(primer3_path=primer3_path, muscle_path=muscle_path)
    
    @staticmethod
    def _find_in_path(executable: str) -> Path:
        """Find executable in system PATH."""
        for path in os.environ.get("PATH", "").split(os.pathsep):
            exe_path = Path(path) / executable
            if exe_path.exists() and exe_path.is_file():
                return exe_path
        
        raise ConfigurationError(f"Could not find {executable} in PATH or bundled binaries")