"""Configuration management for SNP primer pipeline."""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path

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
    primer_product_size_range: tuple[int, int] = (50, 250)
    log_level: str = "INFO"
    show_variant_sites: bool = False  # Hide variant sites by default for cleaner output
    run_specificity: bool = True  # Run off-target specificity assessment
    specificity_warning_threshold: int = 50  # Warning hit count threshold
    specificity_target_window: int = 2000  # Target exclusion window (bp)

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
    def from_yaml(cls, yaml_file: Path) -> PipelineConfig:
        """Load configuration from YAML file."""
        yaml_file = Path(yaml_file)
        if not yaml_file.exists():
            raise ConfigurationError(f"Config file not found: {yaml_file}")

        try:
            with open(yaml_file) as f:
                data = yaml.safe_load(f)

            # Convert paths to Path objects
            if "input_file" in data:
                data["input_file"] = Path(data["input_file"])
            if "reference_file" in data:
                data["reference_file"] = Path(data["reference_file"])
            if "output_dir" in data:
                data["output_dir"] = Path(data["output_dir"])

            return cls(**data)
        except yaml.YAMLError as e:
            raise ConfigurationError(f"Invalid YAML config: {e}") from e
        except TypeError as e:
            raise ConfigurationError(f"Invalid config parameters: {e}") from e

    @classmethod
    def from_args(cls, args: dict) -> PipelineConfig:
        """Create configuration from command-line arguments."""
        # Map command-line argument names to config field names
        arg_mapping = {
            "input": "input_file",
            "reference": "reference_file",
            "output": "output_dir",
            "caps": "design_caps",
            "kasp": "design_kasp",
            "price": "max_price",
            "max_tm": "max_tm",
            "max_size": "max_primer_size",
            "pick_anyway": "pick_anyway",
            "threads": "threads",
            "log_level": "log_level",
        }

        config_args = {}
        for arg_name, config_name in arg_mapping.items():
            if arg_name in args and args[arg_name] is not None:
                config_args[config_name] = args[arg_name]

        # Convert boolean flags (0/1 to bool)
        if "design_caps" in config_args:
            config_args["design_caps"] = bool(config_args["design_caps"])
        if "design_kasp" in config_args:
            config_args["design_kasp"] = bool(config_args["design_kasp"])
        if "pick_anyway" in config_args:
            config_args["pick_anyway"] = bool(config_args["pick_anyway"])

        return cls(**config_args)


@dataclass
class SoftwarePaths:
    """External software paths."""

    primer3_path: Path
    muscle_path: Path

    @classmethod
    def auto_detect(cls) -> SoftwarePaths:
        """Auto-detect software paths based on platform.

        Search order:
        1. System PATH
        2. Package's bin directory (inside the installed package)
        3. Current working directory bin (development mode)
        """
        # 1. Try to find in system PATH first
        primer3_path = cls._find_in_path("primer3_core", raise_error=False)
        muscle_path = cls._find_in_path("muscle", raise_error=False)

        # 2. If not found, look in bundled binaries or local bin
        if primer3_path is None or muscle_path is None:
            # Possible bin directory locations
            search_paths = [
                # Package's bin directory (inside the installed package)
                Path(__file__).parent / "bin",
                # Current working directory bin (useful for development)
                Path.cwd() / "bin",
            ]

            # Platform-specific binary names for bundled files
            if sys.platform.startswith("linux"):
                primer3_bundle_name = "primer3_core"
                muscle_bundle_name = "muscle"
            elif sys.platform == "win32" or sys.platform == "cygwin":
                primer3_bundle_name = "primer3_core.exe"
                muscle_bundle_name = "muscle.exe"
            elif sys.platform == "darwin":  # macOS
                primer3_bundle_name = "primer3_core_darwin64"
                muscle_bundle_name = "muscle3.8.31_i86darwin64"
            else:
                # Unknown platform, default to linux names if we can't determine
                primer3_bundle_name = "primer3_core"
                muscle_bundle_name = "muscle"

            if primer3_path is None:
                primer3_path = cls._find_in_search_paths(search_paths, primer3_bundle_name)

            if muscle_path is None:
                muscle_path = cls._find_in_search_paths(search_paths, muscle_bundle_name)

        # 3. Final check
        if primer3_path is None:
            raise ConfigurationError(
                "Could not find primer3_core in system PATH or bundled binaries.\n"
                "Please install it using: conda install -c bioconda primer3  (or brew install primer3 on macOS)"
            )
        if muscle_path is None:
            raise ConfigurationError(
                "Could not find muscle in system PATH or bundled binaries.\n"
                "Please install it using: conda install -c bioconda muscle"
            )

        return cls(primer3_path=primer3_path, muscle_path=muscle_path)

    @classmethod
    def _find_in_search_paths(cls, search_paths: list[Path], filename: str) -> Path | None:
        """Find a file in a list of search directories."""
        for base_path in search_paths:
            file_path = base_path / filename
            if file_path.exists() and file_path.is_file():
                return file_path
        return None

    @staticmethod
    def _find_in_path(executable: str, raise_error: bool = True) -> Path | None:
        """Find executable in system PATH."""
        for path in os.environ.get("PATH", "").split(os.pathsep):
            try:
                exe_path = Path(path) / executable
                # Check if it exists and is executable
                if exe_path.exists() and exe_path.is_file() and os.access(exe_path, os.X_OK):
                    return exe_path
            except (OSError, ValueError):
                continue

        if raise_error:
            raise ConfigurationError(
                f"Could not find {executable} in system PATH.\n"
                f"Please install it using: conda install -c bioconda primer3  (or brew install primer3 on macOS)"
            )
        return None
