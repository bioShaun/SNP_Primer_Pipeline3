"""Unit tests for PipelineConfig and SoftwarePaths (P2 config items)."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from snp_primer_pipeline.config import PipelineConfig, SoftwarePaths
from snp_primer_pipeline.exceptions import ConfigurationError

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_input_and_ref(tmp_path: Path) -> dict:
    """Create minimal input/reference files and return default kwargs."""
    inp = tmp_path / "input.txt"
    ref = tmp_path / "ref.fa"
    inp.write_text("dummy")
    ref.write_text("dummy")
    return {"input_file": inp, "reference_file": ref}


# ---------------------------------------------------------------------------
# PipelineConfig.__post_init__
# ---------------------------------------------------------------------------


class TestPostInit:
    def test_valid_config(self, tmp_path: Path) -> None:
        cfg = PipelineConfig(**_make_input_and_ref(tmp_path))
        assert cfg.input_file.exists()

    def test_input_file_missing(self, tmp_path: Path) -> None:
        kw = _make_input_and_ref(tmp_path)
        kw["input_file"] = tmp_path / "no_such_file.txt"
        with pytest.raises(ConfigurationError, match="Input file not found"):
            PipelineConfig(**kw)

    @pytest.mark.parametrize("bad_tm", [0, -1, 100.1, 200])
    def test_max_tm_invalid(self, tmp_path: Path, bad_tm: float) -> None:
        kw = _make_input_and_ref(tmp_path)
        kw["max_tm"] = bad_tm
        with pytest.raises(ConfigurationError, match="Invalid max_tm"):
            PipelineConfig(**kw)

    @pytest.mark.parametrize("bad_size", [0, -5, 101])
    def test_max_primer_size_invalid(self, tmp_path: Path, bad_size: int) -> None:
        kw = _make_input_and_ref(tmp_path)
        kw["max_primer_size"] = bad_size
        with pytest.raises(ConfigurationError, match="Invalid max_primer_size"):
            PipelineConfig(**kw)

    @pytest.mark.parametrize("bad_threads", [0, -1])
    def test_threads_invalid(self, tmp_path: Path, bad_threads: int) -> None:
        kw = _make_input_and_ref(tmp_path)
        kw["threads"] = bad_threads
        with pytest.raises(ConfigurationError, match="Invalid threads"):
            PipelineConfig(**kw)


# ---------------------------------------------------------------------------
# PipelineConfig.from_yaml
# ---------------------------------------------------------------------------


class TestFromYaml:
    def test_valid_yaml(self, tmp_path: Path) -> None:
        inp = tmp_path / "input.txt"
        ref = tmp_path / "ref.fa"
        inp.write_text("dummy")
        ref.write_text("dummy")

        yaml_file = tmp_path / "config.yaml"
        yaml_file.write_text(
            f"input_file: {inp}\nreference_file: {ref}\nmax_tm: 65.0\nthreads: 2\n"
        )
        cfg = PipelineConfig.from_yaml(yaml_file)
        assert cfg.max_tm == 65.0
        assert cfg.threads == 2

    def test_yaml_file_missing(self, tmp_path: Path) -> None:
        with pytest.raises(ConfigurationError, match="Config file not found"):
            PipelineConfig.from_yaml(tmp_path / "nonexistent.yaml")

    def test_invalid_yaml_content(self, tmp_path: Path) -> None:
        yaml_file = tmp_path / "bad.yaml"
        yaml_file.write_text(":\n  - :\n  bad: [unterminated")
        with pytest.raises(ConfigurationError, match="Invalid YAML config"):
            PipelineConfig.from_yaml(yaml_file)

    def test_invalid_parameters(self, tmp_path: Path) -> None:
        inp = tmp_path / "input.txt"
        ref = tmp_path / "ref.fa"
        inp.write_text("dummy")
        ref.write_text("dummy")

        yaml_file = tmp_path / "bad_params.yaml"
        yaml_file.write_text(f"input_file: {inp}\nreference_file: {ref}\nnot_a_real_field: 42\n")
        with pytest.raises(ConfigurationError, match="Invalid config parameters"):
            PipelineConfig.from_yaml(yaml_file)


# ---------------------------------------------------------------------------
# SoftwarePaths.auto_detect
# ---------------------------------------------------------------------------


class TestAutoDetect:
    def test_both_found_in_path(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()
        primer3 = bin_dir / "primer3_core"
        muscle = bin_dir / "muscle"
        for f in (primer3, muscle):
            f.write_text("#!/bin/sh\n")
            f.chmod(0o755)

        monkeypatch.setenv("PATH", str(bin_dir))
        paths = SoftwarePaths.auto_detect()
        assert paths.primer3_path == primer3
        assert paths.muscle_path == muscle

    def test_found_in_bundled_bin(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        # Ensure nothing in PATH
        monkeypatch.setenv("PATH", "")

        # Determine platform-specific names used by auto_detect
        if sys.platform.startswith("linux"):
            p3_name, mu_name = "primer3_core", "muscle"
        elif sys.platform in ("win32", "cygwin"):
            p3_name, mu_name = "primer3_core.exe", "muscle.exe"
        elif sys.platform == "darwin":
            p3_name, mu_name = "primer3_core_darwin64", "muscle3.8.31_i86darwin64"
        else:
            p3_name, mu_name = "primer3_core", "muscle"

        # Create fake bundled bin dir and monkey-patch __file__ location
        pkg_dir = tmp_path / "pkg"
        pkg_dir.mkdir()
        bin_dir = pkg_dir / "bin"
        bin_dir.mkdir()
        for name in (p3_name, mu_name):
            f = bin_dir / name
            f.write_text("#!/bin/sh\n")
            f.chmod(0o755)

        # Patch Path(__file__).parent to point to pkg_dir
        import snp_primer_pipeline.config as config_mod

        monkeypatch.setattr(config_mod, "__file__", str(pkg_dir / "config.py"))

        paths = SoftwarePaths.auto_detect()
        assert paths.primer3_path == bin_dir / p3_name
        assert paths.muscle_path == bin_dir / mu_name

    def test_neither_found(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("PATH", "")
        # Point bundled bin to an empty dir
        import snp_primer_pipeline.config as config_mod

        empty = tmp_path / "empty_pkg"
        empty.mkdir()
        monkeypatch.setattr(config_mod, "__file__", str(empty / "config.py"))
        monkeypatch.chdir(tmp_path)

        with pytest.raises(ConfigurationError, match="Could not find primer3_core"):
            SoftwarePaths.auto_detect()


# ---------------------------------------------------------------------------
# SoftwarePaths._find_in_path
# ---------------------------------------------------------------------------


class TestFindInPath:
    def test_found(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        exe = tmp_path / "my_tool"
        exe.write_text("#!/bin/sh\n")
        exe.chmod(0o755)
        monkeypatch.setenv("PATH", str(tmp_path))
        result = SoftwarePaths._find_in_path("my_tool")
        assert result == exe

    def test_not_found_no_raise(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("PATH", "")
        result = SoftwarePaths._find_in_path("nonexistent_tool", raise_error=False)
        assert result is None

    def test_not_found_raise(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("PATH", "")
        with pytest.raises(ConfigurationError, match="Could not find nonexistent_tool"):
            SoftwarePaths._find_in_path("nonexistent_tool", raise_error=True)
