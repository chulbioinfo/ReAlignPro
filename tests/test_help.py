"""Smoke tests for the ReAlignPro CLI.

These tests intentionally avoid running the full pipeline. They only verify that:
  - the package can be imported
  - the CLI entry point responds to --help/--version for each subcommand

Rationale:
  - Fast and stable in CI (no large data, no heavy external tools).
  - Avoids brittle, long-running integration tests in PR checks.
"""

from __future__ import annotations

import subprocess
import sys


def _run(args: list[str]) -> None:
    # Use `python -m realignpro` to avoid PATH ambiguity between conda/bin and ~/.local/bin.
    cmd = [sys.executable, "-m", "realignpro"] + args
    subprocess.run(cmd, check=True, timeout=30)


def test_module_import() -> None:
    subprocess.run([sys.executable, "-c", "import realignpro; import realignpro.cli"], check=True, timeout=30)


def test_top_level_help_and_version() -> None:
    _run(["--help"])
    _run(["--version"])


def test_subcommand_help() -> None:
    _run(["fa2maf", "--help"])
    _run(["maf2bed", "--help"])
    _run(["tsv2fig", "--help"])
