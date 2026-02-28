from __future__ import annotations

import sys
from typing import List, Optional

from . import __version__

# NOTE:
# Subcommands are imported lazily to avoid importing heavy dependencies
# (e.g., matplotlib) when users only request help/version.

HELP_TEXT = f"""ReAlignPro (version {__version__})

Usage:
  realignpro <command> [options]

Commands:
  fa2maf    Find orthologous sequences, align with MUSCLE, and export MAF/TSV outputs.
  maf2bed   Extract target-shared, others-absent positions from MAF and write BED3 intervals.
  tsv2fig   Plot sequence tables from conservation-matrix TSV (PDF outputs).

Run:
  realignpro <command> --help
for command-specific options.
"""


def main(argv: Optional[List[str]] = None) -> int:
    args = sys.argv[1:] if argv is None else list(argv)

    if not args or args[0] in ("-h", "--help", "help"):
        sys.stdout.write(HELP_TEXT)
        return 0

    if args[0] in ("-v", "--version", "version"):
        sys.stdout.write(__version__ + "\n")
        return 0

    cmd = args[0]
    rest = args[1:]

    if cmd == "fa2maf":
        from . import fa2maf
        return int(fa2maf.main(rest))

    if cmd == "maf2bed":
        from . import maf2bed
        return int(maf2bed.main(rest))

    if cmd == "tsv2fig":
        from . import tsv2fig
        return int(tsv2fig.main(rest))

    sys.stderr.write(f"[ERROR] Unknown command: {cmd}\n\n")
    sys.stdout.write(HELP_TEXT)
    return 2
