#!/usr/bin/env python3
"""Root compatibility wrapper for the canonical `Gene_Analysis.cli` entrypoint."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _load_cli_main():
    repo_root = Path(__file__).resolve().parent
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    from Gene_Analysis.cli import main as cli_main
    return cli_main


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compatibility wrapper for full genome analysis.",
    )
    parser.add_argument(
        "genome_file",
        nargs="?",
        default="Genome.txt",
        help="Path to 23andMe/raw genome text file.",
    )
    parser.add_argument(
        "--no-auto-update",
        action="store_true",
        help="Disable automatic evidence refresh when snapshot is stale.",
    )
    parser.add_argument(
        "--json-output",
        default=None,
        help="Optional explicit JSON output path.",
    )
    parser.add_argument(
        "--markdown-output",
        default=None,
        help="Optional explicit Markdown output path.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    cli_args = ["analyze", args.genome_file]
    if args.no_auto_update:
        cli_args.append("--no-auto-update")
    if args.json_output:
        cli_args.extend(["--json-output", args.json_output])
    if args.markdown_output:
        cli_args.extend(["--markdown-output", args.markdown_output])

    cli_main = _load_cli_main()
    return int(cli_main(cli_args))


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

