#!/usr/bin/env python3
"""Compatibility entry point; accepts the same options as gene-analysis analyze."""
from __future__ import annotations
import sys
from Gene_Analysis.cli import main as cli_main


def main(argv: list[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    if arguments and arguments[0] in {'analyze', 'demo', 'update-evidence', 'database-status'}:
        return cli_main(arguments)
    return cli_main(['analyze', *arguments])


if __name__ == '__main__':
    raise SystemExit(main())
