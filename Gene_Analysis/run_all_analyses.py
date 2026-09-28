#!/usr/bin/env python3
"""Compatibility wrapper for the canonical pipeline, including former option names."""
from __future__ import annotations
from pathlib import Path
import sys


def main(argv: list[str] | None = None) -> int:
    if not __package__:
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from Gene_Analysis.cli import main as cli_main
    arguments = list(sys.argv[1:] if argv is None else argv)
    renamed = {'--json': '--json-output', '--markdown': '--markdown-output'}
    arguments = [renamed.get(value, value) for value in arguments]
    if '--genome' in arguments:
        index = arguments.index('--genome')
        if index + 1 >= len(arguments):
            print('Error: --genome requires a file path.', file=sys.stderr)
            return 2
        genome = arguments[index + 1]
        del arguments[index:index + 2]
        arguments.insert(0, genome)
    if arguments and arguments[0] in {'analyze', 'demo', 'update-evidence', 'database-status'}:
        return cli_main(arguments)
    return cli_main(['analyze', *arguments])


if __name__ == '__main__':
    raise SystemExit(main())
