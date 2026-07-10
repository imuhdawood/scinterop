"""Command-line interface for scinterop.

Usage::

    scinterop detect <path>
    scinterop convert <input> <output> [options]
    scinterop run <script> --executor {r,python} [options]
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from . import __version__, read, convert, run_r, run_python, detect, errors

logger = logging.getLogger("scinterop.cli")


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the CLI.

    Subcommands:

    - ``detect`` — Detect the format of a single-cell data file.
    - ``convert`` — Convert between single-cell formats.
    - ``run`` — Execute an external R or Python script.

    Returns:
        A configured ``ArgumentParser``.
    """
    parser = argparse.ArgumentParser(
        prog="scinterop",
        description="Single-cell data interoperability tool",
    )
    parser.add_argument(
        "--version", action="version", version=f"scinterop {__version__}"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable debug logging"
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # detect
    detect_parser = subparsers.add_parser(
        "detect", help="Detect the format of a file or directory"
    )
    detect_parser.add_argument("path", type=str, help="Path to file or directory")

    # convert
    convert_parser = subparsers.add_parser(
        "convert", help="Convert between single-cell formats"
    )
    convert_parser.add_argument("input", type=str, help="Input file or directory")
    convert_parser.add_argument("output", type=str, help="Output file or directory")
    convert_parser.add_argument(
        "--r-exe", type=str, default=None,
        help="Path to Rscript executable (default: SCINTEROP_R_EXE env or 'Rscript')",
    )
    convert_parser.add_argument(
        "--python-exe", type=str, default=None,
        help="Path to Python executable (default: SCINTEROP_PYTHON_EXE env or 'python')",
    )
    convert_parser.add_argument(
        "--seurat", action="store_true",
        help="When converting to RDS, create a Seurat object instead of a plain R list",
    )
    convert_parser.add_argument(
        "--debug", action="store_true",
        help="Keep temporary files for debugging",
    )

    # run
    run_parser = subparsers.add_parser(
        "run", help="Run an external R or Python script"
    )
    run_parser.add_argument("script", type=str, help="Path to script file")
    run_parser.add_argument(
        "--executor", "-e", type=str, required=True,
        choices=["r", "python"],
        help="Script language/executor",
    )
    run_parser.add_argument(
        "--exe", type=str, default=None,
        help="Path to Rscript or Python executable",
    )
    run_parser.add_argument(
        "--conda-env", type=str, default=None,
        help="Conda environment name (Python only)",
    )
    run_parser.add_argument(
        "--args", type=str, nargs=argparse.REMAINDER,
        help="Additional arguments passed to the script",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    """CLI entry point.

    Args:
        argv: Command-line arguments (defaults to ``sys.argv[1:]``).

    Returns:
        Exit code (0 for success, 1 for errors).
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(levelname)s | %(name)s | %(message)s",
        stream=sys.stderr,
    )

    try:
        if args.command == "detect":
            _cmd_detect(args)
        elif args.command == "convert":
            _cmd_convert(args)
        elif args.command == "run":
            _cmd_run(args)
        else:
            parser.print_help()
            return 1
    except errors.ScinteropError as e:
        logger.error(str(e))
        return 1
    except Exception as e:
        logger.exception("Unexpected error: %s", e)
        return 1

    return 0


def _cmd_detect(args) -> None:
    result = detect(args.path)
    print(f"Format:      {result.fmt}")
    print(f"Path:        {result.path}")
    if result.details:
        print("Details:")
        for key, val in result.details.items():
            print(f"  {key}: {val}")


def _cmd_convert(args) -> None:
    if args.debug:
        import os as _os
        _os.environ["SCINTEROP_DEBUG"] = "1"

    result_path = convert(
        args.input,
        args.output,
        r_exe=args.r_exe,
        python_exe=args.python_exe,
        seurat=args.seurat,
    )
    print(f"Output: {result_path}", file=sys.stderr)


def _cmd_run(args) -> None:
    script = Path(args.script)
    if not script.exists():
        logger.error("Script not found: %s", script)
        sys.exit(1)

    if args.executor == "r":
        result = run_r(script, r_exe=args.exe)
    elif args.executor == "python":
        result = run_python(
            script, python_exe=args.exe, conda_env=args.conda_env,
        )

    if result.stdout:
        sys.stdout.write(result.stdout)
    if result.stderr:
        sys.stderr.write(result.stderr)


if __name__ == "__main__":
    sys.exit(main())
