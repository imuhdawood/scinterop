"""Provenance and metadata logging for conversions.

Records input/output paths, formats, timestamps, success status,
and runtime into a JSON-per-line log file or in-memory dict.
"""

from __future__ import annotations

import json
import logging
import os
import platform
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .errors import ProvenanceError

logger = logging.getLogger(__name__)


def write_conversion_record(
    *,
    input_path: str,
    input_format: str,
    output_path: str,
    output_format: str,
    success: bool,
    runtime_s: float,
    error: str | None = None,
    metadata: dict[str, Any] | None = None,
    log_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Write a single conversion provenance record.

    The record is returned as a dict and optionally appended to a
    JSONL log file under ``log_dir``.

    Args:
        input_path: Source file path.
        input_format: Source format string.
        output_path: Output file path.
        output_format: Target format string.
        success: Whether the conversion succeeded.
        runtime_s: Wall-clock runtime in seconds.
        error: Error message if conversion failed.
        metadata: Extra information to attach to the record.
        log_dir: If set, append the record to ``provenance.jsonl``
            in this directory.

    Returns:
        The provenance record dictionary.

    Raises:
        ProvenanceError: If writing the log file fails.
    """
    record = _build_record(
        input_path=input_path,
        input_format=input_format,
        output_path=output_path,
        output_format=output_format,
        success=success,
        runtime_s=runtime_s,
        error=error,
        metadata=metadata,
    )

    if log_dir is not None:
        log_path = Path(log_dir) / "provenance.jsonl"
        try:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            _append_jsonl(record, log_path)
            logger.info("Appended provenance record to %s", log_path)
        except OSError as e:
            raise ProvenanceError(
                f"Failed to write provenance log to '{log_path}': {e}"
            ) from e

    return record


def _build_record(
    *,
    input_path: str,
    input_format: str,
    output_path: str,
    output_format: str,
    success: bool,
    runtime_s: float,
    error: str | None = None,
    metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    from . import __version__

    record = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "package": "scinterop",
        "version": __version__,
        "input_path": input_path,
        "input_format": input_format,
        "output_path": output_path,
        "output_format": output_format,
        "success": success,
        "runtime_seconds": round(runtime_s, 3),
        "system": {
            "platform": platform.system(),
            "release": platform.release(),
            "python": platform.python_version(),
            "hostname": platform.node(),
        },
    }

    if error is not None:
        record["error"] = error
    if metadata is not None:
        record["metadata"] = metadata

    return record


def _append_jsonl(record: dict, path: Path) -> None:
    with open(path, "a") as f:
        f.write(json.dumps(record, default=str) + "\n")


def read_provenance_log(path: str | Path) -> list[dict[str, Any]]:
    """Read a provenance JSONL log file.

    Args:
        path: Path to the JSONL log file.

    Returns:
        List of provenance record dicts.

    Raises:
        ProvenanceError: If the file is missing, unreadable, or
            contains invalid JSON.
    """
    path = Path(path)
    if not path.exists():
        raise ProvenanceError(f"Provenance log not found: {path}")

    records = []
    try:
        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    records.append(json.loads(line))
    except json.JSONDecodeError as e:
        raise ProvenanceError(
            f"Failed to parse provenance log '{path}': {e}"
        ) from e
    except OSError as e:
        raise ProvenanceError(
            f"Failed to read provenance log '{path}': {e}"
        ) from e

    return records
