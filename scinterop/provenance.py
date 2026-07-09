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
