from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .errors import DetectionError

logger = logging.getLogger(__name__)

EXTENSION_MAP: dict[str, str] = {
    ".h5ad": "h5ad",
    ".rds": "rds",
    ".mtx": "mtx",
}

COMPRESSED_EXTENSIONS: dict[str, str] = {
    ".mtx.gz": "mtx",
}

MTX_TRIO = {"matrix.mtx", "matrix.mtx.gz", "barcodes.tsv", "barcodes.tsv.gz",
            "features.tsv", "features.tsv.gz"}


@dataclass
class Detection:
    fmt: str
    path: Path
    details: dict[str, Any]


def detect(path: str | Path) -> Detection:
    p = Path(path)

    details: dict[str, Any] = {}

    if p.is_dir():
        return _detect_directory(p, details)

    if p.exists():
        return _detect_file(p, details)

    if p.suffix:
        result = _detect_by_extension(p, details)
        if result is not None:
            return result

    raise DetectionError(
        f"Cannot determine format from path: {p}. "
        f"Supported extensions: {', '.join(EXTENSION_MAP)}"
    )


def _detect_directory(p: Path, details: dict) -> Detection:
    files = {f.name for f in p.iterdir() if f.is_file()}
    present = MTX_TRIO & files
    n_present = len(present)
    if n_present >= 2:
        details["mtx_trio_found"] = sorted(present)
        details["n_cells_file"] = "barcodes.tsv" in files or "barcodes.tsv.gz" in files
        details["n_genes_file"] = "features.tsv" in files or "features.tsv.gz" in files
        logger.info("Detected 10X MTX format in directory '%s'", p)
        return Detection(fmt="mtx", path=p, details=details)

    raise DetectionError(
        f"Unrecognized directory format at '{p}'. "
        f"Expected a 10X directory with matrix.mtx, barcodes.tsv, features.tsv. "
        f"Found files: {sorted(files)[:20]}"
    )


def _detect_file(p: Path, details: dict) -> Detection:
    result = _detect_by_extension(p, details)
    if result is not None:
        return result

    raise DetectionError(
        f"Unrecognized file format at '{p}'. "
        f"Supported: {', '.join(EXTENSION_MAP)}, {', '.join(COMPRESSED_EXTENSIONS)}."
    )


def _detect_by_extension(p: Path, details: dict) -> Detection | None:
    name_lower = p.name.lower()

    for ext, fmt in COMPRESSED_EXTENSIONS.items():
        if name_lower.endswith(ext):
            details["extension"] = ext
            logger.info("Detected format '%s' from extension '%s'", fmt, ext)
            return Detection(fmt=fmt, path=p, details=details)

    for ext, fmt in EXTENSION_MAP.items():
        if name_lower.endswith(ext):
            details["extension"] = ext
            logger.info("Detected format '%s' from extension '%s'", fmt, ext)
            return Detection(fmt=fmt, path=p, details=details)

    return None


def _peek_h5(p: Path) -> str | None:
    try:
        import h5py
        with h5py.File(p, "r") as f:
            keys = set(f.keys())
            if "X" in keys:
                if "obs" in keys and "var" in keys:
                    return "anndata"
            if "assays" in keys:
                return "seurat_h5"
            return f"keys: {sorted(keys)[:10]}"
    except Exception:
        return None
