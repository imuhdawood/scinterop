import logging

from . import errors
from . import cache
from . import h5ad
from . import mtx
from . import rds
from . import r_runner
from . import python_runner
from . import provenance
from .schema import CanonicalObject
from .detect import detect, Detection
from .validate import validate, assert_valid, ValidationResult

__version__ = "0.1.0"


def read(path: str, **kwargs) -> CanonicalObject:
    detected = detect(path)
    fmt = detected.fmt
    logger = _get_logger()

    logger.info("Detected format '%s' for path '%s'", fmt, path)

    if fmt == "h5ad":
        return h5ad.read_h5ad(path, **kwargs)
    elif fmt == "mtx":
        return mtx.read_mtx(path, **kwargs)
    elif fmt == "rds":
        return rds.read_rds(path, **kwargs)
    else:
        raise errors.DetectionError(
            f"Unsupported format '{fmt}' for path: {path}"
        )


def convert(
    input_path: str,
    output_path: str,
    *,
    r_exe: str = None,
    python_exe: str = None,
    seurat: bool = False,
    **kwargs,
) -> str:
    import time
    from pathlib import Path

    logger = _get_logger()
    t0 = time.time()

    obj = read(input_path, **kwargs)

    out = Path(output_path)
    suffix = _resolve_suffix(out)

    logger.info("Converting to format '%s' at '%s'", suffix, output_path)

    if suffix == "h5ad":
        h5ad.write_h5ad(obj, output_path, **kwargs)
    elif suffix == "mtx":
        mtx.write_mtx(obj, output_path, **kwargs)
    elif suffix == "rds":
        rds.write_rds(obj, output_path, seurat=seurat, r_exe=r_exe, **kwargs)
    else:
        raise errors.FormatAdapterError(
            f"Unsupported output format for extension '{suffix}': {output_path}"
        )

    elapsed = time.time() - t0
    logger.info("Conversion completed in %.2f seconds", elapsed)

    provenance.write_conversion_record(
        input_path=str(input_path),
        input_format=detect(input_path).fmt,
        output_path=str(output_path),
        output_format=suffix,
        success=True,
        runtime_s=elapsed,
    )

    return str(output_path)


def _resolve_suffix(path):
    ext = path.suffix.lower()
    if ext == ".h5ad":
        return "h5ad"
    elif ext == ".rds":
        return "rds"
    elif ext == ".gz":
        stem = path.with_suffix("").suffix.lower()
        if stem == ".mtx":
            return "mtx"
    elif ext == ".mtx":
        return "mtx"
    if not path.suffix or path.is_dir() or not path.exists():
        return "mtx"
    raise errors.FormatAdapterError(
        f"Cannot determine output format from path: {path}"
    )


def run_r(script, r_exe=None, **kwargs):
    return r_runner.run_r(script, r_exe=r_exe, **kwargs)


def run_python(script, python_exe=None, conda_env=None, **kwargs):
    return python_runner.run_python(
        script, python_exe=python_exe, conda_env=conda_env, **kwargs
    )


def _get_logger():
    return logging.getLogger("scinterop")


__all__ = [
    "CanonicalObject",
    "Detection",
    "ValidationResult",
    "read",
    "convert",
    "detect",
    "validate",
    "assert_valid",
    "run_r",
    "run_python",
    "errors",
    "cache",
    "h5ad",
    "mtx",
    "rds",
    "r_runner",
    "python_runner",
    "provenance",
]
