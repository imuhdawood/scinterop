from __future__ import annotations

import gzip
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse as sp
from scipy.io import mmread, mmwrite

from .errors import MtxAdapterError
from .schema import CanonicalObject
from .validate import assert_valid

logger = logging.getLogger(__name__)


def read_mtx(path: str | Path) -> CanonicalObject:
    path = Path(path)

    if path.is_dir():
        return _read_mtx_directory(path)
    elif path.suffix.lower() == ".mtx" or path.name.endswith(".mtx.gz"):
        return _read_mtx_file(path)
    else:
        raise MtxAdapterError(
            f"Expected a .mtx file or a directory with 10X trio, got: {path}"
        )


def write_mtx(obj: CanonicalObject, path: str | Path) -> str:
    assert_valid(obj)
    path = Path(path)

    if not path.exists():
        path.mkdir(parents=True, exist_ok=True)
    elif not path.is_dir():
        raise MtxAdapterError(f"Output path exists and is not a directory: {path}")

    n_cells, n_genes = obj.X.shape

    barcodes = obj.obs.index.tolist() if len(obj.obs) > 0 else [f"cell_{i}" for i in range(n_cells)]
    features = obj.var.index.tolist() if len(obj.var) > 0 else [f"gene_{i}" for i in range(n_genes)]

    mtx_path = path / "matrix.mtx"
    barcodes_path = path / "barcodes.tsv"
    features_path = path / "features.tsv"

    logger.info("Writing MTX to directory: %s", path)

    try:
        # 10X convention: matrix is genes x cells in the file
        X = obj.X.T
        if not sp.issparse(X):
            X = sp.csr_matrix(X)
        mmwrite(mtx_path, X, comment="scinterop-generated MTX")
        _write_tsv(barcodes_path, barcodes)
        _write_feature_tsv(features_path, features)
    except Exception as e:
        raise MtxAdapterError(f"Failed to write MTX files to '{path}': {e}") from e

    logger.info("Wrote MTX with %d cells x %d genes to %s", n_cells, n_genes, path)
    return str(path)


def _read_mtx_directory(path: Path) -> CanonicalObject:
    mtx_file = _find_first(path, "matrix.mtx")
    barcodes_file = _find_first(path, "barcodes.tsv")
    features_file = _find_first(path, "features.tsv")

    missing = []
    if mtx_file is None:
        missing.append("matrix.mtx[.gz]")
    if barcodes_file is None:
        missing.append("barcodes.tsv[.gz]")
    if features_file is None:
        missing.append("features.tsv[.gz]")
    if missing:
        raise MtxAdapterError(
            f"10X MTX directory '{path}' is missing: {', '.join(missing)}"
        )

    logger.info("Reading 10X MTX directory: %s", path)
    logger.debug("  matrix: %s", mtx_file)
    logger.debug("  barcodes: %s", barcodes_file)
    logger.debug("  features: %s", features_file)

    try:
        X = _read_mm_file(mtx_file)
    except Exception as e:
        raise MtxAdapterError(
            f"Failed to read matrix file '{mtx_file}': {e}"
        ) from e

    try:
        barcodes = _read_tsv(barcodes_file)
    except Exception as e:
        raise MtxAdapterError(
            f"Failed to read barcodes file '{barcodes_file}': {e}"
        ) from e

    try:
        features = _read_feature_tsv(features_file)
    except Exception as e:
        raise MtxAdapterError(
            f"Failed to read features file '{features_file}': {e}"
        ) from e

    n_cells = len(barcodes)
    n_genes = len(features)

    # 10X convention: file stores genes x cells → transpose to cells x genes
    X = X.T.copy()

    obs = pd.DataFrame(index=barcodes)
    var = pd.DataFrame(index=features)

    if sp.issparse(X):
        X = X.tocsr()

    obj = CanonicalObject(X=X, obs=obs, var=var)
    logger.info("Loaded MTX: %d cells x %d genes", obj.n_cells, obj.n_genes)
    return obj


def _read_mtx_file(path: Path) -> CanonicalObject:
    logger.info("Reading single MTX file: %s", path)
    try:
        X = _read_mm_file(path)
    except Exception as e:
        raise MtxAdapterError(f"Failed to read MTX file '{path}': {e}") from e

    if sp.issparse(X):
        X = X.tocsr()

    n_cells, n_genes = X.shape
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])

    obj = CanonicalObject(X=X, obs=obs, var=var)
    logger.info("Loaded MTX: %d cells x %d genes", obj.n_cells, obj.n_genes)
    return obj


def _read_mm_file(path: Path) -> sp.spmatrix:
    if str(path).endswith(".gz"):
        with gzip.open(path, "rb") as f:
            return mmread(f)
    return mmread(str(path))


def _find_first(directory: Path, base_name: str) -> Path | None:
    for candidate in [directory / base_name, directory / f"{base_name}.gz"]:
        if candidate.exists():
            return candidate
    return None


def _read_tsv(path: Path) -> list[str]:
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return [line.strip().split("\t")[0] for line in f if line.strip()]
    with open(path, "r") as f:
        return [line.strip().split("\t")[0] for line in f if line.strip()]


def _read_feature_tsv(path: Path) -> list[str]:
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return [line.strip().split("\t")[0] for line in f if line.strip()]
    with open(path, "r") as f:
        return [line.strip().split("\t")[0] for line in f if line.strip()]


def _write_tsv(path: Path, values: list[str]) -> None:
    with open(path, "w") as f:
        for v in values:
            f.write(f"{v}\n")
    logger.debug("Wrote %d lines to %s", len(values), path)


def _write_feature_tsv(path: Path, values: list[str]) -> None:
    with open(path, "w") as f:
        for v in values:
            f.write(f'{v}\t{v}\t"Gene Expression"\n')
    logger.debug("Wrote %d features to %s", len(values), path)
