from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse as sp

from .errors import H5adAdapterError
from .schema import CanonicalObject
from .validate import assert_valid

logger = logging.getLogger(__name__)


def read_h5ad(path: str | Path) -> CanonicalObject:
    path = Path(path)
    if not path.exists():
        raise H5adAdapterError(f"File does not exist: {path}")
    if not path.suffix.lower() in (".h5ad", ".h5"):
        logger.warning("File '%s' does not have .h5ad extension, attempting read", path)

    try:
        import anndata as ad
    except ImportError as e:
        raise H5adAdapterError(
            "anndata package is required to read H5AD files. "
            "Install it with: pip install anndata"
        ) from e

    logger.info("Reading H5AD file: %s", path)
    try:
        adata = ad.read_h5ad(path)
    except Exception as e:
        raise H5adAdapterError(
            f"Failed to read H5AD file '{path}': {e}"
        ) from e

    logger.info(
        "Loaded AnnData with %d cells x %d genes",
        adata.n_obs, adata.n_vars,
    )

    X = _resolve_X(adata)
    obs = adata.obs.copy() if adata.obs is not None else pd.DataFrame()
    var = adata.var.copy() if adata.var is not None else pd.DataFrame()

    obsm = dict(adata.obsm) if adata.obsm is not None else {}
    layers = dict(adata.layers) if adata.layers is not None else {}
    uns = dict(adata.uns) if adata.uns is not None else {}

    raw_obj = None
    if adata.raw is not None:
        try:
            raw_X = adata.raw.X
            if raw_X is None:
                raw_X = adata.raw[:, :].X
            if sp.issparse(raw_X):
                raw_X = raw_X.toarray()
            raw_obj = CanonicalObject(
                X=raw_X,
                obs=obs.copy(),
                var=adata.raw.var.copy() if adata.raw.var is not None else var.copy(),
            )
            logger.info("Loaded raw counts with %d genes", adata.raw.n_vars)
        except Exception as e:
            logger.warning("Failed to load raw counts: %s", e)

    obj = CanonicalObject(
        X=X,
        obs=obs,
        var=var,
        obsm=obsm,
        layers=layers,
        uns=uns,
        raw=raw_obj,
    )

    assert_valid(obj)
    return obj


def write_h5ad(obj: CanonicalObject, path: str | Path) -> str:
    path = Path(path)
    if not path.suffix.lower() == ".h5ad":
        path = path.with_suffix(".h5ad")

    assert_valid(obj)

    try:
        import anndata as ad
    except ImportError as e:
        raise H5adAdapterError(
            "anndata package is required to write H5AD files"
        ) from e

    logger.info("Writing H5AD file: %s", path)

    try:
        path.parent.mkdir(parents=True, exist_ok=True)

        obs = obj.obs
        var = obj.var
        if len(obs) == 0:
            obs = pd.DataFrame(index=[f"cell_{i}" for i in range(obj.n_cells)])
        if len(var) == 0:
            var = pd.DataFrame(index=[f"gene_{i}" for i in range(obj.n_genes)])

        adata = ad.AnnData(X=obj.X, obs=obs, var=var)

        if obj.obsm:
            adata.obsm = obj.obsm
        if obj.layers:
            adata.layers = obj.layers
        if obj.uns:
            adata.uns = obj.uns

        if obj.raw is not None:
            raw_adata = ad.AnnData(
                X=obj.raw.X,
                var=obj.raw.var if len(obj.raw.var) > 0 else None,
            )
            adata.raw = raw_adata

        adata.write(path, compression="gzip")
    except Exception as e:
        raise H5adAdapterError(
            f"Failed to write H5AD file '{path}': {e}"
        ) from e

    logger.info("Successfully wrote H5AD file: %s", path)
    return str(path)


def _resolve_X(adata) -> np.ndarray | sp.spmatrix:
    try:
        X = adata.X
    except Exception as e:
        raise H5adAdapterError(
            f"Failed to access adata.X: {e}. "
            "The file may be corrupted or use a backing store."
        ) from e
    if X is None:
        raise H5adAdapterError("adata.X is None — no expression matrix found")
    return X
