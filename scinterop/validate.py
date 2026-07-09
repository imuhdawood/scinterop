from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse as sp

from .errors import ValidationError
from .schema import CanonicalObject

logger = logging.getLogger(__name__)


@dataclass
class ValidationResult:
    valid: bool
    errors: list[str] = field(default_factory=list)


def validate(obj: Any, require_raw: bool = False) -> ValidationResult:
    errors: list[str] = []

    if not isinstance(obj, CanonicalObject):
        errors.append(f"Expected CanonicalObject, got {type(obj).__name__}")
        return ValidationResult(valid=False, errors=errors)

    if not isinstance(obj.X, (np.ndarray, sp.spmatrix)):
        errors.append(f"X must be numpy array or sparse matrix, got {type(obj.X)}")

    if isinstance(obj.X, np.ndarray):
        if obj.X.ndim != 2:
            errors.append(f"X must be 2-dimensional, got shape {obj.X.shape}")
    elif sp.issparse(obj.X):
        if obj.X.ndim != 2:
            errors.append(f"X must be 2-dimensional, got shape {obj.X.shape}")

    if len(errors) > 0:
        return ValidationResult(valid=False, errors=errors)

    n_cells, n_genes = obj.X.shape

    if not isinstance(obj.obs, pd.DataFrame):
        errors.append(f"obs must be a DataFrame, got {type(obj.obs)}")
    elif len(obj.obs) > 0 and len(obj.obs) != n_cells:
        errors.append(
            f"obs has {len(obj.obs)} rows but X has {n_cells} cells"
        )

    if not isinstance(obj.var, pd.DataFrame):
        errors.append(f"var must be a DataFrame, got {type(obj.var)}")
    elif len(obj.var) > 0 and len(obj.var) != n_genes:
        errors.append(
            f"var has {len(obj.var)} rows but X has {n_genes} genes"
        )

    if not isinstance(obj.obsm, dict):
        errors.append(f"obsm must be a dict, got {type(obj.obsm)}")
    else:
        for key, val in obj.obsm.items():
            if not isinstance(val, np.ndarray):
                errors.append(f"obsm['{key}'] must be a numpy array, got {type(val)}")
            elif val.ndim != 2:
                errors.append(f"obsm['{key}'] must be 2D, got shape {val.shape}")
            elif val.shape[0] != n_cells:
                errors.append(
                    f"obsm['{key}'] has {val.shape[0]} rows but X has {n_cells} cells"
                )

    if not isinstance(obj.layers, dict):
        errors.append(f"layers must be a dict, got {type(obj.layers)}")
    else:
        for key, val in obj.layers.items():
            if not isinstance(val, (np.ndarray, sp.spmatrix)):
                errors.append(
                    f"layers['{key}'] must be array or sparse, got {type(val)}"
                )
            elif val.shape != (n_cells, n_genes):
                errors.append(
                    f"layers['{key}'] shape {val.shape} != X shape {(n_cells, n_genes)}"
                )

    if not isinstance(obj.uns, dict):
        errors.append(f"uns must be a dict, got {type(obj.uns)}")

    if require_raw and obj.raw is None:
        errors.append("raw is required but is None")

    if obj.raw is not None:
        raw_result = validate(obj.raw)
        if not raw_result.valid:
            errors.append(f"raw object is invalid: {raw_result.errors}")

    return ValidationResult(valid=len(errors) == 0, errors=errors)


def assert_valid(obj: Any, require_raw: bool = False) -> None:
    result = validate(obj, require_raw=require_raw)
    if not result.valid:
        msg = "; ".join(result.errors)
        logger.error("Validation failed: %s", msg)
        raise ValidationError(f"CanonicalObject validation failed: {msg}")
    logger.info("CanonicalObject validation passed")
