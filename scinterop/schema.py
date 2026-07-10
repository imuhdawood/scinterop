"""Core data model — the CanonicalObject.

A format-agnostic in-memory representation of a single-cell dataset.
Maps cleanly to AnnData, Seurat, and 10X MTX structures.
"""

from __future__ import annotations

import logging
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse as sp

from .errors import ValidationError

logger = logging.getLogger(__name__)


@dataclass
class CanonicalObject:
    """In-memory container for a single-cell dataset.

    Mirrors the AnnData structure: an expression matrix ``X``,
    cell metadata ``obs``, gene metadata ``var``, multidimensional
    reductions ``obsm``, expression layers ``layers``, unstructured
    metadata ``uns``, and an optional ``raw`` counts object.

    Args:
        X: Expression matrix, cells x genes.
        obs: Cell-level metadata DataFrame.
        var: Gene-level metadata DataFrame.
        obsm: Multidimensional cell embeddings (e.g. PCA, UMAP).
        layers: Named expression layers (same shape as X).
        uns: Unstructured metadata dictionary.
        raw: Raw counts as a nested CanonicalObject.

    Raises:
        ValidationError: If ``obs``/``var`` are not DataFrames,
            their lengths mismatch ``X``, or ``raw`` is not a
            CanonicalObject.
    """

    X: np.ndarray | sp.spmatrix
    obs: pd.DataFrame = field(default_factory=lambda: pd.DataFrame())
    var: pd.DataFrame = field(default_factory=lambda: pd.DataFrame())
    obsm: dict[str, np.ndarray] = field(default_factory=dict)
    layers: dict[str, np.ndarray | sp.spmatrix] = field(default_factory=dict)
    uns: dict[str, Any] = field(default_factory=dict)
    raw: CanonicalObject | None = None

    def __post_init__(self):
        if not isinstance(self.obs, pd.DataFrame):
            raise ValidationError(f"obs must be a DataFrame, got {type(self.obs)}")
        if not isinstance(self.var, pd.DataFrame):
            raise ValidationError(f"var must be a DataFrame, got {type(self.var)}")
        n_cells, n_genes = self.X.shape
        if len(self.obs) > 0 and len(self.obs) != n_cells:
            raise ValidationError(
                f"obs has {len(self.obs)} rows but X has {n_cells} cells"
            )
        if len(self.var) > 0 and len(self.var) != n_genes:
            raise ValidationError(
                f"var has {len(self.var)} rows but X has {n_genes} genes"
            )
        if self.raw is not None and not isinstance(self.raw, CanonicalObject):
            raise ValidationError(
                f"raw must be a CanonicalObject or None, got {type(self.raw)}"
            )

    @property
    def shape(self) -> tuple[int, int]:
        """Return (n_cells, n_genes) tuple."""
        return self.X.shape

    @property
    def n_cells(self) -> int:
        """Return number of cells (rows in X)."""
        return self.X.shape[0]

    @property
    def n_genes(self) -> int:
        """Return number of genes (columns in X)."""
        return self.X.shape[1]

    def copy(self) -> CanonicalObject:
        """Return a deep copy of this object."""
        return deepcopy(self)

    def to_anndata(self, path: str | Path, **kwargs: Any) -> str:
        """Write to H5AD via the h5ad adapter.

        Args:
            path: Output file path.
            **kwargs: Forwarded to :func:`h5ad.write_h5ad`.

        Returns:
            The path the file was written to.
        """
        from .h5ad import write_h5ad
        return write_h5ad(self, path, **kwargs)

    def to_seurat(self, path: str | Path, *, r_exe: str | None = None) -> str:
        """Write to Seurat RDS (via the R bridge).

        Shortcut for ``to_rds(seurat=True)``.

        Args:
            path: Output ``.rds`` file path.
            r_exe: Path to Rscript executable.

        Returns:
            The path the file was written to.
        """
        from .rds import write_rds
        return write_rds(self, path, seurat=True, r_exe=r_exe)

    def to_rds(self, path: str | Path, *, seurat: bool = False, r_exe: str | None = None) -> str:
        """Write to RDS file (plain R list or Seurat object).

        Args:
            path: Output ``.rds`` file path.
            seurat: If True, create a Seurat object instead of a plain list.
            r_exe: Path to Rscript executable.

        Returns:
            The path the file was written to.
        """
        from .rds import write_rds
        return write_rds(self, path, seurat=seurat, r_exe=r_exe)

    def to_mtx(self, path: str | Path) -> str:
        """Write to 10X MTX directory.

        Args:
            path: Output directory path.

        Returns:
            The path the directory was written to.
        """
        from .mtx import write_mtx
        return write_mtx(self, path)
