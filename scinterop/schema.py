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
        return self.X.shape

    @property
    def n_cells(self) -> int:
        return self.X.shape[0]

    @property
    def n_genes(self) -> int:
        return self.X.shape[1]

    def copy(self) -> CanonicalObject:
        return deepcopy(self)

    def to_anndata(self, path: str | Path, **kwargs) -> str:
        from .h5ad import write_h5ad

        return write_h5ad(self, path, **kwargs)

    def to_seurat(self, path: str | Path, *, r_exe: str | None = None) -> str:
        from .rds import write_rds

        return write_rds(self, path, seurat=True, r_exe=r_exe)

    def to_rds(self, path: str | Path, *, seurat: bool = False, r_exe: str | None = None) -> str:
        from .rds import write_rds

        return write_rds(self, path, seurat=seurat, r_exe=r_exe)

    def to_mtx(self, path: str | Path) -> str:
        from .mtx import write_mtx

        return write_mtx(self, path)
