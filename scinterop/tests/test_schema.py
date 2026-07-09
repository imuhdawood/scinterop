from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse as sp

from scinterop import CanonicalObject, errors


class TestCanonicalObject:
    def test_minimal_creation(self):
        obj = CanonicalObject(X=np.random.rand(5, 3))
        assert obj.shape == (5, 3)
        assert obj.n_cells == 5
        assert obj.n_genes == 3

    def test_with_obs_var(self):
        X = np.random.rand(5, 3)
        obs = pd.DataFrame(index=[f"c{i}" for i in range(5)])
        var = pd.DataFrame(index=[f"g{i}" for i in range(3)])
        obj = CanonicalObject(X=X, obs=obs, var=var)
        assert obj.shape == (5, 3)
        assert list(obj.obs.index) == [f"c{i}" for i in range(5)]
        assert list(obj.var.index) == [f"g{i}" for i in range(3)]

    def test_sparse_matrix(self):
        X = sp.csr_matrix(np.random.rand(5, 3))
        obj = CanonicalObject(X=X)
        assert obj.shape == (5, 3)

    def test_shape_mismatch_obs(self):
        X = np.random.rand(5, 3)
        obs = pd.DataFrame(index=[f"c{i}" for i in range(10)])
        try:
            CanonicalObject(X=X, obs=obs)
            assert False, "Should have raised"
        except errors.ValidationError:
            pass

    def test_shape_mismatch_var(self):
        X = np.random.rand(5, 3)
        var = pd.DataFrame(index=[f"g{i}" for i in range(10)])
        try:
            CanonicalObject(X=X, var=var)
            assert False, "Should have raised"
        except errors.ValidationError:
            pass

    def test_with_obsm(self):
        X = np.random.rand(5, 3)
        obsm = {"X_pca": np.random.rand(5, 2)}
        obj = CanonicalObject(X=X, obsm=obsm)
        assert "X_pca" in obj.obsm
        assert obj.obsm["X_pca"].shape == (5, 2)

    def test_with_layers(self):
        X = np.random.rand(5, 3)
        layers = {"counts": np.random.rand(5, 3)}
        obj = CanonicalObject(X=X, layers=layers)
        assert "counts" in obj.layers

    def test_with_uns(self):
        X = np.random.rand(5, 3)
        uns = {"key": "value"}
        obj = CanonicalObject(X=X, uns=uns)
        assert obj.uns["key"] == "value"

    def test_copy(self):
        X = np.random.rand(5, 3)
        obj = CanonicalObject(X=X)
        obj2 = obj.copy()
        assert obj2.shape == obj.shape
        obj2.X[0, 0] = 999
        assert obj.X[0, 0] != 999

    def test_nested_raw(self):
        X = np.random.rand(5, 3)
        raw = CanonicalObject(X=np.random.rand(5, 10))
        obj = CanonicalObject(X=X, raw=raw)
        assert obj.raw is not None
        assert obj.raw.shape == (5, 10)

    def test_invalid_raw_type(self):
        X = np.random.rand(5, 3)
        try:
            CanonicalObject(X=X, raw="not_a_canonical")
            assert False, "Should have raised"
        except errors.ValidationError:
            pass
