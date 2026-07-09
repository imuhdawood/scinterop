from __future__ import annotations

import tempfile

import numpy as np
import pandas as pd
from scipy import sparse as sp

from scinterop import CanonicalObject, h5ad, errors


class TestH5adIO:
    def _require_anndata(self):
        try:
            import anndata  # noqa: F401
        except ImportError:
            import pytest
            pytest.skip("anndata not installed")

    def test_write_and_read_roundtrip(self):
        self._require_anndata()
        np.random.seed(42)
        X = sp.csr_matrix(np.random.poisson(0.5, size=(10, 5)).astype(np.float64))
        obs = pd.DataFrame({"cluster": ["A", "B"] * 5}, index=[f"cell_{i}" for i in range(10)])
        var = pd.DataFrame(index=[f"gene_{i}" for i in range(5)])
        obj = CanonicalObject(X=X, obs=obs, var=var)

        tmpfile = tempfile.mktemp(suffix=".h5ad")
        h5ad.write_h5ad(obj, tmpfile)
        obj2 = h5ad.read_h5ad(tmpfile)

        assert obj2.shape == (10, 5)
        assert np.allclose(obj.X.toarray(), obj2.X.toarray())
        assert list(obj2.obs["cluster"]) == list(obs["cluster"])

    def test_with_obsm_and_uns(self):
        self._require_anndata()
        X = np.random.rand(10, 5)
        obsm = {"X_pca": np.random.rand(10, 3)}
        uns = {"key": "value"}
        obj = CanonicalObject(X=X, obsm=obsm, uns=uns)

        tmpfile = tempfile.mktemp(suffix=".h5ad")
        h5ad.write_h5ad(obj, tmpfile)
        obj2 = h5ad.read_h5ad(tmpfile)

        assert "X_pca" in obj2.obsm
        assert obj2.obsm["X_pca"].shape == (10, 3)
        assert obj2.uns.get("key") == "value"

    def test_with_raw(self):
        self._require_anndata()
        X = np.random.rand(10, 5)
        raw = CanonicalObject(X=np.random.rand(10, 10))
        obj = CanonicalObject(X=X, raw=raw)

        tmpfile = tempfile.mktemp(suffix=".h5ad")
        h5ad.write_h5ad(obj, tmpfile)
        obj2 = h5ad.read_h5ad(tmpfile)

        assert obj2.raw is not None
        assert obj2.raw.shape == (10, 10)

    def test_read_nonexistent(self):
        try:
            h5ad.read_h5ad("/nonexistent/file.h5ad")
            assert False, "Should have raised"
        except errors.H5adAdapterError:
            pass
