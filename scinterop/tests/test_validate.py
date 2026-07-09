from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse as sp

from scinterop import CanonicalObject, validate, assert_valid, errors


class TestValidate:
    def test_valid_object(self):
        obj = CanonicalObject(X=np.random.rand(5, 3))
        result = validate(obj)
        assert result.valid
        assert len(result.errors) == 0

    def test_valid_with_obs_var(self):
        X = np.random.rand(5, 3)
        obj = CanonicalObject(
            X=X,
            obs=pd.DataFrame(index=[f"c{i}" for i in range(5)]),
            var=pd.DataFrame(index=[f"g{i}" for i in range(3)]),
        )
        result = validate(obj)
        assert result.valid

    def test_invalid_type(self):
        result = validate("not_an_object")
        assert not result.valid

    def test_invalid_x_type(self):
        obj = CanonicalObject.__new__(CanonicalObject)
        obj.X = "not_a_matrix"
        obj.obs = pd.DataFrame()
        obj.var = pd.DataFrame()
        obj.obsm = {}
        obj.layers = {}
        obj.uns = {}
        obj.raw = None
        result = validate(obj)
        assert not result.valid

    def test_obsm_shape_mismatch(self):
        X = np.random.rand(5, 3)
        obsm = {"X_pca": np.random.rand(10, 2)}
        obj = CanonicalObject(X=X, obsm=obsm)
        result = validate(obj)
        assert not result.valid
        assert any("obsm" in e for e in result.errors)

    def test_layers_shape_mismatch(self):
        X = np.random.rand(5, 3)
        layers = {"counts": np.random.rand(5, 10)}
        obj = CanonicalObject(X=X, layers=layers)
        result = validate(obj)
        assert not result.valid

    def test_assert_valid_passes(self):
        obj = CanonicalObject(X=np.random.rand(5, 3))
        assert_valid(obj)

    def test_assert_valid_raises(self):
        obj = CanonicalObject.__new__(CanonicalObject)
        obj.X = np.random.rand(5, 3)
        obj.obs = "not_a_dataframe"
        obj.var = pd.DataFrame()
        obj.obsm = {}
        obj.layers = {}
        obj.uns = {}
        obj.raw = None
        try:
            assert_valid(obj)
            assert False
        except errors.ValidationError:
            pass

    def test_require_raw(self):
        obj = CanonicalObject(X=np.random.rand(5, 3))
        result = validate(obj, require_raw=True)
        assert not result.valid

    def test_valid_with_raw(self):
        X = np.random.rand(5, 3)
        raw = CanonicalObject(X=np.random.rand(5, 10))
        obj = CanonicalObject(X=X, raw=raw)
        result = validate(obj, require_raw=True)
        assert result.valid
