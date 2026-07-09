from __future__ import annotations

import tempfile

import numpy as np
import pandas as pd
from scipy import sparse as sp

from scinterop import CanonicalObject, mtx, errors


class TestMtxIO:
    def test_write_and_read_roundtrip(self):
        np.random.seed(42)
        X = sp.csr_matrix(np.random.poisson(0.5, size=(10, 5)).astype(np.float64))
        obs = pd.DataFrame(index=[f"cell_{i}" for i in range(10)])
        var = pd.DataFrame(index=[f"gene_{i}" for i in range(5)])
        obj = CanonicalObject(X=X, obs=obs, var=var)

        tmpdir = tempfile.mkdtemp()
        mtx.write_mtx(obj, tmpdir)
        obj2 = mtx.read_mtx(tmpdir)

        assert obj2.shape == (10, 5)
        assert np.allclose(obj.X.toarray(), obj2.X.toarray())
        assert list(obj2.obs.index) == list(obj.obs.index)
        assert list(obj2.var.index) == list(obj.var.index)

    def test_read_nonexistent_directory(self):
        try:
            mtx.read_mtx("/nonexistent/directory")
            assert False, "Should have raised"
        except errors.MtxAdapterError:
            pass

    def test_read_single_mtx_file(self):
        np.random.seed(42)
        X = sp.csr_matrix(np.random.rand(5, 3))
        tmpdir = tempfile.mkdtemp()
        mtx_path = f"{tmpdir}/matrix.mtx"
        import scipy.io
        scipy.io.mmwrite(mtx_path, X.T)

        obj = mtx._read_mtx_file(tmpdir + "/matrix.mtx")
        assert obj.shape == (3, 5)

    def test_write_to_existing_file_path(self):
        obj = CanonicalObject(X=np.random.rand(5, 3))
        with tempfile.NamedTemporaryFile() as f:
            try:
                mtx.write_mtx(obj, f.name)
                assert False, "Should have raised"
            except errors.MtxAdapterError:
                pass
