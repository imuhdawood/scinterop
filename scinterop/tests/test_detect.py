from __future__ import annotations

import tempfile
from pathlib import Path

from scinterop import detect, errors


class TestDetect:
    def test_detect_h5ad(self):
        d = detect("dummy.h5ad")
        assert d.fmt == "h5ad"

    def test_detect_rds(self):
        d = detect("dummy.rds")
        assert d.fmt == "rds"

    def test_detect_mtx_file(self):
        d = detect("matrix.mtx")
        assert d.fmt == "mtx"

    def test_detect_mtx_gz(self):
        d = detect("matrix.mtx.gz")
        assert d.fmt == "mtx"

    def test_detect_mtx_directory(self):
        tmpdir = Path(tempfile.mkdtemp())
        (tmpdir / "matrix.mtx.gz").write_text("dummy")
        (tmpdir / "barcodes.tsv").write_text("dummy")
        (tmpdir / "features.tsv").write_text("dummy")
        d = detect(str(tmpdir))
        assert d.fmt == "mtx"

    def test_detect_nonexistent(self):
        try:
            detect("/nonexistent/path")
            assert False, "Should have raised"
        except errors.DetectionError:
            pass

    def test_detect_empty_directory(self):
        tmpdir = tempfile.mkdtemp()
        try:
            detect(tmpdir)
            assert False, "Should have raised"
        except errors.DetectionError:
            pass

    def test_detect_unknown_extension(self):
        try:
            detect("file.xyz")
            assert False, "Should have raised"
        except errors.DetectionError:
            pass

    def test_detect_details_present(self):
        d = detect("test.h5ad")
        assert d.details is not None
        assert "extension" in d.details
