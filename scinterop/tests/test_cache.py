from __future__ import annotations

import tempfile
from pathlib import Path

from scinterop import cache, errors


class TestCache:
    def test_tempdir_creation(self):
        mgr = cache.ScratchManager(base=tempfile.mkdtemp())
        ctx = mgr.tempdir("test")
        assert ctx.path.exists()
        assert ctx.path.is_dir()

    def test_tempdir_cleanup(self):
        base = tempfile.mkdtemp()
        mgr = cache.ScratchManager(base=base)
        ctx = mgr.tempdir("test")
        path = ctx.path
        assert path.exists()
        mgr.cleanup(ctx, keep=False)
        assert not path.exists()

    def test_cleanup_keep(self):
        base = tempfile.mkdtemp()
        mgr = cache.ScratchManager(base=base)
        ctx = mgr.tempdir("test")
        path = ctx.path
        mgr.cleanup(ctx, keep=True)
        assert path.exists()
        path.rmdir()

    def test_cleanup_nonexistent(self):
        base = tempfile.mkdtemp()
        mgr = cache.ScratchManager(base=base)
        ctx = cache.ScratchContext(root=Path("/nonexistent_path"), prefix="test", uid="123")
        mgr.cleanup(ctx)

    def test_write_script(self):
        base = tempfile.mkdtemp()
        mgr = cache.ScratchManager(base=base)
        ctx = mgr.tempdir("test")
        script_path = mgr.write_script(ctx, "print('hello')", "test_script")
        assert script_path.exists()
        assert script_path.read_text() == "print('hello')"

    def test_cleanup_all(self):
        base = tempfile.mkdtemp()
        mgr = cache.ScratchManager(base=base)
        ctx1 = mgr.tempdir("test1")
        ctx2 = mgr.tempdir("test2")
        mgr.cleanup_all(keep=False)
        assert not ctx1.path.exists()
        assert not ctx2.path.exists()
