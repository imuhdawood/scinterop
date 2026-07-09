from __future__ import annotations

import logging
import os
import shutil
import uuid
from pathlib import Path
from dataclasses import dataclass

from .errors import CacheError

logger = logging.getLogger(__name__)


@dataclass
class ScratchContext:
    root: Path
    prefix: str
    uid: str

    @property
    def path(self) -> Path:
        return self.root


class ScratchManager:
    def __init__(self, base: str | Path | None = None):
        if base is None:
            base = os.environ.get(
                "SCINTEROP_SCRATCH",
                Path("/tmp") / f"scinterop_{os.environ.get('USER', 'unknown')}",
            )
        self._base = Path(base)
        self._contexts: list[ScratchContext] = []
        self._keep = os.environ.get("SCINTEROP_DEBUG", "").lower() in ("1", "true", "yes")

    def tempdir(self, prefix: str = "scint") -> ScratchContext:
        uid = uuid.uuid4().hex[:12]
        path = self._base / f"{prefix}_{uid}"
        try:
            path.mkdir(parents=True, exist_ok=False)
        except FileExistsError:
            uid = uuid.uuid4().hex[:12]
            path = self._base / f"{prefix}_{uid}"
            path.mkdir(parents=True, exist_ok=False)
        except OSError as e:
            raise CacheError(f"Failed to create scratch dir '{path}': {e}") from e

        ctx = ScratchContext(root=path, prefix=prefix, uid=uid)
        self._contexts.append(ctx)
        logger.info("Created scratch directory: %s", path)
        return ctx

    def cleanup(self, ctx: ScratchContext, *, keep: bool | None = None) -> None:
        if keep is None:
            keep = self._keep
        if keep:
            logger.info("Keeping scratch directory (debug mode): %s", ctx.path)
            return
        if ctx.path.exists():
            try:
                shutil.rmtree(ctx.path)
                logger.info("Removed scratch directory: %s", ctx.path)
            except OSError as e:
                logger.warning("Failed to remove scratch directory '%s': %s", ctx.path, e)

    def cleanup_all(self, *, keep: bool | None = None) -> None:
        for ctx in list(self._contexts):
            self.cleanup(ctx, keep=keep)
        self._contexts.clear()

    def write_script(self, ctx: ScratchContext, content: str, name: str = "script") -> Path:
        path = ctx.path / f"{name}.R"
        try:
            path.write_text(content)
            logger.debug("Wrote script to %s", path)
            return path
        except OSError as e:
            raise CacheError(
                f"Failed to write script '{path}': {e}"
            ) from e
