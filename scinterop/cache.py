"""Scratch-directory management for intermediate files.

The :class:`ScratchManager` creates, caches, and cleans up temporary
directories shared across conversion steps.
"""

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
    """A single scratch directory handle.

    Attributes:
        root: The filesystem path of this scratch directory.
        prefix: Human-readable prefix used in the directory name.
        uid: Unique identifier for this context.
    """

    root: Path
    prefix: str
    uid: str

    @property
    def path(self) -> Path:
        """Return the root path of this scratch directory."""
        return self.root


class ScratchManager:
    """Manage temporary scratch directories for conversion artifacts.

    Args:
        base: Root directory for scratch dirs. Falls back to
            ``SCINTEROP_SCRATCH`` env var or ``/tmp/scinterop_<user>``.
    """

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
        """Create a new scratch subdirectory.

        Args:
            prefix: Name prefix for the directory.

        Returns:
            A :class:`ScratchContext` pointing to the new directory.

        Raises:
            CacheError: If the directory cannot be created.
        """
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
        """Remove a specific scratch directory.

        Args:
            ctx: The context to clean up.
            keep: Override the debug keep setting. If None, uses
                the instance-level ``_keep``.
        """
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
        """Remove all tracked scratch directories.

        Args:
            keep: Override the debug keep setting for all contexts.
        """
        for ctx in list(self._contexts):
            self.cleanup(ctx, keep=keep)
        self._contexts.clear()

    def write_script(self, ctx: ScratchContext, content: str, name: str = "script") -> Path:
        """Write a script file into a scratch directory.

        Args:
            ctx: Scratch context to write into.
            content: Script text content.
            name: Basename (without extension) for the script file.

        Returns:
            Path to the written script file.

        Raises:
            CacheError: If the file cannot be written.
        """
        path = ctx.path / f"{name}.R"
        try:
            path.write_text(content)
            logger.debug("Wrote script to %s", path)
            return path
        except OSError as e:
            raise CacheError(
                f"Failed to write script '{path}': {e}"
            ) from e
