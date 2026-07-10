"""R subprocess execution for Seurat-based conversions.

Runs R scripts (file paths or inline strings) via ``Rscript``
with configurable timeout, working directory, and environment.
"""

from __future__ import annotations

import logging
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from .errors import RExecError

logger = logging.getLogger(__name__)


@dataclass
class RExecResult:
    """Result of an R subprocess execution.

    Attributes:
        returncode: Process exit code.
        stdout: Captured standard output text.
        stderr: Captured standard error text.
        log_path: Path to the log file (if ``log_path`` was set).
    """

    returncode: int
    stdout: str
    stderr: str
    log_path: Path | None = None


def resolve_r_exe(r_exe: str | None = None) -> str:
    """Resolve the Rscript executable path.

    Priority: explicit argument, ``SCINTEROP_R_EXE`` env var,
    ``"Rscript"`` fallback.

    Args:
        r_exe: Explicit path to Rscript.

    Returns:
        Resolved Rscript path.
    """
    if r_exe is not None:
        return r_exe
    env_r = os.environ.get("SCINTEROP_R_EXE")
    if env_r:
        return env_r
    return "Rscript"


def run_r(
    script: str | Path,
    r_exe: str | None = None,
    args: list[str] | None = None,
    log_path: str | Path | None = None,
    workdir: str | Path | None = None,
    timeout: int | None = 600,
    env: dict[str, str] | None = None,
) -> RExecResult:
    """Run an R script or snippet in a subprocess.

    Inline scripts (strings containing newlines) are written to
    a temporary ``.R`` file first.

    Args:
        script: Path to a ``.R`` file, or inline R code.
        r_exe: Rscript path (see :func:`resolve_r_exe`).
        args: Additional arguments passed to the script.
        log_path: If set, write stdout/stderr to this file.
        workdir: Working directory for the subprocess.
        timeout: Maximum wall-clock seconds before killing (default 600).
        env: Extra environment variables for the subprocess.

    Returns:
        An :class:`RExecResult` with the process output.

    Raises:
        RExecError: If the script fails, times out, or the
            executable is not found.
    """
    resolved_r = resolve_r_exe(r_exe)
    logger.info("Running R script with: %s", resolved_r)

    cmd = [resolved_r]

    script_path: Path
    if isinstance(script, str) and "\n" in script:
        import tempfile
        tmp = Path(tempfile.mktemp(suffix=".R", prefix="scinterop_r_"))
        try:
            tmp.write_text(script)
        except OSError as e:
            raise RExecError(f"Failed to write temp R script: {e}") from e
        script_path = tmp
    else:
        script_path = Path(script)

    cmd.append(str(script_path))

    if args:
        cmd.extend(args)

    run_env = os.environ.copy()
    if env:
        run_env.update(env)

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=str(workdir) if workdir else None,
            timeout=timeout,
            env=run_env,
        )
    except subprocess.TimeoutExpired:
        raise RExecError(
            f"R script timed out after {timeout}s: {script_path}"
        ) from None
    except FileNotFoundError:
        raise RExecError(
            f"R executable not found: '{resolved_r}'. "
            "Set SCINTEROP_R_EXE environment variable or ensure Rscript is on PATH."
        ) from None
    except OSError as e:
        raise RExecError(f"Failed to run R script '{script_path}': {e}") from e

    log_p = None
    if log_path:
        log_p = Path(log_path)
        try:
            log_p.parent.mkdir(parents=True, exist_ok=True)
            with open(log_p, "w") as f:
                f.write("STDOUT:\n")
                f.write(result.stdout)
                f.write("\nSTDERR:\n")
                f.write(result.stderr)
        except OSError as e:
            logger.warning("Failed to write R log to '%s': %s", log_p, e)

    if result.returncode != 0:
        msg = (
            f"R script failed with code {result.returncode}.\n"
            f"Script: {script_path}\n"
            f"stderr:\n{result.stderr[:2000]}"
        )
        logger.error(msg)
        raise RExecError(msg)

    logger.info("R script completed successfully (code 0)")
    return RExecResult(
        returncode=result.returncode,
        stdout=result.stdout,
        stderr=result.stderr,
        log_path=log_p,
    )
