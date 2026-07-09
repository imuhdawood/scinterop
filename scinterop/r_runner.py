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
    returncode: int
    stdout: str
    stderr: str
    log_path: Path | None = None


def resolve_r_exe(r_exe: str | None = None) -> str:
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
