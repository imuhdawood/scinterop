from __future__ import annotations

import logging
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .errors import PythonExecError

logger = logging.getLogger(__name__)


@dataclass
class PythonExecResult:
    returncode: int
    stdout: str
    stderr: str
    log_path: Path | None = None


def resolve_python_exe(python_exe: str | None = None) -> str:
    if python_exe is not None:
        return python_exe
    env_py = os.environ.get("SCINTEROP_PYTHON_EXE")
    if env_py:
        return env_py
    return "python"


def run_python(
    script: str | Path,
    python_exe: str | None = None,
    conda_env: str | None = None,
    args: list[str] | None = None,
    log_path: str | Path | None = None,
    workdir: str | Path | None = None,
    timeout: int | None = 600,
    env: dict[str, str] | None = None,
) -> PythonExecResult:
    run_env = os.environ.copy()
    if env:
        run_env.update(env)

    if conda_env:
        conda_exe = os.environ.get("SCINTEROP_CONDA_EXE", "conda")
        resolved_python = conda_exe
        base_cmd = [conda_exe, "run", "-n", conda_env, "python"]
        logger.info("Running Python via conda env '%s'", conda_env)
    else:
        resolved_python = resolve_python_exe(python_exe)
        base_cmd = [resolved_python]
        logger.info("Running Python with: %s", resolved_python)

    cmd = list(base_cmd)

    if isinstance(script, str) and ("\n" in script or not script.endswith(".py")):
        import tempfile
        tmp = Path(tempfile.mktemp(suffix=".py", prefix="scinterop_py_"))
        try:
            tmp.write_text(script)
        except OSError as e:
            raise PythonExecError(f"Failed to write temp Python script: {e}") from e
        cmd.append(str(tmp))
    else:
        cmd.append(str(Path(script)))

    if args:
        cmd.extend(args)

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
        raise PythonExecError(
            f"Python script timed out after {timeout}s"
        ) from None
    except FileNotFoundError:
        exe_hint = conda_env or resolved_python
        raise PythonExecError(
            f"Python executable not found: '{exe_hint}'. "
            "Set SCINTEROP_PYTHON_EXE or ensure python is on PATH."
        ) from None
    except OSError as e:
        raise PythonExecError(f"Failed to run Python script: {e}") from e

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
            logger.warning("Failed to write Python log to '%s': %s", log_p, e)

    if result.returncode != 0:
        msg = (
            f"Python script failed with code {result.returncode}.\n"
            f"stderr:\n{result.stderr[:2000]}"
        )
        logger.error(msg)
        raise PythonExecError(msg)

    logger.info("Python script completed successfully (code 0)")
    return PythonExecResult(
        returncode=result.returncode,
        stdout=result.stdout,
        stderr=result.stderr,
        log_path=log_p,
    )
