# src/open_pacmuci/tools.py
"""Subprocess helpers for running external bioinformatics tools."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path


def _clean_path_for_externals(path_str: str) -> str:
    """Remove Python virtualenv bin directories from PATH.

    External tools like Clair3 internally call ``python3`` and expect the
    system or conda Python, not the project's virtualenv.  When running
    via ``uv run``, the ``.venv/bin`` directory is prepended to PATH,
    causing Clair3 to pick up the wrong Python (e.g. 3.11 instead of
    the conda env's 3.9 where tensorflow is installed).

    Uses two detection strategies (see CPython bugs.python.org/issue42041):

    1. **VIRTUAL_ENV env var:** If set, strip that exact ``bin`` directory.
    2. **Pattern matching:** Strip any PATH entry containing ``.venv/bin``
       or ``venv/bin`` as a fallback (catches cases where VIRTUAL_ENV is
       not set but a venv is on PATH).

    Conda environments are preserved — they use ``envs/`` paths, not
    ``venv`` or ``.venv``.
    """
    venv_dir = os.environ.get("VIRTUAL_ENV", "")
    venv_bin = str(Path(venv_dir) / "bin") if venv_dir else ""

    parts = path_str.split(os.pathsep)
    cleaned = [
        p
        for p in parts
        if not (
            (venv_bin and os.path.normpath(p) == os.path.normpath(venv_bin))
            or os.sep + ".venv" + os.sep in p + os.sep
            or os.sep + "venv" + os.sep in p + os.sep
        )
    ]
    return os.pathsep.join(cleaned)


def run_tool(cmd: list[str], cwd: str | None = None) -> str:
    """Run an external tool and return its stdout.

    Strips virtualenv ``bin`` directories from PATH so that external
    tools (Clair3, bcftools, samtools, etc.) use the system or conda
    Python rather than the project's virtualenv Python.

    Args:
        cmd: Command and arguments as a list.
        cwd: Optional working directory.

    Returns:
        Captured stdout as a string.

    Raises:
        FileNotFoundError: If the command is not found.
        RuntimeError: If the command exits with non-zero status.
    """
    env = os.environ.copy()
    env["PATH"] = _clean_path_for_externals(env.get("PATH", ""))

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
            cwd=cwd,
            env=env,
        )
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"Tool not found: {cmd[0]}") from exc

    if result.returncode != 0:
        raise RuntimeError(
            f"Command {' '.join(cmd)} failed with exit code {result.returncode}.\n"
            f"stderr: {result.stderr}"
        )

    return result.stdout


def check_tools(tools: list[str]) -> bool:
    """Verify that all required tools are available on PATH.

    Args:
        tools: List of tool names to check.

    Returns:
        True if all tools are found.

    Raises:
        RuntimeError: If any tools are missing, listing them all.
    """
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        raise RuntimeError(
            f"Required tools not found: {', '.join(missing)}. "
            f"Install them or activate the conda environment."
        )
    return True
