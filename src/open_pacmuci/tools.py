# src/open_pacmuci/tools.py
"""Subprocess helpers for running external bioinformatics tools."""

from __future__ import annotations

import shutil
import subprocess


def run_tool(cmd: list[str], cwd: str | None = None) -> str:
    """Run an external tool and return its stdout.

    Args:
        cmd: Command and arguments as a list.
        cwd: Optional working directory.

    Returns:
        Captured stdout as a string.

    Raises:
        FileNotFoundError: If the command is not found.
        RuntimeError: If the command exits with non-zero status.
    """
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
            cwd=cwd,
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
