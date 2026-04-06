"""Shared test fixtures and configuration."""

from __future__ import annotations

import shutil

import pytest


def tool_available(name: str) -> bool:
    """Check if a command-line tool is on PATH."""
    return shutil.which(name) is not None


requires_minimap2 = pytest.mark.skipif(
    not tool_available("minimap2"), reason="minimap2 not installed"
)
requires_samtools = pytest.mark.skipif(
    not tool_available("samtools"), reason="samtools not installed"
)
requires_bcftools = pytest.mark.skipif(
    not tool_available("bcftools"), reason="bcftools not installed"
)
requires_clair3 = pytest.mark.skipif(
    not tool_available("run_clair3.sh"), reason="Clair3 not installed"
)
