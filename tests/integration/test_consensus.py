"""Integration tests for bcftools consensus generation."""

from __future__ import annotations

import pytest

from tests.conftest import requires_bcftools, requires_samtools


@requires_bcftools
@pytest.mark.integration
class TestBuildConsensus:
    """Tests for consensus FASTA generation."""

    def test_consensus_placeholder(self) -> None:
        """Placeholder -- tested with real data in e2e tests."""
        pass  # Deferred to e2e tests


@requires_bcftools
@requires_samtools
@pytest.mark.integration
class TestBuildConsensusPerAllele:
    """Tests for per-allele consensus generation."""

    def test_per_allele_placeholder(self) -> None:
        """Placeholder -- tested with real data in e2e tests."""
        pass  # Deferred to e2e tests
