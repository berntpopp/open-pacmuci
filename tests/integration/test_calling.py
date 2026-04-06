"""Integration tests for variant calling with Clair3 + bcftools."""

from __future__ import annotations

import pytest

from tests.conftest import requires_bcftools, requires_clair3, requires_samtools


@requires_samtools
@pytest.mark.integration
class TestExtractAlleleReads:
    """Tests for extracting reads mapped to a specific contig."""

    def test_extract_placeholder(self) -> None:
        """Placeholder -- tested with real data in e2e tests.

        Full integration testing requires a BAM file with reads mapped to
        named contigs, which is provided by the MucOneUp test data set.
        """
        pass  # Deferred to e2e tests


@requires_bcftools
@pytest.mark.integration
class TestFilterVcf:
    """Tests for VCF filtering with bcftools."""

    def test_filter_placeholder(self) -> None:
        """Placeholder -- tested with real data in e2e tests."""
        pass  # Deferred to e2e tests


@requires_clair3
@pytest.mark.integration
class TestRunClair3:
    """Tests for the Clair3 variant caller wrapper."""

    def test_clair3_placeholder(self) -> None:
        """Placeholder -- tested with real data in e2e tests."""
        pass  # Deferred to e2e tests
