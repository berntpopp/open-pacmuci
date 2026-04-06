"""Integration tests for read mapping with minimap2 + samtools."""

from __future__ import annotations

from pathlib import Path

import pytest

from open_pacmuci.mapping import map_reads
from tests.conftest import requires_minimap2, requires_samtools


@requires_minimap2
@requires_samtools
@pytest.mark.integration
class TestMapReads:
    """Tests for the minimap2 mapping pipeline."""

    def test_maps_fastq_to_bam(self, tmp_path: Path) -> None:
        """Mapping a FASTQ input produces a sorted, indexed BAM."""
        ref_fa = tmp_path / "ref.fa"
        ref_fa.write_text(">contig_1\nACGTACGTACGT\n")

        reads_fq = tmp_path / "reads.fq"
        reads_fq.write_text("@read1\nACGTACGT\n+\nIIIIIIII\n")

        out_bam = map_reads(
            input_path=reads_fq,
            reference_path=ref_fa,
            output_dir=tmp_path,
        )
        assert out_bam.exists()
        assert out_bam.suffix == ".bam"
        # Index should also exist
        assert Path(str(out_bam) + ".bai").exists()

    def test_output_filename(self, tmp_path: Path) -> None:
        """Output BAM is named ``mapping.bam``."""
        ref_fa = tmp_path / "ref.fa"
        ref_fa.write_text(">contig_1\nACGTACGTACGT\n")

        reads_fq = tmp_path / "reads.fq"
        reads_fq.write_text("@read1\nACGTACGT\n+\nIIIIIIII\n")

        out_bam = map_reads(
            input_path=reads_fq,
            reference_path=ref_fa,
            output_dir=tmp_path,
        )
        assert out_bam.name == "mapping.bam"


@requires_samtools
@pytest.mark.integration
class TestBamToFastq:
    """Tests for BAM-to-FASTQ conversion."""

    def test_bam_input_converted(self, tmp_path: Path) -> None:
        """If input is BAM, it is converted to FASTQ first.

        This test is a placeholder -- full integration testing happens
        with MucOneUp test data in the e2e test suite.
        """
        pass  # Deferred to e2e tests
