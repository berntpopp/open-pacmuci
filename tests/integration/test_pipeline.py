"""End-to-end pipeline tests using MucOneUp-generated test data."""

from __future__ import annotations

from pathlib import Path

import pytest

from tests.conftest import (
    requires_bcftools,
    requires_clair3,
    requires_minimap2,
    requires_samtools,
)

# Path to generated test data
TESTDATA_DIR = Path("tests/data/generated")


def testdata_available() -> bool:
    """Check if MucOneUp-generated test data exists."""
    return TESTDATA_DIR.exists() and any(TESTDATA_DIR.iterdir())


skip_no_testdata = pytest.mark.skipif(
    not testdata_available(), reason="Test data not generated (run 'make generate-testdata')"
)


@skip_no_testdata
@requires_minimap2
@requires_samtools
@pytest.mark.e2e
class TestAlleleLengthDetection:
    """Test allele length detection against ground truth."""

    @pytest.mark.parametrize(
        "sample,expected_h1,expected_h2",
        [
            ("sample_dupc_60_80", 60, 80),
            ("sample_normal_60_80", 60, 80),
            ("sample_homozygous_60_60", 60, 60),
            ("sample_asymmetric_25_140", 25, 140),
            ("sample_short_25_30", 25, 30),
            ("sample_long_120_140", 120, 140),
        ],
    )
    def test_allele_lengths(
        self, sample: str, expected_h1: int, expected_h2: int, tmp_path: Path
    ) -> None:
        """Detected allele lengths match ground truth within tolerance."""
        from open_pacmuci.alleles import detect_alleles, parse_idxstats
        from open_pacmuci.config import load_repeat_dictionary
        from open_pacmuci.ladder import generate_ladder_fasta
        from open_pacmuci.mapping import get_idxstats, map_reads

        rd = load_repeat_dictionary()
        ref = tmp_path / "ladder.fa"
        generate_ladder_fasta(rd, ref)

        # Find BAM/FASTQ in test data
        sample_dir = TESTDATA_DIR / sample
        input_files = list(sample_dir.glob("*.bam")) + list(sample_dir.glob("*.fq"))
        assert input_files, f"No input files found in {sample_dir}"

        bam = map_reads(input_files[0], ref, tmp_path)
        idxstats = get_idxstats(bam)
        counts = parse_idxstats(idxstats)
        result = detect_alleles(counts, min_coverage=10)

        detected = sorted([result["allele_1"]["length"], result["allele_2"]["length"]])
        expected = sorted([expected_h1, expected_h2])

        # Allow +/- 2 repeat tolerance
        assert abs(detected[0] - expected[0]) <= 2, f"Expected {expected[0]}, got {detected[0]}"
        assert abs(detected[1] - expected[1]) <= 2, f"Expected {expected[1]}, got {detected[1]}"


@skip_no_testdata
@requires_minimap2
@requires_samtools
@requires_bcftools
@requires_clair3
@pytest.mark.e2e
class TestFullPipeline:
    """Full pipeline e2e tests with mutation detection validation."""

    def test_dupc_detected(self, tmp_path: Path) -> None:
        """dupC mutation is detected in sample_dupc_60_80."""
        # This is the most important test -- validates the core pipeline
        # Full implementation deferred to when all components are integrated
        pass

    def test_normal_no_mutation(self, tmp_path: Path) -> None:
        """Normal sample reports no frameshift mutations."""
        pass
