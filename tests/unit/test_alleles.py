"""Tests for allele length detection from idxstats output."""

from __future__ import annotations

import json

import pytest

from open_pacmuci.alleles import (
    PRE_AFTER_REPEAT_COUNT,
    detect_alleles,
    parse_idxstats,
)

# Example idxstats output: contig_name\tcontig_length\tmapped_reads\tunmapped_reads
# contig_N means N canonical X repeats; total allele length = N + 9 (pre+after repeats)
IDXSTATS_TWO_PEAKS = """\
contig_58\t4000\t5\t0
contig_59\t4060\t12\t0
contig_60\t4120\t245\t0
contig_61\t4180\t18\t0
contig_62\t4240\t3\t0
contig_78\t5200\t8\t0
contig_79\t5260\t15\t0
contig_80\t5320\t189\t0
contig_81\t5380\t10\t0
contig_82\t5440\t2\t0
*\t0\t0\t50
"""

IDXSTATS_HOMOZYGOUS = """\
contig_59\t4060\t20\t0
contig_60\t4120\t310\t0
contig_61\t4180\t25\t0
*\t0\t0\t10
"""

IDXSTATS_LOW_COVERAGE = """\
contig_60\t4120\t5\t0
contig_80\t5320\t3\t0
*\t0\t0\t100
"""

# Plateau distribution (like real HiFi data -- reads spread across adjacent contigs)
IDXSTATS_PLATEAU = """\
contig_48\t3000\t30\t0
contig_49\t3060\t115\t0
contig_50\t3120\t115\t0
contig_51\t3180\t115\t0
contig_52\t3240\t115\t0
contig_53\t3300\t113\t0
contig_54\t3360\t85\t0
contig_68\t4200\t15\t0
contig_69\t4260\t53\t0
contig_70\t4320\t53\t0
contig_71\t4380\t53\t0
contig_72\t4440\t53\t0
contig_73\t4500\t49\t0
contig_74\t4560\t38\t0
*\t0\t0\t20
"""


class TestParseIdxstats:
    """Tests for parsing samtools idxstats output."""

    def test_parses_contig_counts(self):
        """Extracts repeat count -> read count mapping."""
        counts = parse_idxstats(IDXSTATS_TWO_PEAKS)
        assert counts[60] == 245
        assert counts[80] == 189

    def test_ignores_unmapped_line(self):
        """The '*' unmapped line is excluded."""
        counts = parse_idxstats(IDXSTATS_TWO_PEAKS)
        # All keys should be integers (no '*' string key)
        assert all(isinstance(k, int) for k in counts)

    def test_extracts_repeat_count_from_name(self):
        """Contig name 'contig_60' maps to repeat count 60."""
        counts = parse_idxstats(IDXSTATS_TWO_PEAKS)
        assert 60 in counts
        assert counts[60] == 245


class TestDetectAlleles:
    """Tests for allele detection via peak finding."""

    def test_two_alleles_detected(self):
        """Detects two distinct allele lengths from two-peak data."""
        counts = parse_idxstats(IDXSTATS_TWO_PEAKS)
        result = detect_alleles(counts, min_coverage=10)
        lengths = sorted([result["allele_1"]["length"], result["allele_2"]["length"]])
        # contig_60 + 9 = 69, contig_80 + 9 = 89
        assert lengths == [60 + PRE_AFTER_REPEAT_COUNT, 80 + PRE_AFTER_REPEAT_COUNT]

    def test_homozygous_detection(self):
        """Single dominant peak is reported as homozygous."""
        counts = parse_idxstats(IDXSTATS_HOMOZYGOUS)
        result = detect_alleles(counts, min_coverage=10)
        assert result["homozygous"] is True
        assert result["allele_1"]["length"] == 60 + PRE_AFTER_REPEAT_COUNT

    def test_low_coverage_raises(self):
        """Raises ValueError when no contig meets minimum coverage."""
        counts = parse_idxstats(IDXSTATS_LOW_COVERAGE)
        with pytest.raises(ValueError, match="coverage"):
            detect_alleles(counts, min_coverage=10)

    def test_allele_reads_reported(self):
        """Reports read counts for each allele."""
        counts = parse_idxstats(IDXSTATS_TWO_PEAKS)
        result = detect_alleles(counts, min_coverage=10)
        # allele_1 should be the one with more reads
        assert result["allele_1"]["reads"] >= result["allele_2"]["reads"]

    def test_output_is_json_serializable(self):
        """Result can be serialized to JSON."""
        counts = parse_idxstats(IDXSTATS_TWO_PEAKS)
        result = detect_alleles(counts, min_coverage=10)
        serialized = json.dumps(result)
        assert isinstance(serialized, str)

    def test_plateau_cluster_detection(self):
        """Detects peaks from plateau distributions (real HiFi-like data)."""
        counts = parse_idxstats(IDXSTATS_PLATEAU)
        result = detect_alleles(counts, min_coverage=10)
        lengths = sorted([result["allele_1"]["length"], result["allele_2"]["length"]])
        # Cluster 1 centers around contig_51 -> 51+9=60
        # Cluster 2 centers around contig_71 -> 71+9=80
        # Allow ±2 tolerance for weighted center rounding
        assert abs(lengths[0] - 60) <= 2, f"Expected ~60, got {lengths[0]}"
        assert abs(lengths[1] - 80) <= 2, f"Expected ~80, got {lengths[1]}"
        assert result["homozygous"] is False
