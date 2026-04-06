"""Allele length detection from samtools idxstats output."""

from __future__ import annotations

import re


def parse_idxstats(idxstats_output: str) -> dict:
    """Parse samtools idxstats output into a mapping of contig info and repeat counts.

    Expects contig names like 'contig_60' where 60 is the repeat count.
    Returns a dict keyed by both the full contig name (str) and the repeat
    count (int), each mapping to the number of mapped reads.

    Args:
        idxstats_output: Raw text output from ``samtools idxstats``.

    Returns:
        Dictionary mapping both contig name (str) and repeat count (int)
        to mapped read count (int).  The '*' unmapped line is excluded.
    """
    counts: dict = {}

    for line in idxstats_output.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue

        contig_name = parts[0]
        if contig_name == "*":
            continue

        mapped_reads = int(parts[2])

        # Store by full contig name
        counts[contig_name] = mapped_reads

        # Also store by repeat count extracted from name (e.g., "contig_60" -> 60)
        match = re.search(r"_(\d+)$", contig_name)
        if match:
            repeat_count = int(match.group(1))
            counts[repeat_count] = mapped_reads

    return counts


def detect_alleles(
    counts: dict,
    min_coverage: int = 10,
) -> dict:
    """Detect allele lengths from read count distribution.

    Finds the two contigs with the highest read counts (peaks) that are
    at least 3 repeat units apart.  If only one such peak exceeds
    ``min_coverage``, reports homozygous.

    Args:
        counts: Mapping produced by :func:`parse_idxstats`.  Both string
            and integer keys are accepted; only integer keys (repeat counts)
            are used for peak detection.
        min_coverage: Minimum mapped reads to consider a peak valid.

    Returns:
        Dictionary with keys ``allele_1``, ``allele_2``, and ``homozygous``.
        Each allele value is a dict with ``length`` (int) and ``reads`` (int).
        ``allele_1`` always has at least as many reads as ``allele_2``.

    Raises:
        ValueError: If no contig meets the minimum coverage threshold.
    """
    # Work only with integer (repeat-count) keys
    int_counts = {k: v for k, v in counts.items() if isinstance(k, int)}

    # Filter to contigs meeting minimum coverage
    passing = {k: v for k, v in int_counts.items() if v >= min_coverage}

    if not passing:
        max_observed = max(int_counts.values()) if int_counts else 0
        raise ValueError(
            f"No contig has >= {min_coverage} mapped reads (minimum coverage). "
            f"Max observed: {max_observed} reads."
        )

    # Sort by read count descending
    sorted_peaks = sorted(passing.items(), key=lambda x: x[1], reverse=True)

    # Primary allele: highest read count
    allele_1_length, allele_1_reads = sorted_peaks[0]

    # Find second peak: must be >= 3 repeat units away to avoid shoulder noise
    allele_2_length: int | None = None
    allele_2_reads = 0

    for length, reads in sorted_peaks[1:]:
        if abs(length - allele_1_length) >= 3:
            allele_2_length = length
            allele_2_reads = reads
            break

    if allele_2_length is None:
        # All other peaks are within shoulder distance -- homozygous
        return {
            "allele_1": {"length": allele_1_length, "reads": allele_1_reads},
            "allele_2": {"length": allele_1_length, "reads": 0},
            "homozygous": True,
        }

    return {
        "allele_1": {"length": allele_1_length, "reads": allele_1_reads},
        "allele_2": {"length": allele_2_length, "reads": allele_2_reads},
        "homozygous": False,
    }
