"""Allele length detection from samtools idxstats output."""

from __future__ import annotations

import re

# Number of fixed repeat units in the ladder reference (pre-repeats 1-5 + after-repeats 6-9).
# Each contig_N has N canonical X repeats plus these fixed repeats,
# so total allele length = N + PRE_AFTER_REPEAT_COUNT.
PRE_AFTER_REPEAT_COUNT = 9


def parse_idxstats(idxstats_output: str) -> dict[int, int]:
    """Parse samtools idxstats output into repeat_count -> read_count mapping.

    Expects contig names like 'contig_60' where 60 is the number of
    canonical X repeats in that contig.

    Args:
        idxstats_output: Raw text output from ``samtools idxstats``.

    Returns:
        Dictionary mapping canonical repeat count (int) to mapped read
        count (int).  The '*' unmapped line is excluded.
    """
    counts: dict[int, int] = {}

    for line in idxstats_output.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue

        contig_name = parts[0]
        if contig_name == "*":
            continue

        mapped_reads = int(parts[2])

        match = re.search(r"_(\d+)$", contig_name)
        if match:
            repeat_count = int(match.group(1))
            counts[repeat_count] = mapped_reads

    return counts


def _find_clusters(
    counts: dict[int, int],
    min_coverage: int,
    min_gap: int = 5,
) -> list[dict]:
    """Identify read-count clusters in the contig distribution.

    Groups contigs that are within ``min_gap`` of each other into clusters,
    then computes the weighted center and total reads for each.

    Args:
        counts: Canonical repeat count -> mapped reads mapping.
        min_coverage: Minimum reads for a contig to be included.
        min_gap: Minimum gap between contigs to start a new cluster.

    Returns:
        List of cluster dicts sorted by total_reads descending.
        Each dict has keys: center (int), total_reads (int),
        contigs (list of (repeat_count, reads) tuples).
    """
    passing = sorted(
        [(k, v) for k, v in counts.items() if v >= min_coverage],
        key=lambda x: x[0],
    )

    if not passing:
        return []

    # Group into clusters by proximity
    clusters: list[list[tuple[int, int]]] = []
    current_cluster: list[tuple[int, int]] = [passing[0]]

    for i in range(1, len(passing)):
        if passing[i][0] - passing[i - 1][0] >= min_gap:
            clusters.append(current_cluster)
            current_cluster = [passing[i]]
        else:
            current_cluster.append(passing[i])
    clusters.append(current_cluster)

    # Compute weighted center and total reads for each cluster
    result: list[dict] = []
    for cluster in clusters:
        total_reads = sum(reads for _, reads in cluster)
        weighted_center = sum(pos * reads for pos, reads in cluster) / total_reads
        result.append(
            {
                "center": round(weighted_center),
                "total_reads": total_reads,
                "contigs": cluster,
            }
        )

    result.sort(key=lambda x: x["total_reads"], reverse=True)
    return result


def detect_alleles(
    counts: dict[int, int],
    min_coverage: int = 10,
) -> dict:
    """Detect allele lengths from read count distribution across ladder contigs.

    Finds two peak clusters in the read distribution. Each cluster represents
    one allele. Reports the total allele length (canonical repeats + 9 fixed
    pre/after repeats).

    Args:
        counts: Canonical repeat count -> mapped reads mapping from
            :func:`parse_idxstats`.
        min_coverage: Minimum mapped reads to include a contig.

    Returns:
        Dictionary with keys ``allele_1``, ``allele_2``, and ``homozygous``.
        Each allele has ``length`` (total repeat units) and ``reads`` (int).
        ``allele_1`` always has at least as many reads as ``allele_2``.

    Raises:
        ValueError: If no contig meets the minimum coverage threshold.
    """
    # Handle mixed-key dicts (legacy compat): only use integer keys
    int_counts = {k: v for k, v in counts.items() if isinstance(k, int)}

    clusters = _find_clusters(int_counts, min_coverage)

    if not clusters:
        max_observed = max(int_counts.values()) if int_counts else 0
        raise ValueError(
            f"No contig has >= {min_coverage} mapped reads (minimum coverage). "
            f"Max observed: {max_observed} reads."
        )

    # Primary allele: cluster with most total reads
    allele_1_canonical = clusters[0]["center"]
    allele_1_reads = clusters[0]["total_reads"]
    allele_1_length = allele_1_canonical + PRE_AFTER_REPEAT_COUNT

    if len(clusters) < 2:
        return {
            "allele_1": {"length": allele_1_length, "reads": allele_1_reads},
            "allele_2": {"length": allele_1_length, "reads": 0},
            "homozygous": True,
        }

    # Second allele: second-largest cluster
    allele_2_canonical = clusters[1]["center"]
    allele_2_reads = clusters[1]["total_reads"]
    allele_2_length = allele_2_canonical + PRE_AFTER_REPEAT_COUNT

    # Homozygous if both clusters center on the same length
    if allele_1_length == allele_2_length:
        return {
            "allele_1": {"length": allele_1_length, "reads": allele_1_reads + allele_2_reads},
            "allele_2": {"length": allele_1_length, "reads": 0},
            "homozygous": True,
        }

    return {
        "allele_1": {"length": allele_1_length, "reads": allele_1_reads},
        "allele_2": {"length": allele_2_length, "reads": allele_2_reads},
        "homozygous": False,
    }
