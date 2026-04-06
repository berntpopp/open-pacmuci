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
        peak_contig (int, contig with most reads),
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

    # Compute weighted center, peak contig, and total reads for each cluster
    result: list[dict] = []
    for cluster in clusters:
        total_reads = sum(reads for _, reads in cluster)
        weighted_center = sum(pos * reads for pos, reads in cluster) / total_reads
        peak_contig = max(cluster, key=lambda x: x[1])[0]
        result.append(
            {
                "center": round(weighted_center),
                "total_reads": total_reads,
                "peak_contig": peak_contig,
                "contigs": cluster,
            }
        )

    result.sort(key=lambda x: x["total_reads"], reverse=True)
    return result


def _build_allele_info(cluster: dict) -> dict:
    """Build allele info dict from a cluster."""
    canonical = cluster["center"]
    peak = cluster["peak_contig"]
    return {
        "length": canonical + PRE_AFTER_REPEAT_COUNT,
        "reads": cluster["total_reads"],
        "canonical_repeats": canonical,
        "contig_name": f"contig_{peak}",
        "cluster_contigs": [f"contig_{c}" for c, _ in cluster["contigs"]],
    }


def detect_alleles(
    counts: dict[int, int],
    min_coverage: int = 10,
) -> dict:
    """Detect allele lengths from read count distribution across ladder contigs.

    Finds two peak clusters in the read distribution. Each cluster represents
    one allele. Reports the total allele length (canonical repeats + 9 fixed
    pre/after repeats) and the contig names needed for downstream processing.

    Args:
        counts: Canonical repeat count -> mapped reads mapping from
            :func:`parse_idxstats`.
        min_coverage: Minimum mapped reads to include a contig.

    Returns:
        Dictionary with keys ``allele_1``, ``allele_2``, and ``homozygous``.
        Each allele has:

        - ``length`` (int): total repeat units including pre/after
        - ``reads`` (int): total mapped reads across the cluster
        - ``canonical_repeats`` (int): number of canonical X repeats
        - ``contig_name`` (str): name of the peak contig (e.g. ``"contig_51"``)
        - ``cluster_contigs`` (list[str]): all contig names in the cluster

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

    allele_1 = _build_allele_info(clusters[0])

    if len(clusters) < 2:
        return {
            "allele_1": allele_1,
            "allele_2": {**allele_1, "reads": 0},
            "homozygous": True,
        }

    allele_2 = _build_allele_info(clusters[1])

    if allele_1["length"] == allele_2["length"]:
        allele_1["reads"] += allele_2["reads"]
        return {
            "allele_1": allele_1,
            "allele_2": {**allele_1, "reads": 0},
            "homozygous": True,
        }

    return {
        "allele_1": allele_1,
        "allele_2": allele_2,
        "homozygous": False,
    }
