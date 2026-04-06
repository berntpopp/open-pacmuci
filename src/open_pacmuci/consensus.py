"""Per-allele consensus sequence generation with bcftools."""

from __future__ import annotations

from pathlib import Path

from open_pacmuci.tools import run_tool


def build_consensus(
    reference_path: Path,
    vcf_path: Path,
    output_path: Path,
) -> Path:
    """Generate a consensus FASTA by applying VCF variants to a reference.

    Runs ``bcftools consensus -f <reference> <vcf>`` and writes the stdout
    to *output_path*.

    Args:
        reference_path: Path to the reference FASTA (typically the matching
            ladder contig extracted with ``samtools faidx``).
        vcf_path: Path to the filtered VCF (may be gzipped).
        output_path: Destination path for the consensus FASTA.

    Returns:
        Path to the written consensus FASTA.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    stdout = run_tool(
        [
            "bcftools",
            "consensus",
            "-f",
            str(reference_path),
            str(vcf_path),
        ]
    )

    output_path.write_text(stdout)
    return output_path


def build_consensus_per_allele(
    reference_path: Path,
    vcf_paths: dict[str, Path],
    alleles: dict,
    output_dir: Path,
) -> dict[str, Path]:
    """Build a consensus FASTA for each detected allele.

    For each allele, the function:

    1. Extracts the matching single contig from the ladder reference with
       ``samtools faidx``.
    2. Indexes the single-contig FASTA.
    3. Calls :func:`build_consensus` to apply the allele's VCF variants.

    Args:
        reference_path: Path to the full ladder reference FASTA.
        vcf_paths: Mapping from allele key (e.g. ``"allele_1"``) to VCF path.
        alleles: Allele detection result from
            :func:`~open_pacmuci.alleles.detect_alleles`.
        output_dir: Base output directory for consensus files.

    Returns:
        Dictionary mapping allele key to the consensus FASTA path.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, Path] = {}

    for allele_key, vcf_path in vcf_paths.items():
        length = alleles[allele_key]["length"]
        contig_name = f"contig_{length}"

        # Extract the single matching contig from the full ladder reference
        contig_fa = output_dir / f"ref_{contig_name}.fa"
        stdout = run_tool(
            [
                "samtools",
                "faidx",
                str(reference_path),
                contig_name,
            ]
        )
        contig_fa.write_text(stdout)

        # Index the single-contig reference so bcftools consensus can use it
        run_tool(["samtools", "faidx", str(contig_fa)])

        # Build consensus
        consensus_path = output_dir / f"consensus_{allele_key}.fa"
        build_consensus(contig_fa, vcf_path, consensus_path)
        results[allele_key] = consensus_path

    return results
