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


def trim_flanking(
    consensus_fasta: Path,
    flank_length: int,
    output_path: Path,
) -> Path:
    """Remove flanking sequences from a consensus FASTA, keeping only the VNTR.

    The ladder contigs have structure:
    ``[left_flank] [pre-repeats] [N * X] [after-repeats] [right_flank]``

    Flanking sequences are genomic sequence outside the VNTR and do not
    contain variants, so their positions are stable after ``bcftools
    consensus``.  Trimming them is safe because Clair3 only calls variants
    within the VNTR region where reads actually align.

    Args:
        consensus_fasta: Path to the full-contig consensus FASTA.
        flank_length: Number of flanking bp on each side (must match the
            value used during ladder generation, typically 500).
        output_path: Destination path for the trimmed FASTA.

    Returns:
        Path to the trimmed FASTA containing only the VNTR region.
    """
    lines = consensus_fasta.read_text().strip().splitlines()
    header = lines[0] if lines and lines[0].startswith(">") else ">consensus"
    sequence = "".join(line for line in lines if not line.startswith(">"))

    # Strip flanking from both ends
    vntr = sequence[flank_length:]
    if flank_length > 0:
        vntr = vntr[:-flank_length]

    output_path.write_text(f"{header}_vntr\n{vntr}\n")
    return output_path


def build_consensus_per_allele(
    reference_path: Path,
    vcf_paths: dict[str, Path],
    alleles: dict,
    output_dir: Path,
    flank_length: int = 500,
) -> dict[str, Path]:
    """Build a consensus FASTA for each detected allele and trim flanking.

    For each allele, the function:

    1. Extracts the matching single contig from the ladder reference with
       ``samtools faidx`` (using the peak contig name from allele detection).
    2. Indexes the single-contig FASTA.
    3. Calls :func:`build_consensus` to apply the allele's VCF variants.
    4. Trims flanking sequences to isolate the VNTR region.

    Args:
        reference_path: Path to the full ladder reference FASTA.
        vcf_paths: Mapping from allele key (e.g. ``"allele_1"``) to VCF path.
        alleles: Allele detection result from
            :func:`~open_pacmuci.alleles.detect_alleles`.
        output_dir: Base output directory for consensus files.
        flank_length: Flanking sequence length used in ladder generation.

    Returns:
        Dictionary mapping allele key to the trimmed consensus FASTA path
        (containing only the VNTR region, ready for repeat classification).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, Path] = {}

    for allele_key, vcf_path in vcf_paths.items():
        allele_info = alleles[allele_key]
        # Use contig_name from allele detection (canonical repeat count),
        # not total length.  Fall back for backwards compatibility.
        contig_name = allele_info.get("contig_name", f"contig_{allele_info['length']}")

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

        # Build full consensus (with flanking)
        full_consensus = output_dir / f"consensus_{allele_key}_full.fa"
        build_consensus(contig_fa, vcf_path, full_consensus)

        # Trim flanking to get VNTR-only sequence for classification
        trimmed = output_dir / f"consensus_{allele_key}.fa"
        trim_flanking(full_consensus, flank_length, trimmed)
        results[allele_key] = trimmed

    return results
