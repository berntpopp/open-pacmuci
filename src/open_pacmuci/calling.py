"""Variant calling with Clair3 and VCF processing with bcftools."""

from __future__ import annotations

from pathlib import Path

from open_pacmuci.tools import run_tool


def extract_allele_reads(
    bam_path: Path,
    contig_name: str,
    output_dir: Path,
) -> Path:
    """Extract reads mapped to a specific contig from a BAM file.

    Uses ``samtools view -b`` to subset the BAM to a single contig and
    indexes the result.

    Args:
        bam_path: Path to the full mapping BAM.
        contig_name: Name of the contig to extract (e.g. ``"contig_60"``).
        output_dir: Directory for output files.

    Returns:
        Path to the extracted, indexed BAM file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    out_bam = output_dir / f"reads_{contig_name}.bam"

    run_tool(
        [
            "samtools",
            "view",
            "-b",
            "-o",
            str(out_bam),
            str(bam_path),
            contig_name,
        ]
    )
    run_tool(["samtools", "index", str(out_bam)])

    return out_bam


def run_clair3(
    bam_path: Path,
    reference_path: Path,
    output_dir: Path,
    model_path: str = "",
    platform: str = "hifi",
    threads: int = 4,
) -> Path:
    """Run the Clair3 variant caller on a BAM file.

    Args:
        bam_path: Path to input BAM (reads for one allele).
        reference_path: Path to the reference FASTA.
        output_dir: Directory for Clair3 output.
        model_path: Path to Clair3 model directory (optional).
        platform: Sequencing platform (default ``"hifi"``).
        threads: Number of threads (default 4).

    Returns:
        Path to the Clair3 output VCF (``merge_output.vcf.gz``).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "run_clair3.sh",
        f"--bam_fn={bam_path}",
        f"--ref_fn={reference_path}",
        f"--output={output_dir}",
        f"--threads={threads}",
        f"--platform={platform}",
        "--sample_name=sample",
    ]
    if model_path:
        cmd.append(f"--model_path={model_path}")

    run_tool(cmd)

    return output_dir / "merge_output.vcf.gz"


def filter_vcf(
    vcf_path: Path,
    reference_path: Path,
    output_dir: Path,
) -> Path:
    """Normalize and filter a VCF file with bcftools.

    Runs ``bcftools norm -f <reference>`` followed by
    ``bcftools view -f PASS`` and indexes the result.

    Args:
        vcf_path: Path to input VCF (may be gzipped).
        reference_path: Path to reference FASTA for left-normalisation.
        output_dir: Directory for output files.

    Returns:
        Path to the filtered, indexed VCF (``variants.vcf.gz``).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    norm_vcf = output_dir / "normalized.vcf.gz"
    filtered = output_dir / "variants.vcf.gz"

    # Left-normalise indels against the reference
    run_tool(
        [
            "bcftools",
            "norm",
            "-f",
            str(reference_path),
            "-o",
            str(norm_vcf),
            "-O",
            "z",
            str(vcf_path),
        ]
    )

    # Keep only PASS variants
    run_tool(
        [
            "bcftools",
            "view",
            "-f",
            "PASS",
            "-o",
            str(filtered),
            "-O",
            "z",
            str(norm_vcf),
        ]
    )
    run_tool(["bcftools", "index", str(filtered)])

    # Remove intermediate normalized VCF
    norm_vcf.unlink(missing_ok=True)

    return filtered


def call_variants_per_allele(
    bam_path: Path,
    reference_path: Path,
    alleles: dict,
    output_dir: Path,
    clair3_model: str = "",
    threads: int = 4,
) -> dict[str, Path]:
    """Run variant calling for each detected allele.

    For each allele key in *alleles* (``"allele_1"`` and optionally
    ``"allele_2"``), this function:

    1. Extracts the reads mapped to the allele's best-matching contig.
    2. Calls variants with Clair3.
    3. Filters the VCF with :func:`filter_vcf`.

    Args:
        bam_path: Path to the full mapping BAM.
        reference_path: Path to the ladder reference FASTA.
        alleles: Allele detection result from :func:`~open_pacmuci.alleles.detect_alleles`.
        output_dir: Base output directory.
        clair3_model: Path to Clair3 model directory (optional).
        threads: Number of threads (default 4).

    Returns:
        Dictionary mapping allele key (``"allele_1"`` / ``"allele_2"``) to
        the filtered VCF path.
    """
    results: dict[str, Path] = {}

    for allele_key in ("allele_1", "allele_2"):
        if alleles.get("homozygous") and allele_key == "allele_2":
            continue

        length = alleles[allele_key]["length"]
        contig_name = f"contig_{length}"
        allele_dir = output_dir / allele_key

        # Extract reads for this allele's contig
        allele_bam = extract_allele_reads(bam_path, contig_name, allele_dir)

        # Run Clair3
        clair3_dir = allele_dir / "clair3"
        vcf = run_clair3(
            allele_bam,
            reference_path,
            clair3_dir,
            model_path=clair3_model,
            threads=threads,
        )

        # Filter VCF
        filtered = filter_vcf(vcf, reference_path, allele_dir)
        results[allele_key] = filtered

    return results
