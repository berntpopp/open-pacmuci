"""Variant calling with Clair3 and VCF processing with bcftools."""

from __future__ import annotations

from pathlib import Path

from open_pacmuci.tools import run_tool


def extract_allele_reads(
    bam_path: Path,
    contig_names: str | list[str],
    output_dir: Path,
) -> Path:
    """Extract reads mapped to one or more contigs from a BAM file.

    Uses ``samtools view -b`` to subset the BAM to the specified contig(s)
    and indexes the result.

    Args:
        bam_path: Path to the full mapping BAM.
        contig_names: Single contig name or list of contig names to extract.
        output_dir: Directory for output files.

    Returns:
        Path to the extracted, indexed BAM file.
    """
    if isinstance(contig_names, str):
        contig_names = [contig_names]

    output_dir.mkdir(parents=True, exist_ok=True)
    out_bam = output_dir / "allele_reads.bam"

    run_tool(
        [
            "samtools",
            "view",
            "-b",
            "-o",
            str(out_bam),
            str(bam_path),
            *contig_names,
        ]
    )
    run_tool(["samtools", "index", str(out_bam)])

    return out_bam


def _extract_and_remap_reads(
    bam_path: Path,
    cluster_contigs: list[str],
    peak_contig: str,
    reference_path: Path,
    output_dir: Path,
    threads: int = 4,
) -> Path:
    """Extract reads from cluster contigs, convert to FASTQ, remap to peak contig.

    Reads from the ladder mapping are spread across multiple contigs in a
    cluster (e.g. contig_48 through contig_54).  Clair3 needs all reads
    aligned to a *single* reference contig to call variants.  This function
    extracts the cluster reads, converts to FASTQ, extracts the single peak
    contig as a mini-reference, and remaps with minimap2.

    Args:
        bam_path: Full ladder mapping BAM.
        cluster_contigs: All contig names in the allele's cluster.
        peak_contig: The single peak contig name to remap against.
        reference_path: Full ladder reference FASTA (for extracting the contig).
        output_dir: Working directory for intermediate files.
        threads: Thread count for minimap2/samtools.

    Returns:
        Path to the remapped, sorted, indexed BAM.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Extract reads from all cluster contigs
    cluster_bam = extract_allele_reads(bam_path, cluster_contigs, output_dir)

    # 2. Convert to FASTQ
    fastq_path = output_dir / "cluster_reads.fq"
    stdout = run_tool(["samtools", "fastq", str(cluster_bam)])
    fastq_path.write_text(stdout)

    # 3. Extract peak contig as mini-reference
    contig_ref = output_dir / f"{peak_contig}.fa"
    stdout = run_tool(
        [
            "samtools",
            "faidx",
            str(reference_path),
            peak_contig,
        ]
    )
    contig_ref.write_text(stdout)
    run_tool(["samtools", "faidx", str(contig_ref)])

    # 4. Remap to peak contig
    sam_path = output_dir / "remapped.sam"
    sam_output = run_tool(
        [
            "minimap2",
            "-a",
            "-x",
            "map-hifi",
            "-t",
            str(threads),
            str(contig_ref),
            str(fastq_path),
        ]
    )
    sam_path.write_text(sam_output)

    # 5. Sort and index
    remapped_bam = output_dir / "allele_reads.bam"
    run_tool(
        [
            "samtools",
            "sort",
            "-@",
            str(threads),
            "-o",
            str(remapped_bam),
            str(sam_path),
        ]
    )
    run_tool(["samtools", "index", str(remapped_bam)])

    # Clean up intermediates
    sam_path.unlink(missing_ok=True)
    fastq_path.unlink(missing_ok=True)

    return remapped_bam


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
        # Our contigs are named contig_N, not chr1..22/X/Y, so Clair3
        # must be told to process all contigs.
        "--include_all_ctgs",
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

        allele_info = alleles[allele_key]
        # Use the peak contig name from allele detection (contig_N where N is
        # canonical X repeat count, NOT total length).  Fall back to computing
        # from length for backwards compatibility with older alleles.json.
        contig_name = allele_info.get("contig_name", f"contig_{allele_info['length']}")
        cluster_contigs = allele_info.get("cluster_contigs", [contig_name])
        allele_dir = output_dir / allele_key

        # Extract reads from ALL contigs in the cluster, then remap to
        # the peak contig.  Reads spread across multiple ladder contigs
        # due to length variation; Clair3 needs them all aligned to a
        # single reference contig to call variants effectively.
        allele_bam = _extract_and_remap_reads(
            bam_path,
            cluster_contigs,
            contig_name,
            reference_path,
            allele_dir,
            threads,
        )

        # Run Clair3 against the single-contig reference (created during
        # remapping).  Using the single contig rather than the full ladder
        # avoids Clair3 scanning 150 empty contigs.
        contig_ref = allele_dir / f"{contig_name}.fa"
        clair3_dir = allele_dir / "clair3"
        vcf = run_clair3(
            allele_bam,
            contig_ref,
            clair3_dir,
            model_path=clair3_model,
            threads=threads,
        )

        # Filter VCF
        filtered = filter_vcf(vcf, contig_ref, allele_dir)
        results[allele_key] = filtered

    return results
