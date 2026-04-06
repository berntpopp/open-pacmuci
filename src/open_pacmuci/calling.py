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
    min_qual: float = 0.0,
    min_dp: int = 0,
) -> Path:
    """Normalize and filter a VCF file with bcftools.

    Runs ``bcftools norm -f <reference>`` followed by
    ``bcftools view -f PASS`` (with optional quality filters) and indexes
    the result.

    Args:
        vcf_path: Path to input VCF (may be gzipped).
        reference_path: Path to reference FASTA for left-normalisation.
        output_dir: Directory for output files.
        min_qual: Minimum QUAL score to keep a variant (0 = no filter).
        min_dp: Minimum INFO/DP to keep a variant (0 = no filter).

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

    # Build view command with PASS filter and optional quality filters
    view_cmd = [
        "bcftools",
        "view",
        "-f",
        "PASS",
    ]
    filters = []
    if min_qual > 0:
        filters.append(f"QUAL>={min_qual}")
    if min_dp > 0:
        filters.append(f"INFO/DP>={min_dp}")
    if filters:
        view_cmd.extend(["-i", " && ".join(filters)])

    view_cmd.extend(
        [
            "-o",
            str(filtered),
            "-O",
            "z",
            str(norm_vcf),
        ]
    )
    run_tool(view_cmd)
    run_tool(["bcftools", "index", str(filtered)])

    # Remove intermediate normalized VCF
    norm_vcf.unlink(missing_ok=True)

    return filtered


def parse_vcf_genotypes(vcf_path: Path) -> list[dict]:
    """Parse VCF to extract variant positions and genotypes.

    Returns list of dicts with keys: chrom, pos, ref, alt, genotype.
    """
    try:
        output = run_tool(
            [
                "bcftools",
                "query",
                "-f",
                "%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n",
                str(vcf_path),
            ]
        )
    except RuntimeError:
        return []

    variants: list[dict] = []
    for line in output.strip().splitlines():
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        variants.append(
            {
                "chrom": fields[0],
                "pos": int(fields[1]),
                "ref": fields[2],
                "alt": fields[3],
                "genotype": fields[4],
            }
        )
    return variants


def parse_vcf_variants(vcf_path: Path) -> list[dict]:
    """Parse VCF to extract variant positions and quality scores."""
    try:
        output = run_tool(["bcftools", "query", "-f", "%POS\\t%QUAL\\n", str(vcf_path)])
    except RuntimeError:
        return []
    variants = []
    for line in output.strip().splitlines():
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            try:
                variants.append({"pos": int(parts[0]), "qual": float(parts[1])})
            except ValueError:
                continue
    return variants


def disambiguate_same_length_alleles(
    bam_path: Path,
    reference_path: Path,
    alleles: dict,
    output_dir: Path,
    clair3_model: str = "",
    threads: int = 4,
    min_qual: float = 15.0,
    min_dp: int = 5,
) -> dict:
    """Disambiguate same-length alleles using Clair3 genotype calls.

    Runs Clair3 on all reads, checks for heterozygous (0/1) variants.
    If found, creates two VCFs: one with only hom-alt variants (WT allele),
    one with all variants (mutant allele).

    Returns dict mapping allele key to filtered VCF path, plus a
    ``"homozygous"`` boolean key indicating whether the alleles are identical.
    """
    allele_info = alleles["allele_1"]
    contig_name = allele_info["contig_name"]
    cluster_contigs = allele_info["cluster_contigs"]
    merged_dir = output_dir / "merged"

    # Remap all cluster reads to peak contig
    merged_bam = _extract_and_remap_reads(
        bam_path,
        cluster_contigs,
        contig_name,
        reference_path,
        merged_dir,
        threads,
    )

    # Run Clair3
    contig_ref = merged_dir / f"{contig_name}.fa"
    clair3_dir = merged_dir / "clair3"
    raw_vcf = run_clair3(
        merged_bam,
        contig_ref,
        clair3_dir,
        model_path=clair3_model,
        threads=threads,
    )
    filtered_vcf = filter_vcf(raw_vcf, contig_ref, merged_dir, min_qual=min_qual, min_dp=min_dp)

    # Check genotypes
    variants = parse_vcf_genotypes(filtered_vcf)
    het_variants = [
        v
        for v in variants
        if any(sep in v["genotype"] for sep in ("/", "|"))
        and v["genotype"] not in ("0/0", "1/1", "./.", "0|0", "1|1", ".|.")
    ]

    results: dict = {}

    if not het_variants:
        # Truly homozygous -- no het variants found
        results["allele_1"] = filtered_vcf
        results["homozygous"] = True
        return results

    # Compound heterozygous -- split into two VCFs

    # allele_1 (WT): exclude het variants, keep only hom-alt
    wt_vcf = output_dir / "allele_1" / "variants.vcf.gz"
    wt_vcf.parent.mkdir(parents=True, exist_ok=True)
    run_tool(
        [
            "bcftools",
            "view",
            "-i",
            'GT="1/1" || GT="1|1"',
            "-o",
            str(wt_vcf),
            "-O",
            "z",
            str(filtered_vcf),
        ]
    )
    run_tool(["bcftools", "index", str(wt_vcf)])
    results["allele_1"] = wt_vcf

    # allele_2 (mutant): include all variants
    mut_vcf = output_dir / "allele_2" / "variants.vcf.gz"
    mut_vcf.parent.mkdir(parents=True, exist_ok=True)
    run_tool(
        [
            "bcftools",
            "view",
            "-o",
            str(mut_vcf),
            "-O",
            "z",
            str(filtered_vcf),
        ]
    )
    run_tool(["bcftools", "index", str(mut_vcf)])
    results["allele_2"] = mut_vcf

    # Note: bcftools consensus applies the ALT allele for het (0/1) calls
    # by default, so the allele_2 VCF with het genotypes will produce the
    # correct mutant consensus without needing to force hom-alt genotypes.
    results["homozygous"] = False
    return results


def call_variants_per_allele(
    bam_path: Path,
    reference_path: Path,
    alleles: dict,
    output_dir: Path,
    clair3_model: str = "",
    threads: int = 4,
    min_qual: float = 15.0,
    min_dp: int = 5,
) -> dict[str, Path]:
    """Run variant calling for each detected allele.

    For each allele key in *alleles* (``"allele_1"`` and optionally
    ``"allele_2"``), this function:

    1. Extracts the reads mapped to the allele's best-matching contig.
    2. Calls variants with Clair3.
    3. Filters the VCF with :func:`filter_vcf`.

    For same-length alleles, delegates to :func:`disambiguate_same_length_alleles`.

    Args:
        bam_path: Path to the full mapping BAM.
        reference_path: Path to the ladder reference FASTA.
        alleles: Allele detection result from :func:`~open_pacmuci.alleles.detect_alleles`.
        output_dir: Base output directory.
        clair3_model: Path to Clair3 model directory (optional).
        threads: Number of threads (default 4).
        min_qual: Minimum QUAL score for VCF filtering (default 15.0).
        min_dp: Minimum INFO/DP for VCF filtering (default 5).

    Returns:
        Dictionary mapping allele key (``"allele_1"`` / ``"allele_2"``) to
        the filtered VCF path.
    """
    if alleles.get("same_length"):
        disambig = disambiguate_same_length_alleles(
            bam_path,
            reference_path,
            alleles,
            output_dir,
            clair3_model,
            threads,
            min_qual,
            min_dp,
        )
        alleles["homozygous"] = disambig.pop("homozygous")
        return disambig

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
        filtered = filter_vcf(vcf, contig_ref, allele_dir, min_qual=min_qual, min_dp=min_dp)
        results[allele_key] = filtered

    return results
