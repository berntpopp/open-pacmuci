"""Read mapping with minimap2 and samtools."""

from __future__ import annotations

from pathlib import Path

from open_pacmuci.tools import run_tool


def bam_to_fastq(bam_path: Path, output_dir: Path) -> Path:
    """Convert BAM to FASTQ using samtools fastq.

    Captures stdout from ``samtools fastq`` and writes it to
    ``extracted_reads.fq`` in *output_dir*.

    Args:
        bam_path: Path to input BAM file.
        output_dir: Directory for output FASTQ.

    Returns:
        Path to the output FASTQ file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    fastq_path = output_dir / "extracted_reads.fq"
    stdout = run_tool(["samtools", "fastq", str(bam_path)])
    fastq_path.write_text(stdout)
    return fastq_path


def map_reads(
    input_path: Path,
    reference_path: Path,
    output_dir: Path,
    threads: int = 4,
) -> Path:
    """Map reads to reference using minimap2 and sort/index with samtools.

    Pipeline: ``minimap2 -a -x map-hifi`` → ``samtools sort`` →
    ``samtools index``.  If *input_path* is a BAM file it is converted to
    FASTQ first with :func:`bam_to_fastq`.

    Args:
        input_path: Path to input FASTQ or BAM file.
        reference_path: Path to reference FASTA.
        output_dir: Directory for output files.
        threads: Number of threads for minimap2/samtools (default 4).

    Returns:
        Path to the sorted, indexed BAM file (``mapping.bam``).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # If input is BAM, convert to FASTQ first
    actual_input = input_path
    if input_path.suffix.lower() == ".bam":
        actual_input = bam_to_fastq(input_path, output_dir)

    sam_path = output_dir / "mapping.sam"
    bam_path = output_dir / "mapping.bam"

    # minimap2 alignment -- stdout is SAM
    sam_output = run_tool(
        [
            "minimap2",
            "-a",
            "-x",
            "map-hifi",
            "-t",
            str(threads),
            str(reference_path),
            str(actual_input),
        ]
    )
    sam_path.write_text(sam_output)

    # samtools sort → BAM
    run_tool(
        [
            "samtools",
            "sort",
            "-@",
            str(threads),
            "-o",
            str(bam_path),
            str(sam_path),
        ]
    )

    # samtools index
    run_tool(["samtools", "index", str(bam_path)])

    # Remove intermediate SAM
    sam_path.unlink(missing_ok=True)

    return bam_path


def get_idxstats(bam_path: Path) -> str:
    """Run ``samtools idxstats`` and return the raw text output.

    Args:
        bam_path: Path to an indexed BAM file.

    Returns:
        Raw idxstats text output (tab-separated).
    """
    return run_tool(["samtools", "idxstats", str(bam_path)])
