"""Smoke tests for the CLI interface."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from open_pacmuci.cli import main


class TestCli:
    """CLI smoke tests."""

    def test_help(self):
        """CLI --help exits with 0."""
        runner = CliRunner()
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "open-pacmuci" in result.output.lower() or "MUC1" in result.output

    def test_version(self):
        """CLI --version shows version."""
        runner = CliRunner()
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0
        assert "0.2.0" in result.output

    def test_ladder_help(self):
        """ladder subcommand has help."""
        runner = CliRunner()
        result = runner.invoke(main, ["ladder", "--help"])
        assert result.exit_code == 0
        assert "ladder" in result.output.lower() or "reference" in result.output.lower()

    def test_map_help(self):
        """map subcommand has help."""
        runner = CliRunner()
        result = runner.invoke(main, ["map", "--help"])
        assert result.exit_code == 0

    def test_alleles_help(self):
        """alleles subcommand has help."""
        runner = CliRunner()
        result = runner.invoke(main, ["alleles", "--help"])
        assert result.exit_code == 0

    def test_call_help(self):
        """call subcommand has help."""
        runner = CliRunner()
        result = runner.invoke(main, ["call", "--help"])
        assert result.exit_code == 0

    def test_consensus_help(self):
        """consensus subcommand has help."""
        runner = CliRunner()
        result = runner.invoke(main, ["consensus", "--help"])
        assert result.exit_code == 0

    def test_classify_help(self):
        """classify subcommand has help."""
        runner = CliRunner()
        result = runner.invoke(main, ["classify", "--help"])
        assert result.exit_code == 0

    def test_run_help(self):
        """run subcommand has help."""
        runner = CliRunner()
        result = runner.invoke(main, ["run", "--help"])
        assert result.exit_code == 0


class TestLadderSubcommand:
    """Tests for the ladder subcommand (no external tools required)."""

    def test_ladder_generates_fasta(self, tmp_path):
        """ladder command writes a FASTA file to the given output path."""
        runner = CliRunner()
        out_file = str(tmp_path / "ladder.fa")
        result = runner.invoke(
            main,
            ["ladder", "--output", out_file, "--min-units", "1", "--max-units", "5"],
        )
        assert result.exit_code == 0, result.output
        assert Path(out_file).exists()

    def test_ladder_output_message(self, tmp_path):
        """ladder command echoes the output path."""
        runner = CliRunner()
        out_file = str(tmp_path / "ladder.fa")
        result = runner.invoke(
            main,
            ["ladder", "--output", out_file, "--min-units", "1", "--max-units", "3"],
        )
        assert result.exit_code == 0, result.output
        assert "Ladder written to" in result.output

    def test_ladder_contig_count_in_message(self, tmp_path):
        """ladder command reports the correct number of contigs."""
        runner = CliRunner()
        out_file = str(tmp_path / "ladder.fa")
        result = runner.invoke(
            main,
            ["ladder", "--output", out_file, "--min-units", "1", "--max-units", "10"],
        )
        assert result.exit_code == 0, result.output
        # 10 - 1 + 1 = 10 contigs
        assert "10 contigs" in result.output

    def test_ladder_fasta_contains_contigs(self, tmp_path):
        """Generated FASTA contains sequence headers."""
        runner = CliRunner()
        out_file = str(tmp_path / "ladder.fa")
        result = runner.invoke(
            main,
            ["ladder", "--output", out_file, "--min-units", "1", "--max-units", "3"],
        )
        assert result.exit_code == 0, result.output
        content = Path(out_file).read_text()
        assert ">" in content

    def test_ladder_flank_length_option(self, tmp_path):
        """ladder command accepts --flank-length option without error."""
        runner = CliRunner()
        out_file = str(tmp_path / "ladder.fa")
        result = runner.invoke(
            main,
            [
                "ladder",
                "--output",
                out_file,
                "--min-units",
                "1",
                "--max-units",
                "2",
                "--flank-length",
                "100",
            ],
        )
        assert result.exit_code == 0, result.output


class TestAllelesSubcommand:
    """Tests for the alleles subcommand with mocked external tools."""

    def test_alleles_writes_json(self, tmp_path):
        """alleles subcommand writes alleles.json to the output directory."""
        # Mock both mapping.run_tool (for idxstats) and alleles.run_tool
        # (for refine_peak_contig's samtools view call)
        with (
            patch("open_pacmuci.mapping.run_tool") as mock_mapping_run,
            patch("open_pacmuci.alleles.run_tool") as mock_alleles_run,
        ):
            mock_mapping_run.return_value = (
                "contig_60\t4120\t200\t0\ncontig_80\t5320\t150\t0\n*\t0\t0\t50\n"
            )
            mock_alleles_run.return_value = ""
            runner = CliRunner()
            bam = tmp_path / "mapping.bam"
            bam.touch()
            out_dir = str(tmp_path / "out")

            result = runner.invoke(
                main,
                [
                    "alleles",
                    "--input",
                    str(bam),
                    "--output-dir",
                    out_dir,
                    "--min-coverage",
                    "10",
                ],
            )

        assert result.exit_code == 0, result.output
        alleles_file = Path(out_dir) / "alleles.json"
        assert alleles_file.exists()
        data = json.loads(alleles_file.read_text())
        assert "allele_1" in data
        assert "allele_2" in data

    def test_alleles_echoes_result(self, tmp_path):
        """alleles subcommand prints the detected alleles."""
        with (
            patch("open_pacmuci.mapping.run_tool") as mock_mapping_run,
            patch("open_pacmuci.alleles.run_tool") as mock_alleles_run,
        ):
            mock_mapping_run.return_value = "contig_60\t4120\t200\t0\n*\t0\t0\t50\n"
            mock_alleles_run.return_value = ""
            runner = CliRunner()
            bam = tmp_path / "mapping.bam"
            bam.touch()

            result = runner.invoke(
                main,
                ["alleles", "--input", str(bam), "--output-dir", str(tmp_path)],
            )

        assert result.exit_code == 0, result.output
        assert "Alleles" in result.output


class TestClassifySubcommand:
    """Tests for the classify subcommand (pure Python, no external tools)."""

    def test_classify_writes_outputs(self, tmp_path):
        """classify subcommand writes repeats.json and repeats.txt."""
        # Build a minimal FASTA with known repeat sequence
        # Use a short valid sequence (will likely produce unknown units but not crash)
        sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        fasta = tmp_path / "consensus.fa"
        fasta.write_text(f">allele_1_vntr\n{sequence}\n")
        out_dir = str(tmp_path / "classify_out")

        runner = CliRunner()
        result = runner.invoke(
            main,
            ["classify", "--input", str(fasta), "--output-dir", out_dir],
        )

        assert result.exit_code == 0, result.output
        assert (Path(out_dir) / "repeats.json").exists()
        assert (Path(out_dir) / "repeats.txt").exists()

    def test_classify_echoes_structure(self, tmp_path):
        """classify subcommand prints the structure string."""
        sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        fasta = tmp_path / "consensus.fa"
        fasta.write_text(f">allele_1_vntr\n{sequence}\n")

        runner = CliRunner()
        result = runner.invoke(
            main,
            ["classify", "--input", str(fasta), "--output-dir", str(tmp_path)],
        )

        assert result.exit_code == 0, result.output
        assert "Structure:" in result.output


class TestMapSubcommand:
    """Tests for the map subcommand with mocked tools.

    The CLI uses lazy imports inside the command body, so we mock at the
    source module (open_pacmuci.mapping / open_pacmuci.tools) rather than
    trying to patch names on open_pacmuci.cli.
    """

    def test_map_calls_map_reads(self, tmp_path):
        """map subcommand calls map_reads and echoes the output path."""
        with (
            patch("open_pacmuci.tools.check_tools", return_value=True),
            patch("open_pacmuci.mapping._run_mapping_pipeline"),
            patch("open_pacmuci.mapping.run_tool", return_value=""),
        ):
            fastq = tmp_path / "reads.fq"
            fastq.touch()
            ref = tmp_path / "ref.fa"
            ref.touch()

            runner = CliRunner()
            result = runner.invoke(
                main,
                [
                    "map",
                    "--input",
                    str(fastq),
                    "--reference",
                    str(ref),
                    "--output-dir",
                    str(tmp_path),
                    "--threads",
                    "2",
                ],
            )

        assert result.exit_code == 0, result.output
        assert "Mapping written to" in result.output

    def test_map_echoes_bam_path(self, tmp_path):
        """map subcommand echoes the BAM path on success."""
        with (
            patch("open_pacmuci.tools.check_tools", return_value=True),
            patch("open_pacmuci.mapping._run_mapping_pipeline"),
            patch("open_pacmuci.mapping.run_tool", return_value=""),
        ):
            fastq = tmp_path / "reads.fq"
            fastq.touch()
            ref = tmp_path / "ref.fa"
            ref.touch()

            runner = CliRunner()
            result = runner.invoke(
                main,
                ["map", "--input", str(fastq), "--reference", str(ref)],
            )

        assert result.exit_code == 0, result.output
        assert "Mapping written to" in result.output


class TestCallSubcommand:
    """Tests for the call subcommand with mocked tools."""

    def test_call_invokes_variant_calling(self, tmp_path):
        """call subcommand runs without error when all tools are mocked."""
        with (
            patch("open_pacmuci.tools.check_tools", return_value=True),
            patch("open_pacmuci.calling.run_tool", return_value=""),
        ):
            bam = tmp_path / "mapping.bam"
            bam.touch()
            ref = tmp_path / "ref.fa"
            ref.touch()
            alleles_data = {
                "homozygous": True,
                "allele_1": {
                    "length": 60,
                    "reads": 200,
                    "canonical_repeats": 51,
                    "contig_name": "contig_51",
                    "cluster_contigs": ["contig_51"],
                },
                "allele_2": {
                    "length": 60,
                    "reads": 0,
                    "canonical_repeats": 51,
                    "contig_name": "contig_51",
                    "cluster_contigs": ["contig_51"],
                },
            }
            alleles_json = tmp_path / "alleles.json"
            alleles_json.write_text(json.dumps(alleles_data))

            runner = CliRunner()
            result = runner.invoke(
                main,
                [
                    "call",
                    "--input",
                    str(bam),
                    "--reference",
                    str(ref),
                    "--alleles-json",
                    str(alleles_json),
                    "--output-dir",
                    str(tmp_path),
                ],
            )

        assert result.exit_code == 0, result.output


class TestConsensusSubcommand:
    """Tests for the consensus subcommand with mocked tools."""

    def test_consensus_invokes_build(self, tmp_path):
        """consensus subcommand runs without error when tools are mocked."""
        flanks = 10
        fake_fasta = f">contig_51\n{'A' * (flanks * 2 + 8)}\n"

        with (
            patch("open_pacmuci.tools.check_tools", return_value=True),
            patch("open_pacmuci.consensus.run_tool", return_value=fake_fasta),
        ):
            ref = tmp_path / "ref.fa"
            ref.touch()
            alleles_data = {
                "homozygous": True,
                "allele_1": {
                    "length": 60,
                    "reads": 200,
                    "canonical_repeats": 51,
                    "contig_name": "contig_51",
                    "cluster_contigs": ["contig_51"],
                },
                "allele_2": {
                    "length": 60,
                    "reads": 0,
                    "canonical_repeats": 51,
                    "contig_name": "contig_51",
                    "cluster_contigs": ["contig_51"],
                },
            }
            alleles_json = tmp_path / "alleles.json"
            alleles_json.write_text(json.dumps(alleles_data))
            # Create a dummy VCF so the subcommand picks it up
            vcf_dir = tmp_path / "allele_1"
            vcf_dir.mkdir()
            (vcf_dir / "variants.vcf.gz").touch()
            bam = tmp_path / "mapping.bam"
            bam.touch()

            runner = CliRunner()
            result = runner.invoke(
                main,
                [
                    "consensus",
                    "--input",
                    str(bam),
                    "--reference",
                    str(ref),
                    "--alleles-json",
                    str(alleles_json),
                    "--output-dir",
                    str(tmp_path),
                ],
            )

        assert result.exit_code == 0, result.output
