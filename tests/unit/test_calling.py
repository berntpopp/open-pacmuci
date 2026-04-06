"""Mocked unit tests for calling module (Clair3 / bcftools wrappers)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from open_pacmuci.calling import (
    _extract_and_remap_reads,
    call_variants_per_allele,
    disambiguate_same_length_alleles,
    extract_allele_reads,
    filter_vcf,
    parse_vcf_genotypes,
    run_clair3,
)


class TestExtractAlleleReads:
    """Tests for extract_allele_reads."""

    @patch("open_pacmuci.calling.run_tool")
    def test_single_contig_string(self, mock_run_tool, tmp_path):
        """Accepts a single contig name as a string and calls samtools view."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()

        extract_allele_reads(bam, "contig_51", tmp_path / "out")

        view_call = mock_run_tool.call_args_list[0][0][0]
        assert view_call[:2] == ["samtools", "view"]
        assert "contig_51" in view_call

    @patch("open_pacmuci.calling.run_tool")
    def test_list_of_contigs(self, mock_run_tool, tmp_path):
        """Accepts a list of contig names and passes them all to samtools view."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        contigs = ["contig_50", "contig_51", "contig_52"]

        extract_allele_reads(bam, contigs, tmp_path / "out")

        view_call = mock_run_tool.call_args_list[0][0][0]
        for c in contigs:
            assert c in view_call

    @patch("open_pacmuci.calling.run_tool")
    def test_returns_allele_bam_path(self, mock_run_tool, tmp_path):
        """Returns path to allele_reads.bam inside output dir."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        out_dir = tmp_path / "out"

        result = extract_allele_reads(bam, "contig_51", out_dir)

        assert result == out_dir / "allele_reads.bam"

    @patch("open_pacmuci.calling.run_tool")
    def test_indexes_result_bam(self, mock_run_tool, tmp_path):
        """samtools index is called on the output BAM."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        out_dir = tmp_path / "out"

        extract_allele_reads(bam, "contig_51", out_dir)

        index_call = mock_run_tool.call_args_list[1][0][0]
        assert index_call[:2] == ["samtools", "index"]

    @patch("open_pacmuci.calling.run_tool")
    def test_creates_output_dir(self, mock_run_tool, tmp_path):
        """Creates the output directory if it does not exist."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        out_dir = tmp_path / "nested" / "out"

        assert not out_dir.exists()
        extract_allele_reads(bam, "contig_51", out_dir)
        assert out_dir.exists()


class TestRunClair3:
    """Tests for run_clair3."""

    @patch("open_pacmuci.calling.run_tool")
    def test_builds_correct_command(self, mock_run_tool, tmp_path):
        """run_clair3 passes required flags to run_clair3.sh."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "allele.bam"
        ref = tmp_path / "ref.fa"
        out_dir = tmp_path / "clair3"

        run_clair3(bam, ref, out_dir, platform="hifi", threads=4)

        cmd = mock_run_tool.call_args[0][0]
        assert cmd[0] == "run_clair3.sh"
        assert f"--bam_fn={bam}" in cmd
        assert f"--ref_fn={ref}" in cmd
        assert "--include_all_ctgs" in cmd
        assert "--platform=hifi" in cmd
        assert "--threads=4" in cmd

    @patch("open_pacmuci.calling.run_tool")
    def test_model_path_appended_when_given(self, mock_run_tool, tmp_path):
        """--model_path flag is added when model_path is non-empty."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "allele.bam"
        ref = tmp_path / "ref.fa"
        out_dir = tmp_path / "clair3"

        run_clair3(bam, ref, out_dir, model_path="/models/hifi")

        cmd = mock_run_tool.call_args[0][0]
        assert "--model_path=/models/hifi" in cmd

    @patch("open_pacmuci.calling.run_tool")
    def test_model_path_omitted_when_empty(self, mock_run_tool, tmp_path):
        """--model_path is not added when model_path is empty string."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "allele.bam"
        ref = tmp_path / "ref.fa"
        out_dir = tmp_path / "clair3"

        run_clair3(bam, ref, out_dir, model_path="")

        cmd = mock_run_tool.call_args[0][0]
        assert not any(a.startswith("--model_path") for a in cmd)

    @patch("open_pacmuci.calling.run_tool")
    def test_returns_vcf_path(self, mock_run_tool, tmp_path):
        """run_clair3 returns path to merge_output.vcf.gz."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "allele.bam"
        ref = tmp_path / "ref.fa"
        out_dir = tmp_path / "clair3"

        result = run_clair3(bam, ref, out_dir)

        assert result == out_dir / "merge_output.vcf.gz"

    @patch("open_pacmuci.calling.run_tool")
    def test_creates_output_dir(self, mock_run_tool, tmp_path):
        """run_clair3 creates the output directory if it does not exist."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "allele.bam"
        ref = tmp_path / "ref.fa"
        out_dir = tmp_path / "clair3" / "nested"

        assert not out_dir.exists()
        run_clair3(bam, ref, out_dir)
        assert out_dir.exists()


class TestFilterVcf:
    """Tests for filter_vcf."""

    @patch("open_pacmuci.calling.run_tool")
    def test_calls_bcftools_norm_then_view_then_index(self, mock_run_tool, tmp_path):
        """filter_vcf runs bcftools norm, bcftools view, then bcftools index."""
        mock_run_tool.return_value = ""
        vcf = tmp_path / "raw.vcf.gz"
        ref = tmp_path / "ref.fa"
        out_dir = tmp_path / "filtered"

        filter_vcf(vcf, ref, out_dir)

        calls = mock_run_tool.call_args_list
        assert calls[0][0][0][:2] == ["bcftools", "norm"]
        assert calls[1][0][0][:2] == ["bcftools", "view"]
        assert calls[2][0][0][:2] == ["bcftools", "index"]

    @patch("open_pacmuci.calling.run_tool")
    def test_norm_uses_reference(self, mock_run_tool, tmp_path):
        """bcftools norm -f flag receives the reference path."""
        mock_run_tool.return_value = ""
        vcf = tmp_path / "raw.vcf.gz"
        ref = tmp_path / "ref.fa"
        out_dir = tmp_path / "filtered"

        filter_vcf(vcf, ref, out_dir)

        norm_cmd = mock_run_tool.call_args_list[0][0][0]
        assert "-f" in norm_cmd
        ref_idx = norm_cmd.index("-f")
        assert norm_cmd[ref_idx + 1] == str(ref)

    @patch("open_pacmuci.calling.run_tool")
    def test_view_filters_pass(self, mock_run_tool, tmp_path):
        """bcftools view uses -f PASS to keep only passing variants."""
        mock_run_tool.return_value = ""
        vcf = tmp_path / "raw.vcf.gz"
        ref = tmp_path / "ref.fa"
        out_dir = tmp_path / "filtered"

        filter_vcf(vcf, ref, out_dir)

        view_cmd = mock_run_tool.call_args_list[1][0][0]
        assert "-f" in view_cmd
        f_idx = view_cmd.index("-f")
        assert view_cmd[f_idx + 1] == "PASS"

    @patch("open_pacmuci.calling.run_tool")
    def test_returns_variants_vcf_path(self, mock_run_tool, tmp_path):
        """filter_vcf returns path to variants.vcf.gz."""
        mock_run_tool.return_value = ""
        vcf = tmp_path / "raw.vcf.gz"
        ref = tmp_path / "ref.fa"
        out_dir = tmp_path / "filtered"

        result = filter_vcf(vcf, ref, out_dir)

        assert result == out_dir / "variants.vcf.gz"


class TestFilterVcfEmptyVcf:
    """Tests for filter_vcf handling empty VCFs (issue #8)."""

    @patch("open_pacmuci.calling.run_tool")
    def test_empty_vcf_skips_quality_filter(self, mock_run_tool, tmp_path):
        """Empty VCF (no records) skips -i filter to avoid INFO/DP crash."""
        # bcftools norm returns empty output, then bcftools view should NOT
        # include -i filter since there are no records to filter
        norm_vcf = tmp_path / "normalized.vcf.gz"

        def side_effect(cmd):
            # After norm runs, create the normalized file with header only
            if cmd[:2] == ["bcftools", "norm"]:
                norm_vcf.write_bytes(b"")  # empty file
            return ""

        mock_run_tool.side_effect = side_effect
        vcf = tmp_path / "input.vcf.gz"
        vcf.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        filter_vcf(vcf, ref, tmp_path, min_qual=15.0, min_dp=5)

        # bcftools view should NOT have -i flag (empty VCF)
        view_calls = [
            c[0][0]
            for c in mock_run_tool.call_args_list
            if c[0][0][:2] == ["bcftools", "view"]
        ]
        assert len(view_calls) == 1
        assert "-i" not in view_calls[0]

    @patch("open_pacmuci.calling.run_tool")
    def test_nonempty_vcf_uses_qual_only_filter(self, mock_run_tool, tmp_path):
        """Non-empty VCF uses QUAL filter only (no INFO/DP which may not exist)."""
        norm_vcf = tmp_path / "normalized.vcf.gz"

        def side_effect(cmd):
            if cmd[:2] == ["bcftools", "norm"]:
                # Write a non-empty file to simulate records
                norm_vcf.write_bytes(b"\x00" * 100)
            return ""

        mock_run_tool.side_effect = side_effect
        vcf = tmp_path / "input.vcf.gz"
        vcf.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        filter_vcf(vcf, ref, tmp_path, min_qual=15.0, min_dp=5)

        view_calls = [
            c[0][0]
            for c in mock_run_tool.call_args_list
            if c[0][0][:2] == ["bcftools", "view"]
        ]
        assert len(view_calls) == 1
        view_cmd = view_calls[0]
        assert "-i" in view_cmd
        i_idx = view_cmd.index("-i")
        expr = view_cmd[i_idx + 1]
        assert "QUAL" in expr
        # Should NOT reference INFO/DP (may not exist in Clair3 output)
        assert "INFO/DP" not in expr


class TestCallVariantsPerAllele:
    """Tests for call_variants_per_allele."""

    def _make_alleles_heterozygous(self):
        return {
            "homozygous": False,
            "allele_1": {
                "length": 60,
                "reads": 200,
                "canonical_repeats": 51,
                "contig_name": "contig_51",
                "cluster_contigs": ["contig_50", "contig_51", "contig_52"],
            },
            "allele_2": {
                "length": 80,
                "reads": 150,
                "canonical_repeats": 71,
                "contig_name": "contig_71",
                "cluster_contigs": ["contig_70", "contig_71", "contig_72"],
            },
        }

    def _make_alleles_homozygous(self):
        return {
            "homozygous": True,
            "allele_1": {
                "length": 60,
                "reads": 400,
                "canonical_repeats": 51,
                "contig_name": "contig_51",
                "cluster_contigs": ["contig_50", "contig_51", "contig_52"],
            },
            "allele_2": {
                "length": 60,
                "reads": 0,
                "canonical_repeats": 51,
                "contig_name": "contig_51",
                "cluster_contigs": ["contig_50", "contig_51", "contig_52"],
            },
        }

    @patch("open_pacmuci.calling.run_tool")
    def test_heterozygous_processes_both_alleles(self, mock_run_tool, tmp_path):
        """For a heterozygous sample, both allele_1 and allele_2 are processed."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()
        alleles = self._make_alleles_heterozygous()

        result = call_variants_per_allele(bam, ref, alleles, tmp_path)

        assert "allele_1" in result
        assert "allele_2" in result

    @patch("open_pacmuci.calling.run_tool")
    def test_homozygous_skips_allele_2(self, mock_run_tool, tmp_path):
        """For a homozygous sample, allele_2 is skipped."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()
        alleles = self._make_alleles_homozygous()

        result = call_variants_per_allele(bam, ref, alleles, tmp_path)

        assert "allele_1" in result
        assert "allele_2" not in result

    @patch("open_pacmuci.calling.run_tool")
    def test_returns_dict_of_vcf_paths(self, mock_run_tool, tmp_path):
        """Results map allele keys to Path objects."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()
        alleles = self._make_alleles_heterozygous()

        result = call_variants_per_allele(bam, ref, alleles, tmp_path)

        for _key, path in result.items():
            assert isinstance(path, Path)

    @patch("open_pacmuci.calling.run_tool")
    def test_fallback_contig_name_from_length(self, mock_run_tool, tmp_path):
        """When contig_name is absent, falls back to contig_<length>."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()
        # Alleles dict without explicit contig_name (legacy format)
        alleles = {
            "homozygous": True,
            "allele_1": {
                "length": 60,
                "reads": 200,
                "canonical_repeats": 51,
                "cluster_contigs": ["contig_51"],
                # no "contig_name" key
            },
            "allele_2": {
                "length": 60,
                "reads": 0,
                "canonical_repeats": 51,
                "cluster_contigs": ["contig_51"],
            },
        }

        # Should not raise
        result = call_variants_per_allele(bam, ref, alleles, tmp_path)
        assert "allele_1" in result


class TestFilterVcfQuality:
    """Tests for VCF quality filter parameters."""

    @patch("open_pacmuci.calling.run_tool")
    def test_filter_vcf_includes_quality_expression(self, mock_run_tool, tmp_path):
        """filter_vcf passes QUAL filter to bcftools view for non-empty VCFs."""
        norm_vcf = tmp_path / "normalized.vcf.gz"

        def side_effect(cmd):
            if cmd[:2] == ["bcftools", "norm"]:
                norm_vcf.write_bytes(b"\x00" * 100)  # non-empty
            return ""

        mock_run_tool.side_effect = side_effect
        vcf = tmp_path / "input.vcf.gz"
        vcf.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        filter_vcf(vcf, ref, tmp_path, min_qual=15.0, min_dp=5)

        # Find the bcftools view call
        view_calls = [
            c[0][0] for c in mock_run_tool.call_args_list if c[0][0][:2] == ["bcftools", "view"]
        ]
        assert len(view_calls) == 1
        view_cmd = view_calls[0]
        assert "-i" in view_cmd
        i_idx = view_cmd.index("-i")
        expr = view_cmd[i_idx + 1]
        assert "QUAL" in expr

    @patch("open_pacmuci.calling.run_tool", return_value="")
    def test_filter_vcf_default_params(self, mock_run_tool, tmp_path):
        """filter_vcf works with default parameters (backward compatible)."""
        vcf = tmp_path / "input.vcf.gz"
        vcf.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        # Should not raise with no extra args
        filter_vcf(vcf, ref, tmp_path)


class TestExtractAndRemapReads:
    """Tests for the private _extract_and_remap_reads helper."""

    @patch("open_pacmuci.calling.run_tool")
    def test_remaps_to_peak_contig(self, mock_run_tool, tmp_path):
        """Remapping pipeline is: extract→fastq→faidx (extract contig)→faidx (index)→minimap2→sort→index."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        _extract_and_remap_reads(
            bam,
            ["contig_50", "contig_51", "contig_52"],
            "contig_51",
            ref,
            tmp_path / "out",
            threads=2,
        )

        tool_names = [c[0][0][0] for c in mock_run_tool.call_args_list]
        assert "minimap2" in tool_names
        # samtools sort should be present
        sort_calls = [
            c[0][0] for c in mock_run_tool.call_args_list if c[0][0][:2] == ["samtools", "sort"]
        ]
        assert sort_calls

    @patch("open_pacmuci.calling.run_tool")
    def test_creates_contig_fasta_for_reference(self, mock_run_tool, tmp_path):
        """faidx is called to extract the peak contig as a mini-reference."""
        mock_run_tool.return_value = ">contig_51\nACGT\n"
        bam = tmp_path / "mapping.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()
        out_dir = tmp_path / "out"

        _extract_and_remap_reads(bam, ["contig_51"], "contig_51", ref, out_dir, threads=1)

        faidx_calls = [
            c[0][0] for c in mock_run_tool.call_args_list if c[0][0][:2] == ["samtools", "faidx"]
        ]
        # At least one faidx call with the peak contig name
        contig_faidx = [c for c in faidx_calls if "contig_51" in c]
        assert contig_faidx


class TestParseVcfGenotypes:
    """Tests for parse_vcf_genotypes."""

    @patch("open_pacmuci.calling.run_tool")
    def test_parses_genotype_fields(self, mock_run_tool):
        mock_run_tool.return_value = "contig_51\t100\tA\tT\t0/1\n"
        result = parse_vcf_genotypes(Path("/fake.vcf"))
        assert len(result) == 1
        assert result[0]["pos"] == 100
        assert result[0]["genotype"] == "0/1"

    @patch("open_pacmuci.calling.run_tool")
    def test_handles_runtime_error(self, mock_run_tool):
        mock_run_tool.side_effect = RuntimeError("fail")
        assert parse_vcf_genotypes(Path("/fake.vcf")) == []


class TestDisambiguateSameLengthAlleles:
    """Tests for disambiguate_same_length_alleles."""

    @patch("open_pacmuci.calling.run_tool", return_value="")
    @patch("open_pacmuci.calling.parse_vcf_genotypes", return_value=[])
    def test_no_het_variants_returns_homozygous(self, mock_geno, mock_run, tmp_path):
        alleles = {
            "allele_1": {"contig_name": "contig_51", "cluster_contigs": ["contig_51"]},
            "allele_2": {"contig_name": "contig_51", "cluster_contigs": ["contig_51"]},
        }
        result = disambiguate_same_length_alleles(
            tmp_path / "bam", tmp_path / "ref.fa", alleles, tmp_path
        )
        assert "allele_1" in result
        assert result.get("homozygous") is True

    @patch("open_pacmuci.calling.run_tool", return_value="")
    @patch(
        "open_pacmuci.calling.parse_vcf_genotypes",
        return_value=[{"chrom": "c", "pos": 100, "ref": "A", "alt": "T", "genotype": "0/1"}],
    )
    def test_het_variants_returns_two_alleles(self, mock_geno, mock_run, tmp_path):
        alleles = {
            "allele_1": {"contig_name": "contig_51", "cluster_contigs": ["contig_51"]},
            "allele_2": {"contig_name": "contig_51", "cluster_contigs": ["contig_51"]},
        }
        result = disambiguate_same_length_alleles(
            tmp_path / "bam", tmp_path / "ref.fa", alleles, tmp_path
        )
        assert "allele_1" in result
        assert "allele_2" in result
        assert result.get("homozygous") is False
