"""Mocked unit tests for mapping module (minimap2/samtools wrappers)."""

from __future__ import annotations

from unittest.mock import patch

from open_pacmuci.mapping import bam_to_fastq, get_idxstats, map_reads


class TestBamToFastq:
    """Tests for bam_to_fastq."""

    @patch("open_pacmuci.mapping.run_tool")
    def test_calls_samtools_fastq(self, mock_run_tool, tmp_path):
        """bam_to_fastq calls samtools fastq with the BAM path."""
        mock_run_tool.return_value = "@read1\nACGT\n+\nIIII\n"
        bam = tmp_path / "input.bam"
        bam.touch()

        bam_to_fastq(bam, tmp_path / "out")

        mock_run_tool.assert_called_once_with(["samtools", "fastq", str(bam)])

    @patch("open_pacmuci.mapping.run_tool")
    def test_writes_fastq_output(self, mock_run_tool, tmp_path):
        """bam_to_fastq writes samtools stdout to extracted_reads.fq."""
        fastq_content = "@read1\nACGT\n+\nIIII\n"
        mock_run_tool.return_value = fastq_content
        bam = tmp_path / "input.bam"
        bam.touch()
        out_dir = tmp_path / "out"

        result = bam_to_fastq(bam, out_dir)

        assert result == out_dir / "extracted_reads.fq"
        assert result.read_text() == fastq_content

    @patch("open_pacmuci.mapping.run_tool")
    def test_creates_output_dir(self, mock_run_tool, tmp_path):
        """bam_to_fastq creates the output directory if it does not exist."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "input.bam"
        bam.touch()
        out_dir = tmp_path / "nested" / "out"

        assert not out_dir.exists()
        bam_to_fastq(bam, out_dir)
        assert out_dir.exists()


class TestMapReads:
    """Tests for map_reads."""

    @patch("open_pacmuci.mapping.run_tool")
    def test_fastq_input_pipeline(self, mock_run_tool, tmp_path):
        """For FASTQ input, minimap2 → samtools sort → samtools index are called."""
        mock_run_tool.return_value = ""
        fastq = tmp_path / "reads.fq"
        fastq.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()
        out_dir = tmp_path / "out"

        map_reads(fastq, ref, out_dir, threads=2)

        calls = mock_run_tool.call_args_list
        # First call: minimap2
        assert calls[0][0][0][0] == "minimap2"
        assert "-x" in calls[0][0][0]
        assert "map-hifi" in calls[0][0][0]
        # Second call: samtools sort
        assert calls[1][0][0][:2] == ["samtools", "sort"]
        # Third call: samtools index
        assert calls[2][0][0][:2] == ["samtools", "index"]

    @patch("open_pacmuci.mapping.run_tool")
    def test_bam_input_converts_to_fastq_first(self, mock_run_tool, tmp_path):
        """For BAM input, samtools fastq is called before minimap2."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "reads.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()
        out_dir = tmp_path / "out"

        map_reads(bam, ref, out_dir, threads=2)

        first_cmd = mock_run_tool.call_args_list[0][0][0]
        assert first_cmd[:2] == ["samtools", "fastq"]

    @patch("open_pacmuci.mapping.run_tool")
    def test_returns_bam_path(self, mock_run_tool, tmp_path):
        """map_reads returns the sorted BAM path."""
        mock_run_tool.return_value = ""
        fastq = tmp_path / "reads.fq"
        fastq.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        result = map_reads(fastq, ref, tmp_path, threads=1)

        assert result == tmp_path / "mapping.bam"

    @patch("open_pacmuci.mapping.run_tool")
    def test_threads_passed_to_minimap2(self, mock_run_tool, tmp_path):
        """The threads parameter is passed to minimap2 -t flag."""
        mock_run_tool.return_value = ""
        fastq = tmp_path / "reads.fq"
        fastq.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        map_reads(fastq, ref, tmp_path, threads=8)

        minimap2_cmd = mock_run_tool.call_args_list[0][0][0]
        t_idx = minimap2_cmd.index("-t")
        assert minimap2_cmd[t_idx + 1] == "8"

    @patch("open_pacmuci.mapping.run_tool")
    def test_sam_file_removed_after_sorting(self, mock_run_tool, tmp_path):
        """Intermediate SAM file is cleaned up after samtools sort."""
        mock_run_tool.return_value = ""
        fastq = tmp_path / "reads.fq"
        fastq.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        # SAM file should not exist afterwards (unlink called)
        map_reads(fastq, ref, tmp_path, threads=1)
        assert not (tmp_path / "mapping.sam").exists()


class TestGetIdxstats:
    """Tests for get_idxstats."""

    @patch("open_pacmuci.mapping.run_tool")
    def test_calls_samtools_idxstats(self, mock_run_tool, tmp_path):
        """get_idxstats calls samtools idxstats with the BAM path."""
        idxstats_output = "contig_60\t4120\t245\t0\n*\t0\t0\t50\n"
        mock_run_tool.return_value = idxstats_output
        bam = tmp_path / "mapping.bam"

        result = get_idxstats(bam)

        mock_run_tool.assert_called_once_with(["samtools", "idxstats", str(bam)])
        assert result == idxstats_output

    @patch("open_pacmuci.mapping.run_tool")
    def test_returns_raw_output(self, mock_run_tool, tmp_path):
        """get_idxstats returns the raw string from run_tool unchanged."""
        raw = "contig_51\t3180\t115\t0\n"
        mock_run_tool.return_value = raw
        bam = tmp_path / "mapping.bam"

        assert get_idxstats(bam) == raw
