"""Smoke tests for the CLI interface."""

from __future__ import annotations

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
        assert "0.1.0" in result.output

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
