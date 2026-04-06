# tests/unit/test_tools.py
"""Tests for subprocess helper utilities."""

from __future__ import annotations

import pytest

from open_pacmuci.tools import check_tools, run_tool


class TestRunTool:
    """Tests for the run_tool function."""

    def test_run_tool_success(self):
        """run_tool returns stdout for a successful command."""
        result = run_tool(["echo", "hello"])
        assert result.strip() == "hello"

    def test_run_tool_failure_raises(self):
        """run_tool raises RuntimeError on non-zero exit."""
        with pytest.raises(RuntimeError, match="failed with exit code"):
            run_tool(["false"])

    def test_run_tool_not_found_raises(self):
        """run_tool raises FileNotFoundError for missing commands."""
        with pytest.raises(FileNotFoundError):
            run_tool(["nonexistent_tool_xyz"])

    def test_run_tool_captures_stderr(self):
        """run_tool includes stderr in error message."""
        with pytest.raises(RuntimeError, match="No such file"):
            run_tool(["ls", "/nonexistent_path_xyz"])


class TestCheckTools:
    """Tests for the check_tools function."""

    def test_check_tools_all_present(self):
        """check_tools returns True when all tools are found."""
        assert check_tools(["echo", "ls"]) is True

    def test_check_tools_missing_tool(self):
        """check_tools raises RuntimeError listing missing tools."""
        with pytest.raises(RuntimeError, match="nonexistent_tool_xyz"):
            check_tools(["echo", "nonexistent_tool_xyz"])

    def test_check_tools_empty_list(self):
        """check_tools with empty list returns True."""
        assert check_tools([]) is True
