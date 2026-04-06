# tests/unit/test_config.py
"""Tests for the config module (repeat dictionary loading and classification)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from open_pacmuci.config import (
    RepeatDictionary,
    classify_repeat_id,
    load_repeat_dictionary,
)


class TestLoadRepeatDictionary:
    """Tests for load_repeat_dictionary."""

    def test_loads_bundled_dictionary(self):
        """Bundled dictionary loads without error and contains >= 34 repeats."""
        rd = load_repeat_dictionary()
        assert isinstance(rd, RepeatDictionary)
        assert len(rd.repeats) >= 34
        assert rd.repeat_length_bp == 60

    def test_canonical_repeat_x_exists(self):
        """Repeat 'X' is present in the bundled dictionary with length 60."""
        rd = load_repeat_dictionary()
        assert "X" in rd.repeats
        assert len(rd.repeats["X"]) == 60

    def test_pre_repeats_present(self):
        """Repeats '1' through '5' are present in the bundled dictionary."""
        rd = load_repeat_dictionary()
        for rid in ["1", "2", "3", "4", "5"]:
            assert rid in rd.repeats, f"Pre-repeat '{rid}' missing from dictionary"

    def test_after_repeats_present(self):
        """Repeats '6' through '9' are present in the bundled dictionary."""
        rd = load_repeat_dictionary()
        for rid in ["6", "7", "8", "9"]:
            assert rid in rd.repeats, f"After-repeat '{rid}' missing from dictionary"

    def test_all_repeats_are_60bp(self):
        """Every repeat sequence in the bundled dictionary is exactly 60 bp."""
        rd = load_repeat_dictionary()
        for rid, seq in rd.repeats.items():
            assert len(seq) == 60, f"Repeat '{rid}' has length {len(seq)}, expected 60"

    def test_flanking_sequences_present(self):
        """Left and right flanking sequences are non-empty strings."""
        rd = load_repeat_dictionary()
        assert isinstance(rd.flanking_left, str) and len(rd.flanking_left) > 0
        assert isinstance(rd.flanking_right, str) and len(rd.flanking_right) > 0

    def test_loads_from_custom_path(self, tmp_path: Path):
        """load_repeat_dictionary accepts a custom Path and loads it correctly."""
        custom_data = {
            "source": "test",
            "repeat_length_bp": 60,
            "repeats": {"X": "A" * 60, "A": "C" * 60},
            "flanking_hg38": {
                "left": "LEFTSEQ",
                "right": "RIGHTSEQ",
                "vntr_region": "chr1:100-200",
            },
            "pre_repeat_ids": [],
            "after_repeat_ids": [],
            "canonical_repeat": "X",
        }
        custom_file = tmp_path / "custom_repeats.json"
        custom_file.write_text(json.dumps(custom_data), encoding="utf-8")

        rd = load_repeat_dictionary(path=custom_file)
        assert len(rd.repeats) == 2
        assert rd.repeats["X"] == "A" * 60
        assert rd.flanking_left == "LEFTSEQ"
        assert rd.flanking_right == "RIGHTSEQ"
        assert rd.vntr_region == "chr1:100-200"
        assert rd.source == "test"

    def test_missing_file_raises(self, tmp_path: Path):
        """load_repeat_dictionary raises FileNotFoundError for missing paths."""
        with pytest.raises(FileNotFoundError):
            load_repeat_dictionary(path=tmp_path / "does_not_exist.json")


class TestClassifyRepeatId:
    """Tests for classify_repeat_id."""

    @pytest.fixture
    def dictionary(self) -> RepeatDictionary:
        return load_repeat_dictionary()

    def test_classify_pre_repeat(self, dictionary: RepeatDictionary):
        """Pre-repeat IDs are classified as 'pre'."""
        for rid in ["1", "2", "3", "4", "4p", "5", "5C"]:
            assert classify_repeat_id(rid, dictionary) == "pre", f"Expected 'pre' for '{rid}'"

    def test_classify_after_repeat(self, dictionary: RepeatDictionary):
        """After-repeat IDs are classified as 'after'."""
        for rid in ["6", "6p", "7", "8", "9"]:
            assert classify_repeat_id(rid, dictionary) == "after", f"Expected 'after' for '{rid}'"

    def test_classify_canonical_repeat(self, dictionary: RepeatDictionary):
        """Canonical repeat ID 'X' is classified as 'canonical'."""
        assert classify_repeat_id("X", dictionary) == "canonical"

    def test_classify_variable_repeat(self, dictionary: RepeatDictionary):
        """Known variant repeat IDs (e.g. 'A', 'B') are classified as 'variable'."""
        for rid in ["A", "B", "C", "D"]:
            assert classify_repeat_id(rid, dictionary) == "variable", (
                f"Expected 'variable' for '{rid}'"
            )

    def test_classify_unknown_id_is_variable(self, dictionary: RepeatDictionary):
        """An unrecognised repeat ID falls back to 'variable'."""
        assert classify_repeat_id("UNKNOWN_ZZZ", dictionary) == "variable"

    def test_classify_uses_bundled_dictionary_by_default(self):
        """classify_repeat_id loads the bundled dictionary when none is provided."""
        assert classify_repeat_id("X") == "canonical"
        assert classify_repeat_id("1") == "pre"
        assert classify_repeat_id("6") == "after"
        assert classify_repeat_id("A") == "variable"
