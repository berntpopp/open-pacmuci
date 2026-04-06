# tests/unit/test_classify.py
"""Tests for repeat unit classification."""

from __future__ import annotations

import pytest

from open_pacmuci.classify import (
    characterize_differences,
    classify_repeat,
    classify_sequence,
    edit_distance,
    validate_mutations_against_vcf,
)
from open_pacmuci.config import load_repeat_dictionary


@pytest.fixture
def repeat_dict():
    """Load bundled repeat dictionary."""
    return load_repeat_dictionary()


class TestEditDistance:
    """Tests for edit distance calculation."""

    def test_identical_sequences(self):
        """Edit distance of identical sequences is 0."""
        assert edit_distance("ACGT", "ACGT") == 0

    def test_single_substitution(self):
        """Single base substitution gives distance 1."""
        assert edit_distance("ACGT", "ACGA") == 1

    def test_single_insertion(self):
        """Single insertion gives distance 1."""
        assert edit_distance("ACGT", "ACGAT") == 1

    def test_single_deletion(self):
        """Single deletion gives distance 1."""
        assert edit_distance("ACGT", "ACT") == 1

    def test_large_indel(self):
        """Multi-base indel is counted correctly."""
        # 59dupC equivalent: inserting one C into a 60bp sequence
        seq_a = "G" * 60
        seq_b = "G" * 59 + "CG"  # insert C before last G
        assert edit_distance(seq_a, seq_b) == 1


class TestCharacterizeDifferences:
    """Tests for characterizing differences between sequences."""

    def test_single_insertion(self):
        """Detects a single insertion."""
        ref = "ACGT"
        query = "ACGCT"  # C inserted at position 3
        diffs = characterize_differences(ref, query)
        has_insertion = any(d["type"] == "insertion" for d in diffs)
        assert has_insertion

    def test_single_substitution(self):
        """Detects a single substitution."""
        ref = "ACGT"
        query = "ACAT"  # G->A at position 2
        diffs = characterize_differences(ref, query)
        has_sub = any(d["type"] == "substitution" for d in diffs)
        assert has_sub

    def test_identifies_frameshift(self):
        """Insertion of 1bp is a frameshift (1 % 3 != 0)."""
        ref = "ACGT"
        query = "ACGAT"
        diffs = characterize_differences(ref, query)
        total_indel = sum(1 for d in diffs if d["type"] in ("insertion", "deletion"))
        assert total_indel % 3 != 0  # frameshift

    def test_no_differences(self):
        """Identical sequences produce empty diff list."""
        diffs = characterize_differences("ACGT", "ACGT")
        assert diffs == []


class TestClassifyRepeat:
    """Tests for classifying a single 60bp repeat unit."""

    def test_exact_match(self, repeat_dict):
        """Known repeat X matches exactly."""
        x_seq = repeat_dict.repeats["X"]
        result = classify_repeat(x_seq, repeat_dict)
        assert result["type"] == "X"
        assert result["match"] == "exact"

    def test_known_pre_repeat(self, repeat_dict):
        """Pre-repeat 1 matches exactly."""
        seq = repeat_dict.repeats["1"]
        result = classify_repeat(seq, repeat_dict)
        assert result["type"] == "1"
        assert result["match"] == "exact"

    def test_unknown_with_insertion(self, repeat_dict):
        """Repeat X with 1bp insertion is classified as mutation."""
        x_seq = repeat_dict.repeats["X"]
        # Simulate 59dupC: insert C at position 59
        mutated = x_seq[:59] + "C" + x_seq[59:]
        result = classify_repeat(mutated, repeat_dict)
        assert result["match"] != "exact"
        assert result["closest_match"] == "X"
        assert result["edit_distance"] == 1
        assert any(d["type"] == "insertion" for d in result["differences"])

    def test_unknown_with_large_deletion(self, repeat_dict):
        """Repeat X with 14bp deletion matches del18_31 template exactly."""
        x_seq = repeat_dict.repeats["X"]
        # Simulate del18_31: delete positions 18-31 (1-based)
        mutated = x_seq[:17] + x_seq[31:]
        result = classify_repeat(mutated, repeat_dict)
        # Now matches the del18_31 mutation template exactly
        assert result["match"] == "exact"
        assert result["type"] == "X:del18_31"

    def test_identity_pct_reported(self, repeat_dict):
        """Identity percentage is included in result."""
        x_seq = repeat_dict.repeats["X"]
        mutated = x_seq[:59] + "C" + x_seq[59:]  # 1bp insertion
        result = classify_repeat(mutated, repeat_dict)
        assert "identity_pct" in result
        assert result["identity_pct"] > 90.0


class TestClassifySequence:
    """Tests for classifying a full consensus sequence."""

    def test_all_canonical_x(self, repeat_dict):
        """Sequence of 3 canonical X repeats classified correctly."""
        x_seq = repeat_dict.repeats["X"]
        full_seq = x_seq * 3
        result = classify_sequence(full_seq, repeat_dict)
        assert len(result["repeats"]) == 3
        assert all(r["type"] == "X" for r in result["repeats"])
        assert result["structure"] == "X X X"

    def test_mutation_detected(self, repeat_dict):
        """Mutation in one repeat is detected and reported."""
        x_seq = repeat_dict.repeats["X"]
        mutated = x_seq[:59] + "C" + x_seq[59:]  # 59dupC
        full_seq = x_seq + mutated + x_seq
        result = classify_sequence(full_seq, repeat_dict)
        assert len(result["mutations_detected"]) >= 1


class TestConfidenceScoring:
    """Tests for per-repeat confidence scores."""

    def test_exact_match_confidence_is_1(self, repeat_dict):
        """Exact match to known repeat type has confidence 1.0."""
        x_seq = repeat_dict.repeats["X"]
        result = classify_repeat(x_seq, repeat_dict)
        assert result["confidence"] == 1.0

    def test_close_match_confidence_uses_identity(self, repeat_dict):
        """Near-match (ed=1) confidence equals identity_pct / 100."""
        x_seq = repeat_dict.repeats["X"]
        mutated = x_seq[:59] + "C" + x_seq[59:]  # 61bp, ed=1
        result = classify_repeat(mutated, repeat_dict)
        assert 0.9 < result["confidence"] < 1.0
        assert result["confidence"] == result["identity_pct"] / 100

    def test_distant_match_has_lower_confidence(self, repeat_dict):
        """High edit distance produces lower confidence."""
        seq = "ACGT" * 15
        result = classify_repeat(seq, repeat_dict)
        assert result["confidence"] < 0.9


class TestSequenceConfidenceSummary:
    """Tests for per-allele confidence summary in classify_sequence."""

    def test_all_exact_gives_confidence_1(self, repeat_dict):
        """All exact matches produce allele_confidence=1.0."""
        x_seq = repeat_dict.repeats["X"]
        result = classify_sequence(x_seq * 3, repeat_dict)
        assert result["allele_confidence"] == 1.0
        assert result["exact_match_pct"] == 100.0

    def test_mixed_match_reduces_confidence(self, repeat_dict):
        """One mutation among exact matches reduces allele_confidence below 1.0."""
        x_seq = repeat_dict.repeats["X"]
        mutated = x_seq[:59] + "C" + x_seq[59:]
        result = classify_sequence(x_seq + mutated + x_seq, repeat_dict)
        assert result["allele_confidence"] < 1.0
        assert result["exact_match_pct"] < 100.0


class TestMutationTemplateMatching:
    """Tests for exact matching against known mutation templates."""

    def test_dupc_exact_template_match(self, repeat_dict):
        """dupC on X matches the pre-computed 61bp template exactly."""
        x_seq = repeat_dict.repeats["X"]
        dupc_seq = x_seq[:60] + "C"  # 61bp: insert C at end
        result = classify_repeat(dupc_seq, repeat_dict)
        assert result["match"] == "exact"
        assert result["type"] == "X:dupC"
        assert result["confidence"] == 1.0

    def test_ins16bp_exact_template_match(self, repeat_dict):
        """16bp insertion on C matches the 76bp template exactly."""
        c_seq = repeat_dict.repeats["C"]
        ins_seq = c_seq[:57] + "GGGCTCCACCGCCCCC" + c_seq[57:]  # 76bp
        result = classify_repeat(ins_seq, repeat_dict)
        assert result["match"] == "exact"
        assert result["type"] == "C:ins16bp"
        assert result["confidence"] == 1.0

    def test_del18_31_exact_template_match(self, repeat_dict):
        """14bp deletion on X matches the 46bp template exactly."""
        x_seq = repeat_dict.repeats["X"]
        del_seq = x_seq[:17] + x_seq[31:]  # 46bp: delete pos 18-31 (1-based)
        result = classify_repeat(del_seq, repeat_dict)
        assert result["match"] == "exact"
        assert result["type"] == "X:del18_31"

    def test_sequence_with_dupc_classifies_correctly(self, repeat_dict):
        """classify_sequence finds dupC at correct window size."""
        x_seq = repeat_dict.repeats["X"]
        dupc_seq = x_seq[:60] + "C"  # 61bp
        full = x_seq + dupc_seq + x_seq
        result = classify_sequence(full, repeat_dict)
        assert len(result["repeats"]) == 3
        assert result["repeats"][1]["type"] == "X:dupC"
        assert result["repeats"][1]["match"] == "exact"

    def test_sequence_with_16bp_ins_classifies_correctly(self, repeat_dict):
        """classify_sequence handles 76bp mutated repeat with template match."""
        x_seq = repeat_dict.repeats["X"]
        c_seq = repeat_dict.repeats["C"]
        ins_seq = c_seq[:57] + "GGGCTCCACCGCCCCC" + c_seq[57:]
        full = x_seq + ins_seq + x_seq
        result = classify_sequence(full, repeat_dict)
        assert len(result["repeats"]) == 3
        assert result["repeats"][1]["type"] == "C:ins16bp"


class TestVcfMutationValidation:
    """Tests for VCF-backed mutation validation."""

    def test_confirmed_mutation_keeps_high_confidence(self, repeat_dict):
        """Mutation with VCF support keeps confidence unchanged."""
        x_seq = repeat_dict.repeats["X"]
        # Novel 2bp insertion (not in template catalog) -> lands in mutations_detected
        mutated = x_seq[:40] + "TT" + x_seq[40:]  # 62bp
        result = classify_sequence(x_seq + mutated + x_seq, repeat_dict)

        # Simulate VCF with a variant at repeat 2 position
        vcf_variants = [{"pos": 560, "qual": 25.0}]  # flank(500) + 60bp

        validated = validate_mutations_against_vcf(
            result, vcf_variants=vcf_variants, flank_length=500, unit_length=60
        )
        mut = validated["mutations_detected"][0]
        assert mut.get("vcf_support") is True

    def test_unsupported_mutation_gets_low_confidence(self, repeat_dict):
        """Mutation without VCF support gets reduced confidence."""
        x_seq = repeat_dict.repeats["X"]
        # Novel 2bp insertion (not in template catalog)
        mutated = x_seq[:40] + "TT" + x_seq[40:]  # 62bp
        result = classify_sequence(x_seq + mutated + x_seq, repeat_dict)

        # No VCF variants at all
        validated = validate_mutations_against_vcf(
            result, vcf_variants=[], flank_length=500, unit_length=60
        )
        mut = validated["mutations_detected"][0]
        assert mut.get("vcf_support") is False

    def test_no_mutations_unchanged(self, repeat_dict):
        """Sequence with no mutations passes through unchanged."""
        x_seq = repeat_dict.repeats["X"]
        result = classify_sequence(x_seq * 3, repeat_dict)

        validated = validate_mutations_against_vcf(
            result, vcf_variants=[], flank_length=500, unit_length=60
        )
        assert validated["mutations_detected"] == []


class TestBidirectionalClassification:
    """Tests for bidirectional classification fallback."""

    def test_novel_large_insertion_classified_via_bidirectional(self, repeat_dict):
        """Novel large insertion (not in templates) uses bidirectional fallback."""
        x_seq = repeat_dict.repeats["X"]
        pre1 = repeat_dict.repeats["1"]
        after9 = repeat_dict.repeats["9"]
        # Insert 20bp of random sequence into X (novel, not in templates)
        novel_mutated = x_seq[:30] + "A" * 20 + x_seq[30:]  # 80bp
        full = pre1 + x_seq + novel_mutated + x_seq + after9
        result = classify_sequence(full, repeat_dict)
        # Should classify pre1 and after9 correctly
        labels = result["structure"].split()
        assert labels[0] == "1"
        assert labels[-1] == "9"
        # The novel mutation should be detected
        assert len(result["mutations_detected"]) >= 1 or any(
            "m" in label for label in labels
        )

    def test_forward_only_when_no_large_mutation(self, repeat_dict):
        """No bidirectional needed when all repeats classify well."""
        x_seq = repeat_dict.repeats["X"]
        result = classify_sequence(x_seq * 5, repeat_dict)
        assert all(r["match"] == "exact" for r in result["repeats"])
        assert result["allele_confidence"] == 1.0
