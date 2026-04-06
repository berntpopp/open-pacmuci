# src/open_pacmuci/classify.py
"""Repeat unit classification and mutation detection.

The classification algorithm handles frameshifted sequences by tracking
the cumulative indel offset.  When a repeat contains an insertion or
deletion, subsequent window boundaries are shifted by the net indel
length so that downstream repeats are correctly framed.

For example, a dupC (1bp insertion) at repeat 25 shifts all windows
after repeat 25 by +1bp.  Without correction, every downstream window
would straddle two repeat boundaries and fail to match any known type.
With correction, the windows realign to the true repeat boundaries and
classify correctly.
"""

from __future__ import annotations

from open_pacmuci.config import RepeatDictionary


def edit_distance(s1: str, s2: str) -> int:
    """Compute Levenshtein edit distance between two sequences.

    Args:
        s1: First sequence.
        s2: Second sequence.

    Returns:
        Minimum number of single-character edits (insert, delete, substitute).
    """
    m, n = len(s1), len(s2)

    # Use single-row optimization for memory efficiency
    prev = list(range(n + 1))
    curr = [0] * (n + 1)

    for i in range(1, m + 1):
        curr[0] = i
        for j in range(1, n + 1):
            if s1[i - 1] == s2[j - 1]:
                curr[j] = prev[j - 1]
            else:
                curr[j] = 1 + min(prev[j], curr[j - 1], prev[j - 1])
        prev, curr = curr, prev

    return prev[n]


def characterize_differences(ref: str, query: str) -> list[dict]:
    """Characterize specific differences between reference and query sequences.

    Uses Needleman-Wunsch style traceback to identify individual
    substitutions, insertions, and deletions with positions.

    Args:
        ref: Reference sequence.
        query: Query sequence.

    Returns:
        List of difference dicts with keys: pos, ref, alt, type.
    """
    if ref == query:
        return []

    m, n = len(ref), len(query)

    # Build full DP matrix for traceback
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if ref[i - 1] == query[j - 1]:
                dp[i][j] = dp[i - 1][j - 1]
            else:
                dp[i][j] = 1 + min(dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1])

    # Traceback
    diffs: list[dict] = []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and ref[i - 1] == query[j - 1]:
            i -= 1
            j -= 1
        elif i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + 1:
            # Substitution
            diffs.append(
                {
                    "pos": i,  # 1-based position in reference
                    "ref": ref[i - 1],
                    "alt": query[j - 1],
                    "type": "substitution",
                }
            )
            i -= 1
            j -= 1
        elif j > 0 and dp[i][j] == dp[i][j - 1] + 1:
            # Insertion in query
            diffs.append(
                {
                    "pos": i + 1,  # position after which insertion occurs
                    "ref": "",
                    "alt": query[j - 1],
                    "type": "insertion",
                }
            )
            j -= 1
        elif i > 0 and dp[i][j] == dp[i - 1][j] + 1:
            # Deletion from reference
            diffs.append(
                {
                    "pos": i,
                    "ref": ref[i - 1],
                    "alt": "",
                    "type": "deletion",
                }
            )
            i -= 1
        else:
            break

    diffs.reverse()
    return diffs


def _compute_net_indel(diffs: list[dict]) -> int:
    """Compute net indel offset from a list of differences.

    Insertions add bases (positive offset), deletions remove bases
    (negative offset).  The net offset tells us how much the sequence
    shifted relative to the reference frame.

    Returns:
        Net indel length: positive = sequence is longer, negative = shorter.
    """
    offset = 0
    for d in diffs:
        if d["type"] == "insertion":
            offset += len(d["alt"])
        elif d["type"] == "deletion":
            offset -= len(d["ref"])
    return offset


def classify_repeat(
    sequence: str,
    repeat_dict: RepeatDictionary,
) -> dict:
    """Classify a single repeat unit against the known dictionary.

    Args:
        sequence: The 60bp (or near-60bp) repeat sequence.
        repeat_dict: The loaded repeat dictionary.

    Returns:
        Classification result dict with type, match, and (for unknowns)
        closest_match, edit_distance, identity_pct, differences.
    """
    # O(1) exact-match lookup via reverse map (sequence -> ID)
    seq_to_id = {seq: rid for rid, seq in repeat_dict.repeats.items()}
    if sequence in seq_to_id:
        return {"type": seq_to_id[sequence], "match": "exact", "confidence": 1.0}

    # Check mutation templates (variable-length exact matches)
    if hasattr(repeat_dict, "mutated_sequences") and sequence in repeat_dict.mutated_sequences:
        parent_repeat, mut_name = repeat_dict.mutated_sequences[sequence]
        return {
            "type": f"{parent_repeat}:{mut_name}",
            "match": "exact",
            "confidence": 1.0,
            "mutation_name": mut_name,
            "parent_repeat": parent_repeat,
        }

    # No exact match -- find closest by edit distance
    best_id = ""
    best_dist = float("inf")

    for repeat_id, ref_seq in repeat_dict.repeats.items():
        dist = edit_distance(ref_seq, sequence)
        if dist < best_dist:
            best_dist = dist
            best_id = repeat_id

    # Characterize the specific differences
    ref_seq = repeat_dict.repeats[best_id]
    diffs = characterize_differences(ref_seq, sequence)

    # Calculate identity percentage based on alignment length
    max_len = max(len(ref_seq), len(sequence))
    identity_pct = round((1 - best_dist / max_len) * 100, 1) if max_len > 0 else 0.0

    # Determine if differences contain indels
    has_indels = any(d["type"] in ("insertion", "deletion") for d in diffs)

    # Calculate total indel length for frameshift check
    indel_bases = sum(
        len(d["alt"]) if d["type"] == "insertion" else len(d["ref"])
        for d in diffs
        if d["type"] in ("insertion", "deletion")
    )
    is_frameshift = has_indels and (indel_bases % 3 != 0)

    result: dict = {
        "type": "unknown",
        "match": "closest",
        "closest_match": best_id,
        "edit_distance": best_dist,
        "identity_pct": identity_pct,
        "confidence": identity_pct / 100,
        "differences": diffs,
    }

    if has_indels:
        result["classification"] = "mutation"
        result["frameshift"] = is_frameshift
    else:
        result["classification"] = "novel_repeat" if best_dist > 2 else "variant"

    return result


def _probe_sizes_generator(
    unit_length: int, max_indel_probe: int, remaining: int
) -> list[int]:
    """Generate probe sizes: canonical first, then small-to-large."""
    sizes = [min(unit_length, remaining)]
    for ps in range(
        max(unit_length - max_indel_probe, unit_length // 2),
        min(unit_length + max_indel_probe + 1, remaining + 1),
    ):
        if ps != unit_length:
            sizes.append(ps)
    return sizes


def classify_sequence(
    sequence: str,
    repeat_dict: RepeatDictionary,
) -> dict:
    """Classify all repeat units in a consensus sequence.

    Uses offset-aware windowing: when a repeat contains an indel, the
    cumulative offset is tracked and subsequent window boundaries are
    shifted accordingly.  This corrects for frameshift propagation --
    a 1bp insertion at repeat 25 would otherwise misalign all downstream
    windows.

    Algorithm:
        1. Start at position 0 with offset = 0
        2. Extract window of ``unit_length + offset`` bases (the mutated
           repeat is longer/shorter than 60bp)
        3. Classify the window
        4. If classification finds indels, compute the net offset and
           accumulate it for subsequent windows
        5. Advance position by ``unit_length + net_indel`` (actual length
           of the repeat in the sequence)
        6. Reset offset to 0 for the next window (each downstream repeat
           is expected to be 60bp again, just starting from the shifted
           position)

    Args:
        sequence: Full consensus sequence (flanking regions should be trimmed).
        repeat_dict: The loaded repeat dictionary.

    Returns:
        Dict with structure string, per-repeat details, and mutation report.
    """
    unit_length = repeat_dict.repeat_length_bp
    # Maximum indel size to probe.  Covers all known MUC1 mutations
    # (largest known: 25bp insertion, 14bp deletion).
    max_indel_probe = 30
    repeats: list[dict] = []
    mutations: list[dict] = []
    labels: list[str] = []

    # Build reverse maps once for O(1) exact-match lookups
    seq_to_id = {seq: rid for rid, seq in repeat_dict.repeats.items()}

    pos = 0
    repeat_index = 0
    cumulative_offset = 0

    while pos < len(sequence):
        repeat_index += 1

        remaining = len(sequence) - pos
        if remaining < unit_length // 2:
            break

        best_result: dict | None = None
        best_dist = float("inf")
        best_window_size = unit_length

        # --- Phase 1: Check ALL probe sizes for exact match ---
        # First check canonical size for standard repeats (common case).
        # Then check ALL sizes for mutation templates (which are non-60bp).
        # Mutation templates take priority over canonical-size standard matches
        # because they explain the actual biological repeat length.
        exact_found = False
        canonical_result: dict | None = None

        for probe_size in _probe_sizes_generator(unit_length, max_indel_probe, remaining):
            window = sequence[pos : pos + probe_size]
            # Check mutation templates first (variable-length exact matches)
            if (
                hasattr(repeat_dict, "mutated_sequences")
                and window in repeat_dict.mutated_sequences
            ):
                parent, mname = repeat_dict.mutated_sequences[window]
                best_result = {
                    "type": f"{parent}:{mname}",
                    "match": "exact",
                    "confidence": 1.0,
                    "mutation_name": mname,
                    "parent_repeat": parent,
                }
                best_window_size = probe_size
                best_dist = 0
                exact_found = True
                break
            # Check standard repeats
            if window in seq_to_id:
                if canonical_result is None:
                    canonical_result = {
                        "type": seq_to_id[window],
                        "match": "exact",
                        "confidence": 1.0,
                    }

        # Use canonical match if no mutation template found
        if not exact_found and canonical_result is not None:
            best_result = canonical_result
            best_window_size = unit_length
            best_dist = 0
            exact_found = True

        # --- Phase 2: Edit distance fallback (only if no exact match) ---
        if not exact_found:
            # Try canonical size first
            if remaining >= unit_length:
                window = sequence[pos : pos + unit_length]
                result = classify_repeat(window, repeat_dict)
                dist = (
                    0 if result["match"] == "exact" else result.get("edit_distance", 999)
                )
                best_dist = dist
                best_result = result
                best_window_size = unit_length

            if best_dist > 0:
                for probe_size in range(
                    max(unit_length - max_indel_probe, unit_length // 2),
                    min(unit_length + max_indel_probe + 1, remaining + 1),
                ):
                    if probe_size == unit_length:
                        continue
                    window = sequence[pos : pos + probe_size]
                    result = classify_repeat(window, repeat_dict)
                    dist = (
                        0
                        if result["match"] == "exact"
                        else result.get("edit_distance", 999)
                    )
                    if dist < best_dist:
                        best_dist = dist
                        best_result = result
                        best_window_size = probe_size
                    if best_dist <= 1:
                        break

        assert best_result is not None
        result = best_result
        advance = best_window_size

        # Track cumulative offset if this window is non-standard size
        if best_window_size != unit_length:
            net_indel = best_window_size - unit_length
            cumulative_offset += net_indel

        result["index"] = repeat_index

        if result["match"] == "exact":
            labels.append(result["type"])
        elif result.get("classification") == "mutation":
            # Use MucOneUp nomenclature: "Xm" = repeat X with mutation
            labels.append(f"{result['closest_match']}m")
            mutations.append(
                {
                    "repeat_index": repeat_index,
                    "closest_type": result["closest_match"],
                    "differences": result["differences"],
                    "frameshift": result.get("frameshift", False),
                }
            )
        else:
            # Variant of known type (substitutions only, no indel)
            labels.append(f"?{result.get('closest_match', '?')}")

        repeats.append(result)
        pos += max(advance, 1)  # always advance at least 1 to avoid infinite loop

    confidences = [r.get("confidence", 1.0) for r in repeats]
    exact_count = sum(1 for r in repeats if r.get("match") == "exact")
    allele_confidence = sum(confidences) / len(confidences) if confidences else 0.0
    exact_match_pct = (exact_count / len(repeats) * 100) if repeats else 0.0

    return {
        "structure": " ".join(labels),
        "repeats": repeats,
        "mutations_detected": mutations,
        "cumulative_offset": cumulative_offset,
        "allele_confidence": round(allele_confidence, 4),
        "exact_match_pct": round(exact_match_pct, 1),
    }


def validate_mutations_against_vcf(
    classification_result: dict,
    vcf_variants: list[dict] | None = None,
    flank_length: int = 500,
    unit_length: int = 60,
) -> dict:
    """Cross-validate detected mutations against VCF variant positions.

    For each mutation, check if a VCF variant overlaps the repeat's
    genomic position.  Adds ``vcf_support`` flag and adjusts confidence.

    Args:
        classification_result: Output from :func:`classify_sequence`.
        vcf_variants: List of dicts with ``pos`` (int) and ``qual`` (float).
            If None, VCF validation is skipped (standalone classify mode).
        flank_length: Flanking bp on each side of the contig.
        unit_length: Expected repeat unit length (60bp).

    Returns:
        Updated classification result with VCF validation annotations.
    """
    result = classification_result.copy()
    result["mutations_detected"] = [m.copy() for m in result.get("mutations_detected", [])]
    result["repeats"] = [r.copy() for r in result.get("repeats", [])]

    if vcf_variants is None:
        return result

    for mutation in result["mutations_detected"]:
        repeat_idx = mutation["repeat_index"]
        # Map repeat index to contig coordinates
        repeat_start = flank_length + (repeat_idx - 1) * unit_length
        repeat_end = repeat_start + unit_length + 30  # allow for indels

        # Check if any VCF variant overlaps this repeat
        supporting = [v for v in vcf_variants if repeat_start <= v["pos"] <= repeat_end]
        mutation["vcf_support"] = len(supporting) > 0
        mutation["vcf_qual"] = max((v["qual"] for v in supporting), default=0.0)

        # Adjust confidence in the corresponding repeat
        if repeat_idx - 1 < len(result["repeats"]):
            repeat_result = result["repeats"][repeat_idx - 1]
            base_confidence = repeat_result.get("confidence", 1.0)
            if supporting:
                best_qual = mutation["vcf_qual"]
                vcf_score = 1.0 if best_qual >= 20 else 0.7
            else:
                vcf_score = 0.3
            repeat_result["confidence"] = round(base_confidence * vcf_score, 4)

    # Recompute allele_confidence
    confidences = [r.get("confidence", 1.0) for r in result["repeats"]]
    result["allele_confidence"] = (
        round(sum(confidences) / len(confidences), 4) if confidences else 0.0
    )

    return result
