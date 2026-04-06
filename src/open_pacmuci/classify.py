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


def split_into_repeats(sequence: str, unit_length: int = 60) -> list[str]:
    """Split a consensus sequence into repeat-sized windows.

    Args:
        sequence: The DNA consensus sequence.
        unit_length: Expected repeat unit length in bp (default 60).

    Returns:
        List of sequence windows. Partial trailing windows (< unit_length)
        are included if they are at least half the unit length.
    """
    if not sequence:
        return []

    windows = []
    for i in range(0, len(sequence), unit_length):
        window = sequence[i : i + unit_length]
        # Include partial windows if at least half length
        if len(window) >= unit_length // 2:
            windows.append(window)

    return windows


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
    # Try exact match first
    for repeat_id, ref_seq in repeat_dict.repeats.items():
        if sequence == ref_seq:
            return {"type": repeat_id, "match": "exact"}

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
        "differences": diffs,
    }

    if has_indels:
        result["classification"] = "mutation"
        result["frameshift"] = is_frameshift
    else:
        result["classification"] = "novel_repeat" if best_dist > 2 else "variant"

    return result


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

    pos = 0
    repeat_index = 0
    cumulative_offset = 0

    while pos < len(sequence):
        repeat_index += 1

        remaining = len(sequence) - pos
        if remaining < unit_length // 2:
            break

        # Probe multiple window sizes around unit_length to find the
        # best-matching window.  A dupC (1bp insertion) makes the repeat
        # 61bp; a 14bp deletion makes it 46bp.  We try all sizes from
        # (unit_length - max_indel) to (unit_length + max_indel) and
        # pick the one with the lowest edit distance.
        best_result: dict | None = None
        best_dist = float("inf")
        best_window_size = unit_length

        for probe_size in range(
            max(unit_length - max_indel_probe, unit_length // 2),
            min(unit_length + max_indel_probe + 1, remaining + 1),
        ):
            window = sequence[pos : pos + probe_size]
            result = classify_repeat(window, repeat_dict)

            dist = 0 if result["match"] == "exact" else result.get("edit_distance", 999)
            if dist < best_dist:
                best_dist = dist
                best_result = result
                best_window_size = probe_size

            # Early exit on exact match
            if dist == 0:
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

    return {
        "structure": " ".join(labels),
        "repeats": repeats,
        "mutations_detected": mutations,
        "cumulative_offset": cumulative_offset,
    }
