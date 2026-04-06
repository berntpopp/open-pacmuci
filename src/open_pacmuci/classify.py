# src/open_pacmuci/classify.py
"""Repeat unit classification and mutation detection."""

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

    Args:
        sequence: Full consensus sequence (flanking regions should be trimmed).
        repeat_dict: The loaded repeat dictionary.

    Returns:
        Dict with structure string, per-repeat details, and mutation report.
    """
    windows = split_into_repeats(sequence, repeat_dict.repeat_length_bp)

    repeats: list[dict] = []
    mutations: list[dict] = []
    labels: list[str] = []

    for i, window in enumerate(windows):
        result = classify_repeat(window, repeat_dict)
        result["index"] = i + 1  # 1-based indexing

        if result["match"] == "exact":
            labels.append(result["type"])
        elif result.get("classification") == "mutation":
            labels.append(f"{result['closest_match']}*")
            mutations.append(
                {
                    "repeat_index": i + 1,
                    "closest_type": result["closest_match"],
                    "differences": result["differences"],
                    "frameshift": result.get("frameshift", False),
                }
            )
        else:
            labels.append(f"?{result.get('closest_match', '?')}")

        repeats.append(result)

    return {
        "structure": " ".join(labels),
        "repeats": repeats,
        "mutations_detected": mutations,
    }
