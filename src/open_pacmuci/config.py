# src/open_pacmuci/config.py
"""Configuration module for repeat dictionary loading and classification."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

# Default path: from src/open_pacmuci/config.py, go up three levels to project root
_BUNDLED_REPEATS = (
    Path(__file__).parent.parent.parent / "data" / "repeats" / "repeats.json"
)


@dataclass
class RepeatDictionary:
    """Holds the MUC1-VNTR repeat type definitions and flanking sequences.

    Attributes:
        repeats: Mapping of repeat ID to 60-bp consensus sequence.
        repeat_length_bp: Expected length in base pairs for each repeat unit.
        pre_repeat_ids: Ordered IDs for repeats that precede the main VNTR array.
        after_repeat_ids: Ordered IDs for repeats that follow the main VNTR array.
        canonical_repeat: ID of the canonical repeat type (usually "X").
        flanking_left: Left flanking sequence (hg38).
        flanking_right: Right flanking sequence (hg38).
        vntr_region: Genomic coordinates of the VNTR region (hg38).
        source: Provenance string for the repeat definitions.
    """

    repeats: dict[str, str]
    repeat_length_bp: int
    pre_repeat_ids: list[str]
    after_repeat_ids: list[str]
    canonical_repeat: str
    flanking_left: str
    flanking_right: str
    vntr_region: str
    source: str


def load_repeat_dictionary(path: Path | None = None) -> RepeatDictionary:
    """Load the repeat dictionary from a JSON file.

    Args:
        path: Path to a JSON file.  If *None*, the bundled
              ``data/repeats/repeats.json`` is used.

    Returns:
        A :class:`RepeatDictionary` instance populated from the file.

    Raises:
        FileNotFoundError: If the JSON file does not exist.
        KeyError: If required keys are missing from the JSON.
    """
    resolved = path if path is not None else _BUNDLED_REPEATS
    if not resolved.exists():
        raise FileNotFoundError(f"Repeat dictionary not found: {resolved}")

    data = json.loads(resolved.read_text(encoding="utf-8"))

    flanking = data["flanking_hg38"]
    return RepeatDictionary(
        repeats=data["repeats"],
        repeat_length_bp=data["repeat_length_bp"],
        pre_repeat_ids=data["pre_repeat_ids"],
        after_repeat_ids=data["after_repeat_ids"],
        canonical_repeat=data["canonical_repeat"],
        flanking_left=flanking["left"],
        flanking_right=flanking["right"],
        vntr_region=flanking["vntr_region"],
        source=data["source"],
    )


def classify_repeat_id(repeat_id: str, dictionary: RepeatDictionary | None = None) -> str:
    """Classify a repeat ID into a structural category.

    Categories:
    - ``"pre"``       - precedes the main VNTR array (e.g. "1" to "5C")
    - ``"after"``     - follows the main VNTR array (e.g. "6" to "9")
    - ``"canonical"`` - the canonical repeat type (usually "X")
    - ``"variable"``  - any other known repeat variant (e.g. "A", "B", ...)

    Args:
        repeat_id: The repeat identifier to classify.
        dictionary: Optional :class:`RepeatDictionary` to use.  If *None*, the
                    bundled dictionary is loaded automatically.

    Returns:
        One of ``"pre"``, ``"after"``, ``"canonical"``, or ``"variable"``.
    """
    if dictionary is None:
        dictionary = load_repeat_dictionary()

    if repeat_id in dictionary.pre_repeat_ids:
        return "pre"
    if repeat_id in dictionary.after_repeat_ids:
        return "after"
    if repeat_id == dictionary.canonical_repeat:
        return "canonical"
    return "variable"
