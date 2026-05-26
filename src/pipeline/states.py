"""
Breathing-state generators for nucleosome free energy calculations.

Three styles are supported:
  - "b_index"    → index pairs with length 14  (fully wrapped = (0, 13))
  - "ph_index"   → index pairs with length 28  (fully wrapped = (0, 27))
  - "open_sites" → open-site pairs that sum ≤ length
"""

from __future__ import annotations

from typing import List, Tuple


def total_states_index(length: int) -> List[Tuple[int, int]]:
    """
    Generate all (left, right) binding-index pairs where
    0 <= left <= right < length.

    For "b_index" style, length=14  → fully wrapped = (0, 13).
    For "ph_index" style, length=28 → fully wrapped = (0, 27).
    """
    states: List[Tuple[int, int]] = []
    for left in range(length):
        for right in range(left, length):
            states.append((left, right))
    return states


def total_open_states(length: int = 14) -> List[Tuple[int, int]]:
    """
    Generate all (left, right) open-site pairs where
    left >= 0, right >= 0, and left + right <= length.
    """
    open_states: List[Tuple[int, int]] = []
    for left in range(length + 1):
        for right in range(length + 1):
            if left + right <= length:
                open_states.append((left, right))
    return open_states


# ── convenience map ───────────────────────────────────────────────────────────
_STYLE_DEFAULTS = {
    "b_index":    (total_states_index, 14),
    "ph_index":   (total_states_index, 28),
    "open_sites": (total_open_states,  14),
}

_FULLY_BOUND = {
    "b_index":    [(0, 13)],
    "ph_index":   [(0, 27)],
    "open_sites": [(0, 0)],
}


def get_states(
    style: str,
    only_fullbound: bool = False,
) -> List[Tuple[int, int]]:
    """
    Return the list of breathing states for a given style.

    Parameters
    ----------
    style : str
        One of "b_index", "ph_index", "open_sites".
    only_fullbound : bool
        If True, return only the single fully-wrapped state.

    Returns
    -------
    List[Tuple[int, int]]
    """
    if style not in _STYLE_DEFAULTS:
        raise ValueError(
            f"Unknown style '{style}'. Choose from {list(_STYLE_DEFAULTS)}"
        )
    if only_fullbound:
        return _FULLY_BOUND[style]

    gen_fn, length = _STYLE_DEFAULTS[style]
    return gen_fn(length)
