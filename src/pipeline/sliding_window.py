"""
Sliding-window generation and batching for nucleosome energy calculations.
"""

from __future__ import annotations

import itertools
from typing import Generator, NamedTuple


class Subsequence(NamedTuple):
    """A 147 bp window extracted from a larger sequence."""
    id: str          # sequence header / identifier
    sub_id: int      # global sequential index
    start: int       # start position in the reference coordinate system
    end: int         # end position (exclusive)
    sequence: str    # the DNA window string


# ── Terminal MN reversion ────────────────────────────────────────────────────

def _is_terminal_step(i: int, L: int) -> bool:
    """True if dinucleotide step *i* is within the first or last 2 steps."""
    return i == 0 or i == 1 or i == L - 3 or i == L - 2


def _revert_terminal_mn(subseq: str) -> str:
    """
    Revert ``MN → CG`` at terminal positions of a window, including orphaned
    ``M`` or ``N`` bases that arise when a window boundary splits an MN pair.

    Three cases:
      1. Complete MN pair at a terminal step   → both reverted.
      2. Orphaned N at positions 0 or 1        → N → G.
      3. Orphaned M at positions L-2 or L-1    → M → C.
    """
    L = len(subseq)
    seq = list(subseq)

    # case 1 — full MN pair
    for i in (i for i in range(L - 1) if seq[i] == "M" and seq[i + 1] == "N"):
        if _is_terminal_step(i, L):
            seq[i] = "C"
            seq[i + 1] = "G"

    # case 2 — orphaned N at start
    for i in (0, 1):
        if seq[i] == "N":
            seq[i] = "G"

    # case 3 — orphaned M at end
    for i in (L - 2, L - 1):
        if seq[i] == "M":
            seq[i] = "C"

    return "".join(seq)


# ── Window generators ────────────────────────────────────────────────────────

def sliding_window_sequence(
    sequence: str,
    seq_id: str,
    window_size: int = 147,
    step_size: int = 1,
) -> Generator[Subsequence, None, None]:
    """
    Yield sliding-window :class:`Subsequence` tuples from a plain string.

    If the sequence contains MN-encoded methylated CpGs, terminal MN
    dinucleotides are reverted to CG per window before yielding.
    """
    seq_len = len(sequence)
    if seq_len < window_size:
        return

    n_windows = (seq_len - window_size) // step_size + 1
    for i in range(n_windows):
        start = i * step_size
        end = start + window_size
        subseq = _revert_terminal_mn(sequence[start:end])
        yield Subsequence(id=seq_id, sub_id=i, start=start, end=end, sequence=subseq)


def nuc_window_generator(
    full_sequence: str,
    segments: list[tuple[int, int]],
    seq_id: str,
    window_size: int = 147,
    step_size: int = 1,
) -> Generator[Subsequence, None, None]:
    """
    Yield sliding-window :class:`Subsequence` tuples across multiple
    nucleosome-occupied segments.

    ``sub_id`` is a **global** counter across all segments.  ``start`` /
    ``end`` are in full-sequence coordinates.
    """
    global_idx = 0
    for seg_start, seg_end in segments:
        seg_seq = full_sequence[seg_start:seg_end]
        if len(seg_seq) < window_size:
            continue
        for win in sliding_window_sequence(
            seg_seq, seq_id=seq_id, window_size=window_size, step_size=step_size
        ):
            yield Subsequence(
                id=win.id,
                sub_id=global_idx,
                start=seg_start + win.start,
                end=seg_start + win.end,
                sequence=win.sequence,
            )
            global_idx += 1


# ── Batching ─────────────────────────────────────────────────────────────────

def batcher(
    it: Generator,
    size: int,
) -> Generator[list[Subsequence], None, None]:
    """Batch a generator into successive lists of length *size*."""
    it = iter(it)
    for first in it:
        yield list(itertools.chain([first], itertools.islice(it, size - 1)))
