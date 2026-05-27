"""Sequence engineering for targeted methylation response.

See ``docs/superpowers/specs/2026-05-26-sequence-engineering-design.md``.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Iterable, Literal, Optional

import numpy as np
import pandas as pd

from .methylation import apply_methylation_to_sequence, find_cpg_positions
from .sliding_window import _revert_terminal_mn

WINDOW_SIZE = 147
DYAD = 73
DEFAULT_MOTIF = "GGCCCAATTT"  # GC half + AT half, no CpG steps anywhere


def position_map(
    df: pd.DataFrame,
    stab_quantile: float = 0.25,
    destab_quantile: float = 0.75,
    neutral_quantile: float = 0.25,
    ddg_col: str = "ddG",
    pos_col: str = "cpg_pos",
) -> dict:
    """Per-position mean ΔΔG plus stabilizing / destabilizing / neutral pools.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain columns ``pos_col`` and ``ddg_col``.
    stab_quantile, destab_quantile : float
        Quantile thresholds on the per-position mean ΔΔG.
    neutral_quantile : float
        Quantile cut on the per-position |mean ΔΔG|; positions in the
        bottom quantile (closest-to-zero mean) are the neutral pool.

    Returns
    -------
    dict with keys
        - ``'mean_ddg'``: np.ndarray, length WINDOW_SIZE - 1, per-position mean
          (NaN where no CpG was observed at that position).
        - ``'stab'``: np.ndarray of positions in the bottom `stab_quantile`.
        - ``'destab'``: np.ndarray of positions in the top
          `1 - destab_quantile` (i.e., positions with mean above the
          `destab_quantile` cut).
        - ``'neutral'``: np.ndarray of positions in the bottom
          `neutral_quantile` of |mean ΔΔG|.
    """
    mean = (
        df.groupby(pos_col)[ddg_col].mean()
        .reindex(range(WINDOW_SIZE - 1))
    )
    mean_arr = mean.to_numpy()
    finite = mean.dropna()
    stab_cut = finite.quantile(stab_quantile)
    destab_cut = finite.quantile(destab_quantile)
    abs_cut = finite.abs().quantile(neutral_quantile)
    stab = finite[finite <= stab_cut].index.to_numpy(dtype=int)
    destab = finite[finite >= destab_cut].index.to_numpy(dtype=int)
    neutral = finite[finite.abs() <= abs_cut].index.to_numpy(dtype=int)
    return {
        "mean_ddg": mean_arr,
        "stab": np.sort(stab),
        "destab": np.sort(destab),
        "neutral": np.sort(neutral),
    }


def widom_backbone(
    length: int = WINDOW_SIZE,
    dyad: int = DYAD,
    motif: str = DEFAULT_MOTIF,
    seed: Optional[int] = None,
) -> str:
    """Periodic CpG-free nucleosome-favoring backbone.

    The 10-bp ``motif`` is repeated with a phase chosen so the GC-rich
    first half lands centered at ``dyad`` (major groove faces histone at
    SHL 0). ``seed`` performs A↔T swaps within the AT half (the GC half
    is left fixed to avoid introducing CG steps).

    Backbone rule references (cite in the notebook):
      - Zuiddam & Schiessel, Phys. Rev. E 99, 012422 (2019).
      - Pérez-Pulido et al., Nucleic Acids Res. 50, 1864 (2022).
      - Collings et al., Sci. Reports 3, 2121 (2013).
    """
    if len(motif) != 10:
        raise ValueError("motif must be 10 bp")
    if "CG" in (motif + motif):  # check across tile junction too
        raise ValueError("motif must contain no CG dinucleotide")

    # Phase: place midpoint of GC half (motif index 2) at the dyad.
    phase = (dyad - 2) % 10
    seq = [motif[(p - phase) % 10] for p in range(length)]

    if seed is not None:
        rng = np.random.default_rng(seed)
        for p in range(length):
            m_idx = (p - phase) % 10
            if 5 <= m_idx < 10 and rng.random() < 0.5:  # AT half only
                seq[p] = "T" if seq[p] == "A" else "A"
    return "".join(seq)


def inject_cpgs(backbone: str, positions: Iterable[int]) -> str:
    """Return ``backbone`` with ``seq[p:p+2] = 'CG'`` at each ``p``.

    Raises ``ValueError`` if any two positions are adjacent (the CGs
    would overlap and corrupt each other).
    """
    pos = sorted(set(int(p) for p in positions))
    for a, b in zip(pos, pos[1:]):
        if b - a < 2:
            raise ValueError(f"positions {a} and {b} are adjacent")
    seq = list(backbone)
    for p in pos:
        if p + 1 >= len(seq):
            raise ValueError(f"position {p} too close to sequence end")
        seq[p] = "C"
        seq[p + 1] = "G"
    return "".join(seq)


def score_full_meth(seq: str, nb) -> tuple[float, float, float]:
    """Return (G_unmeth, G_meth_full, ddG) for one sequence.

    ``nb`` is a ``NucleosomeBreathModular`` instance. Methylation pattern
    is every CpG dinucleotide in ``seq``.
    """
    cpgs = set(find_cpg_positions(seq))
    seq_m = _revert_terminal_mn(apply_methylation_to_sequence(seq, cpgs))

    res_u = nb.calculate_free_energy(sequence=seq, left=0, right=13, style="b_index")
    res_m = nb.calculate_free_energy(sequence=seq_m, left=0, right=13, style="b_index")

    g_un = float(res_u.F - res_u.F_freedna)
    g_me = float(res_m.F - res_m.F_freedna)
    return g_un, g_me, g_me - g_un


def _build_candidate(args):
    """Worker: build one candidate, score it, return a row dict."""
    backbone, class_name, positions, cgna_config = args
    # NucleosomeBreathModular is non-picklable per-process; lazily build one.
    from femodules.nucleosome_breath_modular import NucleosomeBreathModular
    global _WORKER_NB  # noqa: PLW0603
    try:
        nb = _WORKER_NB
    except NameError:
        nb = NucleosomeBreathModular(cgna_config)
        _WORKER_NB = nb  # type: ignore[name-defined]

    seq = backbone if not positions else inject_cpgs(backbone, positions)
    g_un, g_me, dd = score_full_meth(seq, nb)
    return {
        "class": class_name,
        "n_cpg": len(positions) if positions else len(find_cpg_positions(seq)),
        "cpg_positions": tuple(int(p) for p in positions),
        "G_un": g_un,
        "G_meth": g_me,
        "ddG": dd,
        "seq": seq,
    }


def design_library(
    class_name: Literal["stable", "unstable", "inert"],
    pos_map: dict,
    cgna_config,
    n_candidates: int = 100,
    cpg_count_range: tuple[int, int] = (4, 12),
    n_backbones: int = 20,
    n_workers: int = 8,
    seed: int = 0,
) -> pd.DataFrame:
    """Generate and score ``n_candidates`` sequences for one class.

    For all three classes, CpG positions are sampled from the corresponding
    pool (``pos_map['stab']`` / ``pos_map['destab']`` / ``pos_map['neutral']``)
    with ``n_cpg`` drawn uniformly from ``cpg_count_range`` and respecting
    non-adjacency (resamples until a valid set is drawn). For 'inert',
    the pool is the near-zero mean ΔΔG positions so CpGs are present and
    methylated but contribute roughly zero net effect.
    """
    rng = np.random.default_rng(seed)
    backbones = [widom_backbone(seed=int(s)) for s in rng.integers(1, 10**9, size=n_backbones)]

    pool_key = {"stable": "stab", "unstable": "destab", "inert": "neutral"}[class_name]
    position_pool = pos_map[pool_key]
    if len(position_pool) == 0:
        raise ValueError(f"empty position pool for class {class_name!r}")

    tasks = []
    for _ in range(n_candidates):
        backbone = backbones[int(rng.integers(0, n_backbones))]
        lo, hi = cpg_count_range
        for _attempt in range(50):
            k = int(rng.integers(lo, hi + 1))
            k = min(k, len(position_pool))
            chosen = np.sort(rng.choice(position_pool, size=k, replace=False))
            if all(b - a >= 2 for a, b in zip(chosen, chosen[1:])):
                positions = tuple(int(x) for x in chosen)
                break
        else:
            raise RuntimeError("could not sample non-adjacent positions; "
                               "loosen cpg_count_range or check position pool")
        tasks.append((backbone, class_name, positions, cgna_config))

    rows = []
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        for fut in as_completed(pool.submit(_build_candidate, t) for t in tasks):
            rows.append(fut.result())
    return pd.DataFrame(rows)
