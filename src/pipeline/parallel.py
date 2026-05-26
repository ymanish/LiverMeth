"""
Parallel free energy calculation with progress tracking.

Uses ``concurrent.futures.ProcessPoolExecutor`` so each worker creates and
caches its own ``NucleosomeBreathModular`` once.  Terminal MN reversion is
applied in the worker before energy calculation.  Progress is displayed with
``tqdm``, showing active workers and pending batches.
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

from .sliding_window import Subsequence, _revert_terminal_mn

# ── per-process worker state (set by _init_worker) ──────────────────────────
_WORKER_CALC = None
_WORKER_STYLE: str = "b_index"
_WORKER_STATES: List[Tuple[int, int]] = []
_WORKER_SKIP_NO_METH: bool = False


def _init_worker(
    parameterization: str,
    cgna_parameter_set: str,
    cgna_group_split: bool,
    nuc_method: str,
    free_dna_method: Optional[str],
    style: str,
    style_states: List[Tuple[int, int]],
    skip_if_no_meth: bool,
) -> None:
    """Pool initializer — create the calculator once per worker process."""
    global _WORKER_CALC, _WORKER_STYLE, _WORKER_STATES, _WORKER_SKIP_NO_METH

    from femodules.nucleosome_breath_modular import NucleosomeBreathModular
    from femodules.config import RBPConfig, CgnaConfig

    if parameterization == "cgna":
        config = CgnaConfig(
            parameter_set_name=cgna_parameter_set,
            group_split=cgna_group_split,
        )
    else:
        config = RBPConfig(
            nuc_method=nuc_method,
            free_dna_method=free_dna_method,
        )

    _WORKER_CALC = NucleosomeBreathModular(config=config)
    _WORKER_STYLE = style
    _WORKER_STATES = style_states
    _WORKER_SKIP_NO_METH = skip_if_no_meth


_ENERGY_COLS = ("F", "F_entropy", "F_enthalpy", "F_freedna", "dF")
_NAN_ENERGY = {c: float("nan") for c in _ENERGY_COLS}
_BAR_FMT = "{desc} |{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}] {postfix}"


def _process_batch(batch: List[Subsequence]) -> List[Dict[str, Any]]:
    """Worker: terminal MN reversion → energy calc for each window × state."""
    results: List[Dict[str, Any]] = []

    for sub in batch:
        seq = _revert_terminal_mn(sub.sequence)
        has_meth = "M" in seq or "N" in seq

        for left_s, right_s in _WORKER_STATES:
            base = {
                "id": sub.id,
                "subid": sub.sub_id,
                "start": sub.start,
                "end": sub.end,
                "left": left_s,
                "right": right_s,
            }
            if has_meth or not _WORKER_SKIP_NO_METH:
                res = _WORKER_CALC.calculate_free_energy(
                    sequence=seq,
                    left=left_s,
                    right=right_s,
                    id=sub.id,
                    subid=sub.sub_id,
                    style=_WORKER_STYLE,
                )
                base.update({
                    "F": float(res.F),
                    "F_entropy": float(res.F_entropy),
                    "F_enthalpy": float(res.F_enthalpy),
                    "F_freedna": float(res.F_freedna),
                    "dF": float(res.F) - float(res.F_freedna),
                })
            else:
                base.update(_NAN_ENERGY)
            results.append(base)

    return results


# ── public helpers ───────────────────────────────────────────────────────────

def filter_unchanged_windows(
    windows: List[Subsequence],
    ref_windows: List[Subsequence],
) -> tuple[list[Subsequence], list[Subsequence]]:
    """
    Partition *windows* into **changed** / **unchanged** relative to
    *ref_windows*.  Comparison is by ``sub_id`` on the ``.sequence`` field
    (which has already been through terminal MN reversion during generation).

    Returns ``(changed, unchanged)``.
    """
    ref_seqs = {w.sub_id: w.sequence for w in ref_windows}

    changed: list[Subsequence] = []
    unchanged: list[Subsequence] = []

    for w in windows:
        ref_seq = ref_seqs.get(w.sub_id)
        if ref_seq is not None and w.sequence == ref_seq:
            unchanged.append(w)
        else:
            changed.append(w)

    return changed, unchanged


def copy_energies_from_reference(
    windows: List[Subsequence],
    df_ref: pd.DataFrame,
    style_states: List[Tuple[int, int]],
) -> pd.DataFrame:
    """
    Build a result DataFrame for *windows* by copying energy values from
    *df_ref*, matched on ``(_region, subid, left, right)``.

    Use for windows known to be identical to a previously computed condition.
    """
    cols = ["id", "subid", "start", "end", "left", "right"] + list(_ENERGY_COLS)
    if not windows:
        return pd.DataFrame(columns=cols)

    # Metadata rows (one per window × state)
    rows: list[dict[str, Any]] = []
    for w in windows:
        for left_s, right_s in style_states:
            rows.append({
                "id": w.id, "subid": w.sub_id,
                "start": w.start, "end": w.end,
                "left": left_s, "right": right_s,
            })

    df = pd.DataFrame(rows)

    # Lookup energies from reference
    energy_cols = list(_ENERGY_COLS)
    df["_region"] = _extract_region(df["id"])

    ref = df_ref.copy()
    ref["_region"] = _extract_region(ref["id"])
    key_cols = ["_region", "subid", "left", "right"]
    ref_lookup = ref.set_index(key_cols)[energy_cols]

    idx = pd.MultiIndex.from_frame(df[key_cols])
    filled = ref_lookup.reindex(idx)
    filled.index = df.index
    for col in energy_cols:
        df[col] = filled[col].values

    df.drop(columns="_region", inplace=True)
    return df


# ── core compute ─────────────────────────────────────────────────────────────

def run_pool_energy(
    batches: List[List[Subsequence]],
    style: str,
    style_states: List[Tuple[int, int]],
    parameterization: str = "cgna",
    cgna_parameter_set: str = "Di_hmethyl_methylated-hemi_combine",
    cgna_group_split: bool = True,
    nuc_method: str = "crystal",
    free_dna_method: Optional[str] = None,
    n_workers: Optional[int] = None,
    skip_if_no_meth: bool = False,
    label: str = "",
) -> List[Dict[str, Any]]:
    """
    Parallel free energy calculation with a live progress bar.

    Parameters
    ----------
    batches : List[List[Subsequence]]
        Pre-batched subsequences.
    style, style_states
        Breathing style and state list.
    parameterization, cgna_parameter_set, cgna_group_split,
    nuc_method, free_dna_method
        Calculator configuration forwarded to worker initializer.
    n_workers : int | None
        Pool size (``None`` → ``cpu_count``).
    skip_if_no_meth : bool
        When *True*, windows without M/N emit NaN (back-fill later).
    label : str
        Short label for the progress bar (e.g. ``"Un"``, ``"d0"``).

    Returns
    -------
    List[Dict[str, Any]]
    """
    if not batches:
        return []

    total = len(batches)

    if n_workers is None:
        n_workers = os.cpu_count() or 1
    n_workers = min(n_workers, total)

    init_args = (
        parameterization, cgna_parameter_set, cgna_group_split,
        nuc_method, free_dna_method, style, style_states, skip_if_no_meth,
    )

    all_results: List[Dict[str, Any]] = []
    desc = f"  {label:>3s}" if label else " calc"

    if n_workers <= 1:
        # ── serial ──
        _init_worker(*init_args)
        with tqdm(total=total, desc=desc, unit="batch",
                  bar_format=_BAR_FMT) as pbar:
            for i, batch in enumerate(batches):
                all_results.extend(_process_batch(batch))
                rem = total - i - 1
                pbar.set_postfix_str(f"active=1/1 pending={rem}")
                pbar.update(1)
    else:
        # ── parallel ──
        with ProcessPoolExecutor(
            max_workers=n_workers,
            initializer=_init_worker,
            initargs=init_args,
        ) as executor:
            futures = [executor.submit(_process_batch, b) for b in batches]
            with tqdm(total=total, desc=desc, unit="batch",
                      bar_format=_BAR_FMT) as pbar:
                for future in as_completed(futures):
                    all_results.extend(future.result())
                    done = pbar.n + 1
                    remaining = total - done
                    active = min(n_workers, remaining)
                    pending = max(0, remaining - active)
                    pbar.set_postfix_str(
                        f"active={active}/{n_workers} pending={pending}"
                    )
                    pbar.update(1)

    return all_results


# ── DataFrame helpers ────────────────────────────────────────────────────────

def results_to_dataframe(results: List[Dict[str, Any]]) -> pd.DataFrame:
    """Flat result list → sorted DataFrame."""
    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values(["id", "subid"]).reset_index(drop=True)
    return df


def _extract_region(id_col: pd.Series) -> pd.Series:
    """``'chr17:100-200_d0'`` → ``'chr17:100-200'``."""
    return id_col.str.split("_").str[0]


def fill_from_reference(
    df: pd.DataFrame,
    df_ref: pd.DataFrame,
) -> pd.DataFrame:
    """
    Back-fill NaN energy rows in *df* from *df_ref*.

    Matching on ``(_region, subid, left, right)`` where ``_region`` is the
    chromosomal coordinate prefix of the ``id`` column.

    Returns a new DataFrame; inputs are not modified.
    """
    if df.empty or df_ref.empty:
        return df.copy()

    energy_cols = list(_ENERGY_COLS)

    ref = df_ref.copy()
    ref["_region"] = _extract_region(ref["id"])
    key_cols = ["_region", "subid", "left", "right"]
    ref_lookup = ref.set_index(key_cols)[energy_cols]

    out = df.copy()
    out["_region"] = _extract_region(out["id"])

    nan_mask = out["F"].isna()
    if not nan_mask.any():
        out.drop(columns="_region", inplace=True)
        return out

    nan_keys = out.loc[nan_mask, key_cols]
    fill_idx = pd.MultiIndex.from_frame(nan_keys)
    filled = ref_lookup.reindex(fill_idx)
    filled.index = nan_keys.index

    out.loc[nan_mask, energy_cols] = filled.values
    out.drop(columns="_region", inplace=True)
    return out


# Backward-compatible alias
fill_from_unmethylated = fill_from_reference