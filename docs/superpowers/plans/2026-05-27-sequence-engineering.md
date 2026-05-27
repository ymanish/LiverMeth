# Sequence Engineering for Targeted Methylation Response — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a library of 147 bp DNA sequences in three classes (stable / unstable / inert under full CpG methylation) using a Widom-style CpG-free backbone plus rule-driven CpG injection, and validate by computing ΔΔG = G_meth_full − G_unmeth for each candidate.

**Architecture:** A new module `src/pipeline/sequence_design.py` holds five focused functions (`position_map`, `widom_backbone`, `inject_cpgs`, `score_full_meth`, `design_library`). A new notebook `notebooks/Sequence_engineering.ipynb` loads the existing scan output, runs the library design + scoring (parallel), and plots class-separation.

**Tech Stack:** `NucleosomeBreathModular` (CGNA+ free-energy soft-binding), `apply_methylation_to_sequence`, `find_cpg_positions`, `ProcessPoolExecutor`, pandas, matplotlib, nbformat.

**Branch:** `bp-level-prop`. Module commits go to git; notebooks and data files are gitignored.

---

## Conventions

- **ΔΔG sign**: `ddG = G_meth_full − G_unmeth`. Negative → stabilizing, positive → destabilizing.
- **Free-energy quantity**: `dF_el = result.F − result.F_freedna` (binding free energy minus free-DNA reference), matching the existing notebook's convention.
- **Binding spec**: `left=0, right=13, style="b_index"` (fully wrapped nucleosome).
- **Window**: `L = 147`, `dyad = 73`.
- **DOF / sequence canonicalisation**: same as existing scan — methylation via `apply_methylation_to_sequence` (M/N encoding); no `_revert_terminal_mn` needed here because `NucleosomeBreathModular.calculate_free_energy` handles its own canonicalisation.

---

## Task 1: Persist scan outputs from the existing notebook

The new notebook needs `df` (per-CpG ΔΔG of the random-sequence scan) and `sequences` (the 50 random sequences with labels) loaded from disk. Append a save cell to the existing notebook and run it once.

**Files:**
- Modify: `notebooks/Modelvalidation_meth_exp+analysis.ipynb` (append one cell — gitignored, no commit)
- Create: `data/random_scan/df.csv`, `data/random_scan/sequences.json` (gitignored)

- [ ] **Step 1: Append a save cell via nbformat**

Run this from the repo root:

```bash
python - <<'PY'
import nbformat, pathlib
nb_path = pathlib.Path("notebooks/Modelvalidation_meth_exp+analysis.ipynb")
nb = nbformat.read(nb_path, as_version=4)
src = '''# Persist scan outputs for downstream sequence-engineering notebook
import json, pathlib
out = pathlib.Path("data/random_scan"); out.mkdir(parents=True, exist_ok=True)
df.to_csv(out / "df.csv", index=False)
with open(out / "sequences.json", "w") as f:
    json.dump(sequences, f)
print(f"wrote {out/'df.csv'}  ({len(df)} rows)")
print(f"wrote {out/'sequences.json'}  ({len(sequences)} sequences)")
'''
nb.cells.append(nbformat.v4.new_code_cell(src))
nbformat.write(nb, nb_path)
print("appended save cell")
PY
```

- [ ] **Step 2: Run the new cell in the existing notebook**

In JupyterLab/VS Code, open `notebooks/Modelvalidation_meth_exp+analysis.ipynb`, scroll to the last cell, and run it. Expected output:

```
wrote data/random_scan/df.csv  (~784 rows)
wrote data/random_scan/sequences.json  (50 sequences)
```

- [ ] **Step 3: Verify files exist**

Run:

```bash
ls -la data/random_scan/
```

Expected: `df.csv` and `sequences.json` listed.

No commit — notebooks and `data/` are gitignored.

---

## Task 2: Create `sequence_design` module — utility functions

Create the module with three pure functions that have no dependency on `NucleosomeBreathModular`: position-map extraction, backbone construction, and CpG injection.

**Files:**
- Create: `src/pipeline/sequence_design.py`
- Modify: `src/pipeline/__init__.py`

- [ ] **Step 1: Create the module file**

Write `src/pipeline/sequence_design.py`:

```python
"""Sequence engineering for targeted methylation response.

See ``docs/superpowers/specs/2026-05-26-sequence-engineering-design.md``.
"""
from __future__ import annotations

from typing import Iterable, Literal, Optional

import numpy as np
import pandas as pd

WINDOW_SIZE = 147
DYAD = 73
DEFAULT_MOTIF = "GGCCCAATTT"  # GC half + AT half, no CpG steps anywhere


def position_map(
    df: pd.DataFrame,
    stab_quantile: float = 0.25,
    destab_quantile: float = 0.75,
    ddg_col: str = "ddG",
    pos_col: str = "cpg_pos",
) -> dict:
    """Per-position mean ΔΔG plus stabilizing / destabilizing position sets.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain columns ``pos_col`` and ``ddg_col``.
    stab_quantile, destab_quantile : float
        Quantile thresholds on the per-position mean ΔΔG.

    Returns
    -------
    dict with keys
        - ``'mean_ddg'``: np.ndarray, length WINDOW_SIZE - 1, per-position mean
          (NaN where no CpG was observed at that position).
        - ``'stab'``: np.ndarray of positions in the bottom `stab_quantile`.
        - ``'destab'``: np.ndarray of positions in the top
          `1 - destab_quantile` (i.e., positions with mean above the
          `destab_quantile` cut).
    """
    mean = (
        df.groupby(pos_col)[ddg_col].mean()
        .reindex(range(WINDOW_SIZE - 1))
    )
    mean_arr = mean.to_numpy()
    finite = mean.dropna()
    stab_cut = finite.quantile(stab_quantile)
    destab_cut = finite.quantile(destab_quantile)
    stab = finite[finite <= stab_cut].index.to_numpy(dtype=int)
    destab = finite[finite >= destab_cut].index.to_numpy(dtype=int)
    return {"mean_ddg": mean_arr, "stab": np.sort(stab), "destab": np.sort(destab)}


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
```

- [ ] **Step 2: Export from `src/pipeline/__init__.py`**

Append to `src/pipeline/__init__.py`:

```python
from .sequence_design import (
    position_map,
    widom_backbone,
    inject_cpgs,
)
```

- [ ] **Step 3: Smoke-check the utilities**

Run:

```bash
python - <<'PY'
from src.pipeline.sequence_design import widom_backbone, inject_cpgs
import numpy as np

bb = widom_backbone()
assert len(bb) == 147
assert "CG" not in bb, "backbone must be CG-free"
print(f"backbone[68:80]: {bb[68:80]}  (GGCCC should center at pos 73)")

bb_seed = widom_backbone(seed=42)
assert "CG" not in bb_seed
assert bb != bb_seed, "seed should jitter"
print("seed=42 differs from default:", sum(a != b for a, b in zip(bb, bb_seed)), "positions")

seq = inject_cpgs(bb, [10, 30, 50])
assert seq[10:12] == "CG" and seq[30:32] == "CG" and seq[50:52] == "CG"
print("inject_cpgs OK")
PY
```

Expected: `backbone[68:80]` shows `GGCCC` centered around position 73, jitter changes ~7-ish positions, inject_cpgs prints OK.

- [ ] **Step 4: Commit**

```bash
git add src/pipeline/sequence_design.py src/pipeline/__init__.py
git commit -m "Add sequence_design module utilities (position_map, widom_backbone, inject_cpgs)"
```

---

## Task 3: Add scoring + library functions to the module

Add the two functions that depend on `NucleosomeBreathModular` for free-energy evaluation.

**Files:**
- Modify: `src/pipeline/sequence_design.py`

- [ ] **Step 1: Append scoring + design_library**

Append to `src/pipeline/sequence_design.py`:

```python
from concurrent.futures import ProcessPoolExecutor, as_completed

from .methylation import apply_methylation_to_sequence, find_cpg_positions


def score_full_meth(seq: str, nb) -> tuple[float, float, float]:
    """Return (G_unmeth, G_meth_full, ddG) for one sequence.

    ``nb`` is a ``NucleosomeBreathModular`` instance. Methylation pattern
    is every CpG dinucleotide in ``seq``.
    """
    cpgs = set(find_cpg_positions(seq))
    seq_m = apply_methylation_to_sequence(seq, cpgs)

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
        "n_cpg": len(positions),
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

    For ``'stable'`` / ``'unstable'``: CpG positions are sampled from
    ``pos_map['stab']`` / ``pos_map['destab']`` respectively, with
    ``n_cpg`` drawn uniformly from ``cpg_count_range`` and respecting
    non-adjacency (resamples until a valid set is drawn).

    For ``'inert'``: no CpGs are injected; ``n_candidates`` bare-backbone
    candidates are scored to capture backbone-to-backbone variation.
    """
    rng = np.random.default_rng(seed)

    backbones = [widom_backbone(seed=int(s)) for s in rng.integers(1, 10**9, size=n_backbones)]

    if class_name == "stable":
        position_pool = pos_map["stab"]
    elif class_name == "unstable":
        position_pool = pos_map["destab"]
    elif class_name == "inert":
        position_pool = np.array([], dtype=int)
    else:
        raise ValueError(f"unknown class_name {class_name!r}")

    tasks = []
    for _ in range(n_candidates):
        backbone = backbones[int(rng.integers(0, n_backbones))]
        if class_name == "inert":
            positions: tuple[int, ...] = ()
        else:
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
```

- [ ] **Step 2: Smoke-check `score_full_meth` only**

Run (the worker-pool path is exercised in Task 5; here we just check the single-call scorer):

```bash
python - <<'PY'
from femodules.nucleosome_breath_modular import NucleosomeBreathModular
from src.config.config_loader import load_cgna_config  # adjust if loader name differs
from src.pipeline.sequence_design import widom_backbone, inject_cpgs, score_full_meth

# Build nb from the project's CGNA config; if loader path differs, copy the
# notebook's CGNA_CONFIG construction here.
cfg = load_cgna_config()
nb = NucleosomeBreathModular(cfg)

bb  = widom_backbone()
seq = inject_cpgs(bb, [40, 60, 80])
g_un, g_me, dd = score_full_meth(seq, nb)
print(f"G_un={g_un:.3f}  G_meth={g_me:.3f}  ddG={dd:+.3f} kT")
PY
```

If `load_cgna_config` does not exist with that name, copy the `CGNA_CONFIG = ...` line from the top of `notebooks/Modelvalidation_meth_exp+analysis.ipynb` into this snippet instead. Expected: three finite floats.

- [ ] **Step 3: Commit**

```bash
git add src/pipeline/sequence_design.py
git commit -m "Add score_full_meth and design_library to sequence_design module"
```

---

## Task 4: Create the new notebook with imports + data load + pos_map

Build `notebooks/Sequence_engineering.ipynb` from scratch via nbformat. Initial cells: imports, load persisted scan outputs, build `pos_map`, quick visual confirmation.

**Files:**
- Create: `notebooks/Sequence_engineering.ipynb` (gitignored, no commit)

- [ ] **Step 1: Generate the notebook skeleton + first 4 cells**

Run from the repo root:

```bash
python - <<'PY'
import nbformat, pathlib

nb = nbformat.v4.new_notebook()
nb.cells = []

# ── markdown: title ─────────────────────────────────────────────────
nb.cells.append(nbformat.v4.new_markdown_cell(
    "# Sequence engineering: targeted methylation response\n\n"
    "Designs three libraries (stable / unstable / inert under full CpG methylation)\n"
    "using the rule-based approach from\n"
    "`docs/superpowers/specs/2026-05-26-sequence-engineering-design.md`.\n"
    "Requires the persisted scan outputs at `data/random_scan/df.csv` and\n"
    "`data/random_scan/sequences.json` (Task 1)."
))

# ── code: imports + data load ───────────────────────────────────────
nb.cells.append(nbformat.v4.new_code_cell('''
import json, pathlib, sys
sys.path.insert(0, str(pathlib.Path.cwd()))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from femodules.nucleosome_breath_modular import NucleosomeBreathModular
from src.pipeline.sequence_design import (
    position_map, widom_backbone, inject_cpgs, score_full_meth, design_library,
    WINDOW_SIZE, DYAD,
)

# Reuse the same CGNA_CONFIG construction as Modelvalidation_meth_exp+analysis.ipynb.
# Paste here the CGNA_CONFIG = ... line(s) from that notebook before running.
from src.config.config_loader import load_cgna_config  # adjust if needed
CGNA_CONFIG = load_cgna_config()
N_WORKERS   = 8

scan_dir = pathlib.Path("data/random_scan")
df        = pd.read_csv(scan_dir / "df.csv")
sequences = json.loads((scan_dir / "sequences.json").read_text())
print(f"df: {len(df)} rows  sequences: {len(sequences)}  (WINDOW_SIZE={WINDOW_SIZE}, DYAD={DYAD})")
'''))

# ── code: build pos_map and visualize ──────────────────────────────
nb.cells.append(nbformat.v4.new_code_cell('''
pos_map = position_map(df)
print(f"|stab|   = {len(pos_map['stab'])}  positions")
print(f"|destab| = {len(pos_map['destab'])} positions")

fig, ax = plt.subplots(figsize=(10, 3))
positions = np.arange(WINDOW_SIZE - 1)
ax.plot(positions, pos_map["mean_ddg"], color="steelblue", lw=1)
ax.scatter(pos_map["stab"],   pos_map["mean_ddg"][pos_map["stab"]],   color="seagreen", s=20, label="stab pool", zorder=3)
ax.scatter(pos_map["destab"], pos_map["mean_ddg"][pos_map["destab"]], color="crimson",  s=20, label="destab pool", zorder=3)
ax.axhline(0, color="k", lw=0.5, ls="--")
ax.axvline(DYAD, color="k", lw=0.5, ls=":")
ax.set_xlabel("CpG position (bp)"); ax.set_ylabel("mean ΔΔG (kT)")
ax.legend(loc="upper right", fontsize=8)
ax.set_title("Per-position mean ΔΔG from random-sequence scan")
plt.tight_layout(); plt.show()
'''))

# ── code: smoke-check backbone + injection ─────────────────────────
nb.cells.append(nbformat.v4.new_code_cell('''
bb = widom_backbone()
print(f"backbone[68:80] = {bb[68:80]}   (expect GGCCC centered at 73)")
print(f"backbone has CG? {('CG' in bb)}")

example = inject_cpgs(bb, pos_map["stab"][:5])
print(f"example with 5 stabilizing CpGs: ...{example[40:80]}...")
print(f"  CpG count = {example.count('CG')}")
'''))

out = pathlib.Path("notebooks/Sequence_engineering.ipynb")
nbformat.write(nb, out)
print(f"wrote {out}  ({len(nb.cells)} cells)")
PY
```

- [ ] **Step 2: Open the notebook and run the four cells**

In JupyterLab/VS Code, open `notebooks/Sequence_engineering.ipynb` and Run All. Expected:
- Cell 2 prints `df: ~784 rows  sequences: 50 sequences  ...`
- Cell 3 plots mean ΔΔG vs position with green and red markers; prints `|stab| ≈ 36` and `|destab| ≈ 36` (roughly the bottom/top quartile of finite positions).
- Cell 4 prints `backbone[68:80] = ...GGCCC...` containing `GGCCC` centered at position 73, `backbone has CG? False`, and the example shows 5 CGs.

If `load_cgna_config` is not the correct name, copy the existing notebook's `CGNA_CONFIG = ...` construction into Cell 2 in place of the loader import. No commit (notebook is gitignored).

---

## Task 5: Run the three libraries + baseline, write `df_lib`

Generate 100 candidates × 3 classes and also score the 50 random sequences as a baseline (full-methylation ΔΔG).

**Files:**
- Modify: `notebooks/Sequence_engineering.ipynb` (append two cells)

- [ ] **Step 1: Append the library + baseline cell**

Run:

```bash
python - <<'PY'
import nbformat, pathlib
nb_path = pathlib.Path("notebooks/Sequence_engineering.ipynb")
nb = nbformat.read(nb_path, as_version=4)

nb.cells.append(nbformat.v4.new_code_cell('''
# ── Run all three libraries ────────────────────────────────────────
import time
t0 = time.time()
df_stable   = design_library("stable",   pos_map, CGNA_CONFIG, n_candidates=100,
                             cpg_count_range=(4, 12), n_workers=N_WORKERS, seed=1)
df_unstable = design_library("unstable", pos_map, CGNA_CONFIG, n_candidates=100,
                             cpg_count_range=(4, 12), n_workers=N_WORKERS, seed=2)
df_inert    = design_library("inert",    pos_map, CGNA_CONFIG, n_candidates=100,
                             n_workers=N_WORKERS, seed=3)
df_lib = pd.concat([df_stable, df_unstable, df_inert], ignore_index=True)
print(f"df_lib: {len(df_lib)} rows  ({time.time()-t0:.1f}s)")
df_lib.groupby("class")["ddG"].agg(["mean", "std", "min", "max"])
'''))

nb.cells.append(nbformat.v4.new_code_cell('''
# ── Baseline: full-methylation ΔΔG for the 50 random scan sequences ──
from concurrent.futures import ProcessPoolExecutor, as_completed
from src.pipeline.sequence_design import _build_candidate

baseline_tasks = [(s["sequence"], "baseline_random", tuple(), CGNA_CONFIG) for s in sequences]
rows = []
with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
    futures = [pool.submit(_build_candidate, t) for t in baseline_tasks]
    for fut in as_completed(futures):
        rows.append(fut.result())
df_baseline = pd.DataFrame(rows)
# Override n_cpg / cpg_positions to reflect what is actually in the sequence
from src.pipeline.methylation import find_cpg_positions
df_baseline["n_cpg"]         = df_baseline["seq"].map(lambda s: len(find_cpg_positions(s)))
df_baseline["cpg_positions"] = df_baseline["seq"].map(lambda s: tuple(find_cpg_positions(s)))
print(f"baseline: {len(df_baseline)} rows")
df_baseline[["n_cpg", "G_un", "G_meth", "ddG"]].describe()
'''))

nbformat.write(nb, nb_path)
print("appended 2 cells")
PY
```

- [ ] **Step 2: Run both new cells**

In the notebook, run the two new cells. Expected for the first: 300 rows, runtime in the order of minutes (depending on `N_WORKERS`), and the per-class `ddG` summary shows `stable.mean < 0`, `unstable.mean > 0`, `inert.mean ≈ 0`. Expected for the second: 50 rows, `ddG.describe()` printed.

If `stable.mean` is not clearly negative (or `unstable.mean` not clearly positive), the additivity assumption is weaker than hoped — flag it and continue; the plot in Task 6 still shows the actual outcome.

No commit (notebook gitignored).

---

## Task 6: Plot class separation + print top-10 picks per class

**Files:**
- Modify: `notebooks/Sequence_engineering.ipynb` (append two cells)

- [ ] **Step 1: Append plotting + top-picks cells**

Run:

```bash
python - <<'PY'
import nbformat, pathlib
nb_path = pathlib.Path("notebooks/Sequence_engineering.ipynb")
nb = nbformat.read(nb_path, as_version=4)

nb.cells.append(nbformat.v4.new_code_cell('''
# ── Class-separation histogram + random baseline overlay ───────────
fig, ax = plt.subplots(figsize=(9, 5))
bins = np.linspace(
    min(df_lib["ddG"].min(), df_baseline["ddG"].min()) - 0.5,
    max(df_lib["ddG"].max(), df_baseline["ddG"].max()) + 0.5,
    50,
)
palette = {"stable": "seagreen", "unstable": "crimson", "inert": "slategray"}
for cls, color in palette.items():
    ax.hist(df_lib[df_lib["class"] == cls]["ddG"], bins=bins,
            alpha=0.55, color=color, label=f"{cls} (N={(df_lib['class']==cls).sum()})")
ax.hist(df_baseline["ddG"], bins=bins, histtype="step",
        color="black", lw=1.5, label=f"random baseline (N={len(df_baseline)})")
ax.axvline(0, color="k", lw=0.6, ls="--")
ax.set_xlabel("ΔΔG = G(full-meth) − G(unmeth)  (kT)")
ax.set_ylabel("count")
ax.legend()
ax.set_title("Methylation-response distribution per designed class")
plt.tight_layout(); plt.show()
'''))

nb.cells.append(nbformat.v4.new_code_cell('''
# ── Top-10 picks per class ─────────────────────────────────────────
cols = ["class", "n_cpg", "G_un", "G_meth", "ddG", "seq"]
print("=== TOP 10 stable_on_meth (most negative ΔΔG) ===")
print(df_lib[df_lib["class"]=="stable"].nsmallest(10, "ddG")[cols].to_string(index=False))
print("\\n=== TOP 10 unstable_on_meth (most positive ΔΔG) ===")
print(df_lib[df_lib["class"]=="unstable"].nlargest(10, "ddG")[cols].to_string(index=False))
print("\\n=== TOP 10 inert_to_meth (smallest |ΔΔG|) ===")
inert = df_lib[df_lib["class"]=="inert"].copy()
inert["abs_ddG"] = inert["ddG"].abs()
print(inert.nsmallest(10, "abs_ddG")[cols].to_string(index=False))
'''))

nbformat.write(nb, nb_path)
print("appended 2 cells")
PY
```

- [ ] **Step 2: Run both cells**

In the notebook, run the new cells. Expected:
- Histogram shows three colored bars (green = stable shifted left, red = unstable shifted right, gray = inert near zero) with the black step curve (random baseline) typically straddling zero.
- Top-10 tables print three blocks, each with 10 rows including the full 147 bp sequence.

No commit (notebook gitignored).

---

## Final verification

- [ ] In a fresh kernel, run the existing scan notebook through to the save cell (Task 1 Step 2) to ensure `data/random_scan/{df.csv,sequences.json}` are present.
- [ ] In a fresh kernel of `notebooks/Sequence_engineering.ipynb`, Run All. All cells execute; the class-separation figure renders.
- [ ] `git status` shows only the two module-related commits on `bp-level-prop` (sequence_design.py and __init__.py), nothing tracked from notebooks/ or data/.

---

## Self-review

- **Spec coverage**: position_map (T2) + widom_backbone (T2) + inject_cpgs (T2) + score_full_meth (T3) + design_library (T3) + notebook with the 5-stage pipeline (T4–T6) + class-separation plot + top-10 picks (T6) — all spec sections covered.
- **Placeholders**: none. Every code step has full code; the only branch is the `load_cgna_config` import which is explicitly noted to be swapped for the notebook's existing `CGNA_CONFIG = ...` block if the loader name differs.
- **Type consistency**: column names (`ddG`, `cpg_pos`, `G_un`, `G_meth`, `n_cpg`, `cpg_positions`, `seq`, `class`) are used identically across `position_map`, `_build_candidate`, `design_library`, plotting, and top-picks; `pos_map` keys (`stab`, `destab`, `mean_ddg`) are used consistently.
- **Out-of-scope items deliberately omitted** (per the spec): graph-based optimization, iterative refinement, partial methylation, positional CpG-density sanity panel, diversity metrics.
