# Sequence Engineering for Targeted Methylation Response — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a library of 147 bp DNA sequences in three classes (stable / unstable / inert under full CpG methylation) using a Widom-style CpG-free backbone plus rule-driven CpG injection, and validate by computing ΔΔG = G_meth_full − G_unmeth for each candidate.

**Architecture:** A new module `src/pipeline/sequence_design.py` holds five focused functions (`position_map`, `widom_backbone`, `inject_cpgs`, `score_full_meth`, `design_library`). A new notebook `notebooks/Sequence_engineering.ipynb` loads the pre-computed single-site scan CSV, runs the library design + scoring (parallel), and plots class-separation.

**Tech Stack:** `NucleosomeBreathModular` (CGNA+ free-energy soft-binding), `apply_methylation_to_sequence`, `find_cpg_positions`, `ProcessPoolExecutor`, pandas, matplotlib, nbformat.

**Branch:** `bp-level-prop`. Module commits go to git; notebooks and data files are gitignored.

---

## Conventions

- **ΔΔG sign**: `ddG = G_meth_full − G_unmeth`. Negative → stabilizing, positive → destabilizing.
- **Free-energy quantity**: `dF_el = result.F − result.F_freedna` (binding free energy minus free-DNA reference), matching the existing notebook's convention.
- **Binding spec**: `left=0, right=13, style="b_index"` (fully wrapped nucleosome).
- **Window**: `L = 147`, `dyad = 73`.
- **Input data**: `notebooks/files/random_seq_single_site_ddG.csv` — pre-computed single-site scan, columns include `seq_label`, `cpg_pos`, `ddG`.

---

## Task 1: Create `sequence_design` module — utility functions

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

Expected: `backbone[68:80]` contains `GGCCC` centered around position 73, jitter changes ~7-ish positions, inject_cpgs prints OK.

- [ ] **Step 4: Commit**

```bash
git add src/pipeline/sequence_design.py src/pipeline/__init__.py
git commit -m "Add sequence_design module utilities (position_map, widom_backbone, inject_cpgs)"
```

---

## Task 2: Add scoring + library functions to the module

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
```

- [ ] **Step 2: Smoke-check `score_full_meth` only**

Run (the worker-pool path is exercised in Task 4; here we just check the single-call scorer):

```bash
python - <<'PY'
from femodules.nucleosome_breath_modular import NucleosomeBreathModular
from src.pipeline.sequence_design import widom_backbone, inject_cpgs, score_full_meth

# Build nb from the project's CGNA config. Copy the CGNA_CONFIG = ...
# construction here from the top of
# notebooks/Modelvalidation_meth_exp+analysis.ipynb (the line(s) that
# build the CgnaConfig dataclass before NucleosomeBreathModular is
# instantiated).
from src.config.custom_types import CgnaConfig
CGNA_CONFIG = CgnaConfig()  # replace with the actual fields used in the notebook
nb = NucleosomeBreathModular(CGNA_CONFIG)

bb  = widom_backbone()
seq = inject_cpgs(bb, [40, 60, 80])
g_un, g_me, dd = score_full_meth(seq, nb)
print(f"G_un={g_un:.3f}  G_meth={g_me:.3f}  ddG={dd:+.3f} kT")
PY
```

Expected: three finite floats. (If `CgnaConfig()` with no args fails, copy the exact construction from the existing notebook.)

- [ ] **Step 3: Commit**

```bash
git add src/pipeline/sequence_design.py
git commit -m "Add score_full_meth and design_library to sequence_design module"
```

---

## Task 3: Create the new notebook — references, imports, data load, pos_map

Build `notebooks/Sequence_engineering.ipynb` from scratch via nbformat. Cells: references / motivation markdown, imports + CSV load, build `pos_map` and plot it, smoke-check backbone + injection.

**Files:**
- Create: `notebooks/Sequence_engineering.ipynb` (gitignored, no commit)

- [ ] **Step 1: Generate the notebook skeleton + first 5 cells**

Run from the repo root:

```bash
python - <<'PY'
import nbformat, pathlib

nb = nbformat.v4.new_notebook()
nb.cells = []

# ── markdown: title + references ───────────────────────────────────
nb.cells.append(nbformat.v4.new_markdown_cell(
    "# Sequence engineering: targeted methylation response\n\n"
    "Designs three libraries (stable / unstable / inert under full CpG methylation)\n"
    "using the rule-based approach from\n"
    "`docs/superpowers/specs/2026-05-26-sequence-engineering-design.md`.\n"
    "Single-site reference data is loaded from\n"
    "`notebooks/files/random_seq_single_site_ddG.csv`.\n\n"
    "## References for the Widom-style backbone\n\n"
    "The CpG-free 10-bp tile `GGCCCAATTT` (GC-rich half centered on the\n"
    "dyad, AT-rich half centered on SHL ±0.5) reflects the following\n"
    "nucleosome-positioning rules:\n\n"
    "1. **Zuiddam & Schiessel, Phys. Rev. E 99, 012422 (2019).**\n"
    "   Shortest-path optimization over the trinucleotide nucleosome\n"
    "   energy model. Lowest-energy 147 bp sequences are ~80 % GC with\n"
    "   CC/GG/GC steps where the major groove faces the histone (integer\n"
    "   SHL positions); highest-energy sequences are A/T-rich with\n"
    "   A-tracts in minor-groove-in positions.\n"
    "2. **Pérez-Pulido et al., Nucleic Acids Res. 50, 1864 (2022).**\n"
    "   DNA methylation increases roll and decreases twist at CpG steps;\n"
    "   methylated CpGs prefer minor-groove-in positions whereas\n"
    "   unmethylated CpGs prefer minor-groove-out.\n"
    "3. **Collings et al., Sci. Reports 3, 2121 (2013).**\n"
    "   Ten-bp periodicity of (m)CpG dinucleotides along nucleosomal\n"
    "   DNA, dictating the rotational frame of methylation effects.\n\n"
    "These rules motivate (a) the GC-vs-AT half split of the tile\n"
    "and (b) the periodic CpG-free framing — so methylation response is\n"
    "controlled by where CpGs are injected, not by the backbone itself."
))

# ── code: imports + data load ───────────────────────────────────────
nb.cells.append(nbformat.v4.new_code_cell('''
import pathlib, sys
sys.path.insert(0, str(pathlib.Path.cwd()))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from femodules.nucleosome_breath_modular import NucleosomeBreathModular
from src.pipeline.sequence_design import (
    position_map, widom_backbone, inject_cpgs, score_full_meth, design_library,
    WINDOW_SIZE, DYAD,
)

# Reuse the same CGNA_CONFIG construction as
# notebooks/Modelvalidation_meth_exp+analysis.ipynb. Paste the
# `CGNA_CONFIG = ...` line(s) from that notebook here before running.
from src.config.custom_types import CgnaConfig
CGNA_CONFIG = CgnaConfig()  # <-- replace with the exact construction used in the existing notebook
N_WORKERS   = 8

df = pd.read_csv("notebooks/files/random_seq_single_site_ddG.csv")
print(f"loaded df: {len(df)} rows, columns: {list(df.columns)}")
df.head(3)
'''))

# ── code: build pos_map and visualize ──────────────────────────────
nb.cells.append(nbformat.v4.new_code_cell('''
pos_map = position_map(df)
print(f"|stab|    = {len(pos_map['stab'])}  positions  (mean ΔΔG most negative)")
print(f"|destab|  = {len(pos_map['destab'])} positions  (mean ΔΔG most positive)")
print(f"|neutral| = {len(pos_map['neutral'])} positions  (|mean ΔΔG| smallest)")

fig, ax = plt.subplots(figsize=(10, 3))
positions = np.arange(WINDOW_SIZE - 1)
ax.plot(positions, pos_map["mean_ddg"], color="steelblue", lw=1)
ax.scatter(pos_map["stab"],    pos_map["mean_ddg"][pos_map["stab"]],    color="seagreen", s=20, label="stab",    zorder=3)
ax.scatter(pos_map["destab"],  pos_map["mean_ddg"][pos_map["destab"]],  color="crimson",  s=20, label="destab",  zorder=3)
ax.scatter(pos_map["neutral"], pos_map["mean_ddg"][pos_map["neutral"]], color="slategray",s=20, label="neutral", zorder=3)
ax.axhline(0, color="k", lw=0.5, ls="--")
ax.axvline(DYAD, color="k", lw=0.5, ls=":")
ax.set_xlabel("CpG position (bp)"); ax.set_ylabel("mean ΔΔG (kT)")
ax.legend(loc="upper right", fontsize=8)
ax.set_title("Per-position mean ΔΔG from single-site scan (3 design pools)")
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
- Cell 2 prints `loaded df: 4175 rows, columns: [...]`.
- Cell 3 plots mean ΔΔG vs position with green / red / gray markers; prints `|stab|`, `|destab|`, `|neutral|` (roughly a quarter of finite positions each).
- Cell 4 prints `backbone[68:80] = ...GGCCC...` containing `GGCCC` centered at position 73, `backbone has CG? False`, and the example shows 5 CGs.

If `CgnaConfig()` with no args fails, paste the exact `CGNA_CONFIG = ...` construction from `notebooks/Modelvalidation_meth_exp+analysis.ipynb` into Cell 2. No commit (notebook is gitignored).

---

## Task 4: Run the three libraries + a random baseline

Generate 100 candidates × 3 classes, plus a 50-sequence random baseline scored under full methylation.

**Files:**
- Modify: `notebooks/Sequence_engineering.ipynb` (append two cells)

- [ ] **Step 1: Append the library + baseline cells**

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
                             cpg_count_range=(4, 12), n_workers=N_WORKERS, seed=3)
df_lib = pd.concat([df_stable, df_unstable, df_inert], ignore_index=True)
print(f"df_lib: {len(df_lib)} rows  ({time.time()-t0:.1f}s)")
df_lib.groupby("class")["ddG"].agg(["mean", "std", "min", "max"])
'''))

nb.cells.append(nbformat.v4.new_code_cell('''
# ── Baseline: 50 fresh random sequences scored under full methylation ──
from concurrent.futures import ProcessPoolExecutor, as_completed
from src.pipeline.sequence_design import _build_candidate

rng_base = np.random.default_rng(7)
random_seqs = ["".join(rng_base.choice(list("ACGT"), size=WINDOW_SIZE)) for _ in range(50)]
baseline_tasks = [(s, "baseline_random", tuple(), CGNA_CONFIG) for s in random_seqs]

rows = []
with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
    for fut in as_completed(pool.submit(_build_candidate, t) for t in baseline_tasks):
        rows.append(fut.result())
df_baseline = pd.DataFrame(rows)
print(f"baseline: {len(df_baseline)} rows")
df_baseline[["n_cpg", "G_un", "G_meth", "ddG"]].describe()
'''))

nbformat.write(nb, nb_path)
print("appended 2 cells")
PY
```

- [ ] **Step 2: Run both new cells**

Expected for the first cell: 300 rows; runtime in the order of minutes (depending on `N_WORKERS`); per-class `ddG` summary shows `stable.mean < 0`, `unstable.mean > 0`, `inert.mean ≈ 0`. Expected for the second cell: 50 rows with `ddG.describe()` printed.

No commit (notebook gitignored).

---

## Task 5: Plot class separation + print top-10 picks per class

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

Expected: histogram with three colored bars (green = stable shifted left, red = unstable shifted right, gray = inert near zero) plus the black step curve (random baseline) typically straddling zero; three top-10 tables print, each with 10 rows including the 147 bp sequence.

No commit (notebook gitignored).

---

## Final verification

- [ ] In a fresh kernel of `notebooks/Sequence_engineering.ipynb`, Run All. All cells execute; the class-separation figure renders.
- [ ] `git status` shows only the two module-related commits on `bp-level-prop` (sequence_design.py and __init__.py), nothing tracked from notebooks/ or data/.

---

## Self-review

- **Spec coverage**: position_map with three pools (T1) + widom_backbone (T1) + inject_cpgs (T1) + score_full_meth (T2) + design_library handling stable/unstable/inert via neutral pool (T2) + notebook with references + data load + pos_map + libraries + baseline + plot + top picks (T3–T5) — every spec section covered.
- **Placeholders**: none. Every code step has full code; the only branch is the `CGNA_CONFIG = CgnaConfig()` line, explicitly flagged to be swapped for the notebook's existing construction if defaults don't match.
- **Type consistency**: column names (`ddG`, `cpg_pos`, `G_un`, `G_meth`, `n_cpg`, `cpg_positions`, `seq`, `class`) used identically across `position_map`, `_build_candidate`, `design_library`, plotting, and top-picks; `pos_map` keys (`stab`, `destab`, `neutral`, `mean_ddg`) used consistently.
- **Out-of-scope items deliberately omitted** (per the spec): graph-based optimization, iterative refinement, partial methylation patterns, positional CpG-density sanity panel, diversity metrics.
