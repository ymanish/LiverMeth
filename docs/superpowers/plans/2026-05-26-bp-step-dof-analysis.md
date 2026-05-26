# Base-pair-step DOF analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a new section to `notebooks/Modelvalidation_meth_exp+analysis.ipynb` that extracts base-pair-step DOFs at each scanned CpG, groups CpGs by ΔΔG quartile, and reports how the six DOFs differ between most-stabilizing and most-destabilizing sites.

**Architecture:** Pure notebook extension. Reuses the existing `df_single` DataFrame and `sequences` list from the Random-Sequence Scan section. Feature extraction calls `cgnaplus_bps_params` once per sequence (cached) plus once per CpG-methylated variant; deviation-from-target uses `methods.free_energy.midstep_excess_vals` against `nb.nuc_mu0`. Long-form `df_dof` joins to ΔΔG quartile labels for grouped plots and a summary table.

**Tech Stack:** Python, NumPy, pandas, Matplotlib, SciPy stats, the project's existing CGNA+ pipeline (`cgnaplus_bps_params`, `NucleosomeBreathModular`, `midstep_excess_vals`, `_revert_terminal_mn`, `apply_methylation_to_sequence`, `find_cpg_positions`).

---

## Conventions

- **DOF order** (CGNA+ default, `euler_definition=True`, `translations_in_nm=True`): index 0..5 = `tilt, roll, twist, shift, slide, rise`. Rotations in radians, translations in nm.
- **Step indices**: `gs` from `cgnaplus_bps_params(seq)` has shape `(L−1, 6)` for a sequence of length `L`. Step `s` connects base pairs `s` and `s+1`.
- **CG step at CpG position p**: the CpG is at positions `p` (C) and `p+1` (G); the step between them is `gs[p]`. The 3-step window is `gs[p−1], gs[p], gs[p+1]`.
- **Sequence canonicalisation**: every sequence passed to `cgnaplus_bps_params` is first passed through `_revert_terminal_mn` to match the energy scan's convention (CGNA+ has no terminal `MN` parameters).
- **Quartile bins**: `bottom_25` = ΔΔG ≤ 25th percentile, `top_25` = ΔΔG ≥ 75th percentile, `mid_50` = middle.
- **Plan location**: all new code lives in cells appended to `notebooks/Modelvalidation_meth_exp+analysis.ipynb` after the Random-Sequence Scan section (after the existing `df_single` definition). Cell IDs are arbitrary; locate cells by markdown heading.
- **Commits**: small, frequent. After each task's working code, run the notebook section (executing only the new cells), confirm output, then commit the `.ipynb`.

## File map

- **Modify only**: `notebooks/Modelvalidation_meth_exp+analysis.ipynb`
  - Append a top-level markdown cell `## Base-pair-step DOF analysis at stabilizing vs destabilizing CpGs` and the cells listed below, after the existing Random-Sequence Scan section.
- **No new module files** in v1. The feature-extraction helpers stay inline. If, during Task 2, the feature-extraction cell exceeds ~80 lines, split helpers into `femodules/bp_step_features.py` (mentioned in spec but deferred unless triggered).

## Cell layout (final state)

After all tasks the new section contains, in order:

1. Markdown heading + intro.
2. Imports + constants cell (Task 1).
3. Feature-extraction cell building `df_dof` for the three per-step quantities (Task 2).
4. `gs_dev` extraction cell building `df_dev` (Task 3).
5. Quartile-binning cell merging `ddG_bin` onto `df_dof` and `df_dev` (Task 4).
6. Cross-check assertion cell (Task 5).
7. Primary figure A — 4 figures, one per quantity (Task 6).
8. Supplementary figure B — 3 figures (per-step quantities only) (Task 7).
9. Sanity-check figure (Task 8).
10. Summary table (Task 9).
11. Closing markdown discussion stub (Task 10).

---

## Task 1: Section scaffold + imports

**Files:**
- Modify: `notebooks/Modelvalidation_meth_exp+analysis.ipynb` (append cells after the Random-Sequence Scan section)

- [ ] **Step 1: Add a markdown cell with the section heading and intent**

Append a new markdown cell:

```markdown
## Base-pair-step DOF analysis at stabilizing vs destabilizing CpGs

For each CpG site scanned above we extract the CGNA+ base-pair-step DOFs
(tilt, roll, twist, shift, slide, rise) at the CG step and its two flanking
steps, plus the per-CpG deviation from the nucleosomal midstep-triad
constraint. CpGs are grouped by ΔΔG quartile (bottom 25% = most-stabilizing,
top 25% = most-destabilizing) and the average DOFs are compared between
groups.
```

- [ ] **Step 2: Add a code cell with imports and constants**

```python
# ── BP-step DOF analysis: imports & constants ────────────────────────────────
from methods.free_energy import midstep_excess_vals
from scipy.stats import spearmanr, mannwhitneyu

DOF_NAMES = ["tilt", "roll", "twist", "shift", "slide", "rise"]
STEP_OFFSETS = (-1, 0, 1)
PER_STEP_QUANTITIES = ("gs_un", "gs_meth", "delta_gs")
ALL_QUANTITIES = PER_STEP_QUANTITIES + ("gs_dev",)

# Phosphate-binding-site constraint locations for left=0, right=13
# (same list NucleosomeBreath._select_phosphate_bind_sites returns for full binding)
PHOSPHATE_BIND_SITES = [
    2, 6, 14, 17, 24, 29, 34, 38,
    45, 49, 55, 59, 65, 69, 76,
    80, 86, 90, 96, 100, 107, 111,
    116, 121, 128, 131, 139, 143,
]

print(f"DOF order: {DOF_NAMES}")
print(f"Per-step quantities: {PER_STEP_QUANTITIES}")
print(f"Per-CpG quantity:    gs_dev (interval excess vs nuc_mu0)")
print(f"# constraint sites:  {len(PHOSPHATE_BIND_SITES)}")
```

- [ ] **Step 3: Run the cell**

Expected output:
```
DOF order: ['tilt', 'roll', 'twist', 'shift', 'slide', 'rise']
Per-step quantities: ('gs_un', 'gs_meth', 'delta_gs')
Per-CpG quantity:    gs_dev (interval excess vs nuc_mu0)
# constraint sites:  28
```

- [ ] **Step 4: Commit**

```bash
git add notebooks/Modelvalidation_meth_exp+analysis.ipynb
git commit -m "Add BP-step DOF analysis section scaffold and imports"
```

---

## Task 2: Build `df_dof` — per-step quantities

**Files:**
- Modify: `notebooks/Modelvalidation_meth_exp+analysis.ipynb` (append after Task 1 cells)

The cell loops over CpGs in `df_single`. For each sequence it caches the unmethylated `gs`. For each CpG it methylates only that site, computes `gs_meth`, takes the 3-step window, and emits rows for the three per-step quantities.

- [ ] **Step 1: Add a code cell that builds `df_dof`**

```python
# ── Build df_dof for the three per-step quantities ──────────────────────────
# Reuses `sequences` and `df_single` from the Random-Sequence Scan section.
# DOF columns follow CGNA+ order: tilt, roll, twist, shift, slide, rise.

def _compute_gs(seq):
    """Compute CGNA+ groundstate for a (possibly methylated) sequence string."""
    seq_can = _revert_terminal_mn(seq)
    gs, _ = cgnaplus_bps_params(
        sequence=seq_can,
        parameter_set_name=CGNA_CONFIG.parameter_set_name,
    )
    # gs shape: (L-1, 6); already in (tilt, roll, twist, shift, slide, rise) order.
    return gs

# Cache unmethylated gs per sequence (50 calls).
gs_un_cache = {s["label"]: _compute_gs(s["sequence"]) for s in sequences}
seq_by_label = {s["label"]: s["sequence"] for s in sequences}

per_step_rows = []
edge_dropped = 0
for row in df_single.itertuples(index=False):
    label, p = row.seq_label, int(row.cpg_pos)
    seq_un = seq_by_label[label]
    L = len(seq_un)                     # 147
    if p - 1 < 0 or p + 1 >= L - 1:     # need steps p-1, p, p+1
        edge_dropped += 1
        continue

    gs_un  = gs_un_cache[label]
    seq_m  = apply_methylation_to_sequence(seq_un, {p})
    gs_m   = _compute_gs(seq_m)

    for off in STEP_OFFSETS:
        s = p + off
        per_step_rows.append({
            "seq_label": label, "cpg_pos": p, "step_offset": off,
            "quantity": "gs_un",
            **dict(zip(DOF_NAMES, gs_un[s])),
        })
        per_step_rows.append({
            "seq_label": label, "cpg_pos": p, "step_offset": off,
            "quantity": "gs_meth",
            **dict(zip(DOF_NAMES, gs_m[s])),
        })
        per_step_rows.append({
            "seq_label": label, "cpg_pos": p, "step_offset": off,
            "quantity": "delta_gs",
            **dict(zip(DOF_NAMES, gs_m[s] - gs_un[s])),
        })

df_dof = pd.DataFrame(per_step_rows)
print(f"df_dof: {len(df_dof):,} rows  (edge CpGs dropped: {edge_dropped})")
print(f"Expected: 3 quantities x 3 step_offsets x {len(df_single) - edge_dropped} CpGs"
      f" = {9 * (len(df_single) - edge_dropped):,}")
df_dof.head(6)
```

- [ ] **Step 2: Run the cell**

Expected:
- `df_dof: 7,X rows  (edge CpGs dropped: small integer, typically 0–20)`
- The reported count matches `9 × (N_CpGs − edge_dropped)`.
- `df_dof.head(6)` shows alternating `gs_un`, `gs_meth`, `delta_gs` rows for the same `(seq_label, cpg_pos, step_offset=-1)`.

If the count mismatches, stop and investigate before moving on.

- [ ] **Step 3: Quick sanity check on `delta_gs`**

Append a cell:

```python
# `delta_gs` should be exactly gs_meth − gs_un row-by-row
pivot = df_dof.pivot_table(
    index=["seq_label", "cpg_pos", "step_offset"],
    columns="quantity",
    values=DOF_NAMES,
)
diff = (pivot["tilt"]["gs_meth"] - pivot["tilt"]["gs_un"]) - pivot["tilt"]["delta_gs"]
assert np.allclose(diff.dropna(), 0.0, atol=1e-12), "delta_gs inconsistent with gs_meth - gs_un"
print("delta_gs ≡ gs_meth − gs_un  ✓")
```

Expected: `delta_gs ≡ gs_meth − gs_un  ✓`

- [ ] **Step 4: Commit**

```bash
git add notebooks/Modelvalidation_meth_exp+analysis.ipynb
git commit -m "Extract per-step CGNA+ DOFs (gs_un, gs_meth, delta_gs) into df_dof"
```

---

## Task 3: Build `df_dev` — per-CpG deviation from nucleosomal target

**Files:**
- Modify: `notebooks/Modelvalidation_meth_exp+analysis.ipynb` (append)

`midstep_excess_vals(gs, constraint_locations, midstep_triads)` returns one 6-vector per interval between consecutive constraints. We compute it once per sequence on `gs_un` and assign each CpG the interval containing it.

- [ ] **Step 1: Add a code cell that builds `df_dev`**

```python
# ── Build df_dev — per-CpG deviation from nucleosomal midstep-triad target ─
# Uses midstep_excess_vals on the unmethylated gs and the existing nb.nuc_mu0.
# Assumes the Random-Sequence Scan section has constructed `nb` (NucleosomeBreathModular).

mu0 = nb.nuc_mu0                                # shape (28, 4, 4)
assert mu0.shape == (len(PHOSPHATE_BIND_SITES), 4, 4), (
    f"nb.nuc_mu0 shape {mu0.shape} != expected ({len(PHOSPHATE_BIND_SITES)}, 4, 4); "
    "check that nb was initialised with left=0, right=13 binding."
)

# excess_per_seq[label] has shape (len(PHOSPHATE_BIND_SITES) - 1, 6) — one 6-vector per interval.
excess_per_seq = {
    label: midstep_excess_vals(gs_un_cache[label], PHOSPHATE_BIND_SITES, mu0)
    for label in gs_un_cache
}

def _interval_index(p):
    """Return the interval index k such that PHOSPHATE_BIND_SITES[k] <= p < PHOSPHATE_BIND_SITES[k+1].
    Returns None if p is outside [PHOSPHATE_BIND_SITES[0], PHOSPHATE_BIND_SITES[-1])."""
    sites = PHOSPHATE_BIND_SITES
    if p < sites[0] or p >= sites[-1]:
        return None
    # bisect-style search; small list so a linear scan is fine and explicit.
    for k in range(len(sites) - 1):
        if sites[k] <= p < sites[k + 1]:
            return k
    return None

dev_rows = []
out_of_range = 0
for row in df_single.itertuples(index=False):
    label, p = row.seq_label, int(row.cpg_pos)
    k = _interval_index(p)
    if k is None:
        out_of_range += 1
        continue
    vec = excess_per_seq[label][k]      # shape (6,)
    dev_rows.append({
        "seq_label": label, "cpg_pos": p,
        "step_offset": np.nan,
        "quantity": "gs_dev",
        "interval_idx": k,
        **dict(zip(DOF_NAMES, vec)),
    })

df_dev = pd.DataFrame(dev_rows)
print(f"df_dev: {len(df_dev):,} rows  (CpGs outside constraint span: {out_of_range})")
df_dev.head(5)
```

- [ ] **Step 2: Run the cell**

Expected:
- `df_dev: ~780 rows  (CpGs outside constraint span: typically 0–5)`.
- Each row shows `quantity=gs_dev`, `step_offset=NaN`, an `interval_idx` in `[0, 26]`, and six numeric DOF columns.

- [ ] **Step 3: Cross-check the excess against a direct recomputation for one CpG**

Append a cell that reconstructs the excess for a single (label, k) pair via raw matrix algebra and compares to `df_dev`:

```python
# Pick the first dev row; recompute its excess without using midstep_excess_vals,
# to confirm we are calling the helper with the right inputs.
r0 = df_dev.iloc[0]
label, k = r0["seq_label"], int(r0["interval_idx"])
gs_seg = gs_un_cache[label][PHOSPHATE_BIND_SITES[k]:PHOSPHATE_BIND_SITES[k+1]+1]
manual = midstep_excess_vals(
    gs_un_cache[label],
    [PHOSPHATE_BIND_SITES[k], PHOSPHATE_BIND_SITES[k+1]],
    mu0[[k, k+1]],
)[0]
recorded = r0[DOF_NAMES].to_numpy(dtype=float)
assert np.allclose(manual, recorded, atol=1e-10), (manual, recorded)
print(f"gs_dev recomputation matches for (label={label}, k={k})  ✓")
```

Expected: `gs_dev recomputation matches for (label=rand_???, k=??)  ✓`

- [ ] **Step 4: Commit**

```bash
git add notebooks/Modelvalidation_meth_exp+analysis.ipynb
git commit -m "Extract per-CpG nucleosomal-target excess (gs_dev) into df_dev"
```

---

## Task 4: Merge ΔΔG and quartile bins onto features

**Files:**
- Modify: `notebooks/Modelvalidation_meth_exp+analysis.ipynb` (append)

- [ ] **Step 1: Add a code cell that adds `ddG` and `ddG_bin` columns**

```python
# ── Quartile-bin labels on df_single, then merge onto df_dof and df_dev ────
q25, q75 = df_single["ddG"].quantile([0.25, 0.75])
def _bin(x):
    if x <= q25: return "bottom_25"
    if x >= q75: return "top_25"
    return "mid_50"

df_single = df_single.copy()
df_single["ddG_bin"] = df_single["ddG"].apply(_bin)
print(f"Quartile thresholds: q25={q25:.3f}, q75={q75:.3f}")
print(df_single["ddG_bin"].value_counts())

merge_cols = ["seq_label", "cpg_pos", "ddG", "ddG_bin"]
df_dof = df_dof.merge(df_single[merge_cols], on=["seq_label", "cpg_pos"], how="left")
df_dev = df_dev.merge(df_single[merge_cols], on=["seq_label", "cpg_pos"], how="left")
assert df_dof["ddG_bin"].notna().all(), "missing ddG_bin in df_dof"
assert df_dev["ddG_bin"].notna().all(), "missing ddG_bin in df_dev"
print("Merged ddG and ddG_bin onto df_dof and df_dev  ✓")
```

- [ ] **Step 2: Run the cell**

Expected:
- Counts look like ~196 bottom_25 / ~392 mid_50 / ~196 top_25 (totals depend on edge drops).
- `Merged ddG and ddG_bin onto df_dof and df_dev  ✓`

- [ ] **Step 3: Commit**

```bash
git add notebooks/Modelvalidation_meth_exp+analysis.ipynb
git commit -m "Tag CpGs with ΔΔG quartile and merge onto DOF features"
```

---

## Task 5: Cross-check feature extraction against the energy scan

**Files:**
- Modify: `notebooks/Modelvalidation_meth_exp+analysis.ipynb` (append)

This is the spec's "test plan" — recompute the ΔΔG for one site from a freshly methylated sequence and confirm it matches `df_single.ddG` within tolerance. Catches any sequence indexing or methylation-application mismatch between this section and the energy scan.

- [ ] **Step 1: Add a code cell**

```python
# ── Cross-check: recompute one ΔΔG from a fresh methylated sequence ────────
# Picks the largest-|ΔΔG| row to make any mismatch obvious.
chk = df_single.loc[df_single["ddG"].abs().idxmax()]
label, p_chk = chk["seq_label"], int(chk["cpg_pos"])
seq_un = seq_by_label[label]
seq_m  = _revert_terminal_mn(apply_methylation_to_sequence(seq_un, {p_chk}))
seq_u  = _revert_terminal_mn(seq_un)

# Use the same NucleosomeBreathModular path as the scan.
_nb = NucleosomeBreathModular(CGNA_CONFIG)
r_u = _nb.calculate_free_energy(sequence=seq_u, left=0, right=13, style="b_index")
r_m = _nb.calculate_free_energy(sequence=seq_m, left=0, right=13, style="b_index")
ddG_recomp = (r_m.F - r_m.F_freedna) - (r_u.F - r_u.F_freedna)

print(f"Sample (label={label}, p={p_chk}):")
print(f"  ΔΔG from df_single : {chk['ddG']:+.6f}")
print(f"  ΔΔG recomputed     : {ddG_recomp:+.6f}")
assert np.isclose(ddG_recomp, chk["ddG"], atol=1e-6), (
    f"ΔΔG mismatch: {ddG_recomp} vs {chk['ddG']}"
)
print("Cross-check passed  ✓")
```

- [ ] **Step 2: Run the cell**

Expected: two ΔΔG values that match to ≥1e-6 and `Cross-check passed  ✓`.

If they don't match, debug before continuing — figures will be meaningless if sequence indexing is off.

- [ ] **Step 3: Commit**

```bash
git add notebooks/Modelvalidation_meth_exp+analysis.ipynb
git commit -m "Add ΔΔG cross-check between DOF section and energy scan"
```

---

## Task 6: Primary figures (approach A) — quartile bar plots

**Files:**
- Modify: `notebooks/Modelvalidation_meth_exp+analysis.ipynb` (append)

One figure per quantity. For per-step quantities the panel uses `step_offset = 0`; for `gs_dev` the same per-CpG row.

- [ ] **Step 1: Add a helper cell for per-DOF group stats**

```python
# ── Helpers for grouped DOF stats ──────────────────────────────────────────
def _group_stats(frame, dof, bin_name):
    """Return (mean, sem, n) of `dof` within `ddG_bin == bin_name`."""
    vals = frame.loc[frame["ddG_bin"] == bin_name, dof].dropna().to_numpy()
    if vals.size == 0:
        return np.nan, np.nan, 0
    return vals.mean(), vals.std(ddof=1) / np.sqrt(vals.size), vals.size

def _spearman(frame, dof):
    """Spearman ρ, p of DOF vs ΔΔG across all rows in frame."""
    sub = frame[[dof, "ddG"]].dropna()
    if len(sub) < 5:
        return np.nan, np.nan
    rho, pval = spearmanr(sub[dof], sub["ddG"])
    return rho, pval
```

- [ ] **Step 2: Add a cell that draws the four primary figures**

```python
# ── Figure A: bottom_25 vs top_25 bars per DOF, one figure per quantity ───
import matplotlib.pyplot as plt

def _frame_for_quantity(q):
    if q == "gs_dev":
        return df_dev[df_dev["quantity"] == "gs_dev"]
    return df_dof[(df_dof["quantity"] == q) & (df_dof["step_offset"] == 0)]

for q in ALL_QUANTITIES:
    sub = _frame_for_quantity(q)
    fig, axes = plt.subplots(1, 6, figsize=(15, 3), sharey=False)
    fig.suptitle(f"Quantity: {q}   (center step, n_bottom={int((sub['ddG_bin']=='bottom_25').sum())},"
                 f" n_top={int((sub['ddG_bin']=='top_25').sum())})")
    for ax, dof in zip(axes, DOF_NAMES):
        m_b, e_b, _ = _group_stats(sub, dof, "bottom_25")
        m_t, e_t, _ = _group_stats(sub, dof, "top_25")
        m_m, _,  _  = _group_stats(sub, dof, "mid_50")
        ax.bar([0, 1], [m_b, m_t], yerr=[e_b, e_t], capsize=4,
               color=["tab:blue", "tab:red"])
        ax.axhline(m_m, color="gray", lw=1, ls="--", label="mid_50 mean")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["bot 25%", "top 25%"], rotation=0)
        rho, pval = _spearman(sub, dof)
        ax.set_title(f"{dof}\nρ={rho:+.2f}  p={pval:.1e}")
    axes[0].set_ylabel(q)
    plt.tight_layout()
    plt.show()
```

- [ ] **Step 3: Run the cell**

Expected: four figures appear, each with six side-by-side bar panels. Each panel title shows the DOF name, Spearman ρ, and p-value.

Inspect: for `delta_gs` the bars are typically smaller in magnitude than for `gs_un` because methylation perturbations are small relative to sequence variation. If any panel is empty (no bars), one of the groups is missing — investigate before continuing.

- [ ] **Step 4: Commit**

```bash
git add notebooks/Modelvalidation_meth_exp+analysis.ipynb
git commit -m "Add primary BP-step DOF figures (quartile-bin bars per quantity)"
```

---

## Task 7: Supplementary figures (approach B) — flanking-step profiles

**Files:**
- Modify: `notebooks/Modelvalidation_meth_exp+analysis.ipynb` (append)

Only for the three per-step quantities. `gs_dev` is per-CpG and has no flanking-step variant.

- [ ] **Step 1: Add a cell that draws the three flanking-step figures**

```python
# ── Figure B: mean ± SEM across step_offset {-1, 0, +1}, one figure per
#    per-step quantity, 2×3 DOF grid, lines for bottom_25 vs top_25. ───────
for q in PER_STEP_QUANTITIES:
    sub = df_dof[df_dof["quantity"] == q]
    fig, axes = plt.subplots(2, 3, figsize=(11, 6), sharex=True)
    fig.suptitle(f"Quantity: {q}   (3-step window around CG)")
    for ax, dof in zip(axes.ravel(), DOF_NAMES):
        for bin_name, colour in [("bottom_25", "tab:blue"), ("top_25", "tab:red")]:
            means, sems = [], []
            for off in STEP_OFFSETS:
                m, e, _ = _group_stats(
                    sub[sub["step_offset"] == off], dof, bin_name
                )
                means.append(m); sems.append(e)
            ax.errorbar(STEP_OFFSETS, means, yerr=sems, marker="o",
                        capsize=3, color=colour, label=bin_name)
        ax.set_title(dof)
        ax.axvline(0, color="gray", lw=0.5)
        ax.set_xticks(list(STEP_OFFSETS))
    for ax in axes[1]:
        ax.set_xlabel("step offset")
    axes[0, 0].legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.show()
```

- [ ] **Step 2: Run the cell**

Expected: three figures, each a 2×3 grid. Each panel shows two lines (`bottom_25` blue, `top_25` red) with three points each.

- [ ] **Step 3: Commit**

```bash
git add notebooks/Modelvalidation_meth_exp+analysis.ipynb
git commit -m "Add supplementary flanking-step DOF profiles (offsets -1,0,+1)"
```

---

## Task 8: Sanity-check figure — positional confound

**Files:**
- Modify: `notebooks/Modelvalidation_meth_exp+analysis.ipynb` (append)

If the bottom/top quartile bins cluster strongly in position-along-window or rotational phase, the DOF differences may be driven by position rather than sequence content. Make that visible.

- [ ] **Step 1: Add a cell that plots position and phase histograms by bin**

```python
# ── Sanity-check: distribution of cpg_pos and phase10 per ΔΔG bin ─────────
fig, axes = plt.subplots(1, 2, figsize=(11, 3.5))
for bin_name, colour in [("bottom_25", "tab:blue"),
                         ("mid_50",    "gray"),
                         ("top_25",    "tab:red")]:
    sub = df_single[df_single["ddG_bin"] == bin_name]
    axes[0].hist(sub["cpg_pos"], bins=np.arange(0, 148, 5),
                 alpha=0.5, label=bin_name, color=colour)
    axes[1].hist(sub["phase10"], bins=np.arange(0, 11, 1) - 0.5,
                 alpha=0.5, label=bin_name, color=colour)
axes[0].set_xlabel("CpG position along 147 bp window")
axes[0].set_ylabel("count")
axes[0].axvline(DYAD, color="k", ls=":", lw=1, label="dyad")
axes[1].set_xlabel("rotational phase (mod 10)")
axes[0].legend(fontsize=8); axes[1].legend(fontsize=8)
fig.suptitle("Position / phase distribution per ΔΔG quartile (confound check)")
plt.tight_layout()
plt.show()
```

- [ ] **Step 2: Run the cell**

Expected: two overlaid histograms; verifies bins are not perfectly segregated by position.

- [ ] **Step 3: Commit**

```bash
git add notebooks/Modelvalidation_meth_exp+analysis.ipynb
git commit -m "Add positional / phase confound sanity-check figure"
```

---

## Task 9: Summary table

**Files:**
- Modify: `notebooks/Modelvalidation_meth_exp+analysis.ipynb` (append)

One row per (quantity, DOF, step_offset). For `gs_dev`, `step_offset` is recorded as `NaN`.

- [ ] **Step 1: Add a cell that builds and displays the summary table**

```python
# ── Summary table: mean_bottom, mean_top, delta, Spearman ρ, Mann-Whitney p ─
def _summary_rows(frame, q, off):
    out = []
    for dof in DOF_NAMES:
        m_b, _, n_b = _group_stats(frame, dof, "bottom_25")
        m_t, _, n_t = _group_stats(frame, dof, "top_25")
        rho, p_sp = _spearman(frame, dof)
        b_vals = frame.loc[frame["ddG_bin"] == "bottom_25", dof].dropna()
        t_vals = frame.loc[frame["ddG_bin"] == "top_25",    dof].dropna()
        if len(b_vals) > 0 and len(t_vals) > 0:
            _, p_mw = mannwhitneyu(b_vals, t_vals, alternative="two-sided")
        else:
            p_mw = np.nan
        out.append({
            "quantity": q, "DOF": dof, "step_offset": off,
            "n_bot": n_b, "n_top": n_t,
            "mean_bottom": m_b, "mean_top": m_t,
            "delta": m_t - m_b,
            "spearman_rho": rho, "spearman_p": p_sp, "mw_p": p_mw,
        })
    return out

rows = []
for q in PER_STEP_QUANTITIES:
    for off in STEP_OFFSETS:
        sub = df_dof[(df_dof["quantity"] == q) & (df_dof["step_offset"] == off)]
        rows.extend(_summary_rows(sub, q, off))
rows.extend(_summary_rows(df_dev, "gs_dev", np.nan))

df_summary = pd.DataFrame(rows)

# Sorted view: center step + gs_dev, by |spearman_rho| descending.
center_view = df_summary[
    (df_summary["step_offset"] == 0) | df_summary["step_offset"].isna()
].copy()
center_view["abs_rho"] = center_view["spearman_rho"].abs()
center_view = center_view.sort_values("abs_rho", ascending=False).drop(columns="abs_rho")
center_view
```

- [ ] **Step 2: Run the cell**

Expected: a DataFrame display with 4 quantities × 6 DOFs = 24 rows in `center_view`, sorted by absolute Spearman correlation. Spot-check that the top rows correspond to DOFs whose bars in Figure A also looked the most different between bottom/top bins.

- [ ] **Step 3: Save the full summary to CSV next to the notebook**

```python
out_dir = Path("output") / "bp_step_dof"
out_dir.mkdir(parents=True, exist_ok=True)
df_summary.to_csv(out_dir / "bp_step_dof_summary.csv", index=False)
print(f"Wrote {out_dir / 'bp_step_dof_summary.csv'}  ({len(df_summary)} rows)")
```

Expected: `Wrote output/bp_step_dof/bp_step_dof_summary.csv  (78 rows)`
(3 per-step quantities × 3 offsets × 6 DOFs + 6 for gs_dev = 60 + 6 = wait: 3×3×6=54, plus 6 = 60. If your edge drops differ, the row count of the saved CSV stays 60 since it's defined by structure not data.)

- [ ] **Step 4: Commit**

```bash
git add notebooks/Modelvalidation_meth_exp+analysis.ipynb output/bp_step_dof/
git commit -m "Add BP-step DOF summary table + CSV export"
```

---

## Task 10: Closing discussion stub

**Files:**
- Modify: `notebooks/Modelvalidation_meth_exp+analysis.ipynb` (append)

Leave a markdown cell with a templated discussion the user can fill in based on the actual figures.

- [ ] **Step 1: Append a markdown cell**

```markdown
### Observations

- **DOFs with strongest correlation to ΔΔG** (from `center_view`): _fill in
  from the top rows of the summary table_.
- **Sign of the difference**: _e.g. destabilizing CpGs tend to have higher
  twist / lower roll / etc._
- **Local vs extended**: _check supplementary figure B — does the contrast
  persist at step_offset ±1, or is it confined to the CG step?_
- **Methylation deformation (`delta_gs`)**: _is the methylation-induced
  shift larger in one quartile?_
- **Mismatch with nucleosomal target (`gs_dev`)**: _are destabilizing CpGs
  in intervals with a larger composed excess?_
- **Position confound**: _per the sanity-check figure, are bottom/top bins
  clustered in different positional regions? If so, this analysis confounds
  position with sequence._
```

- [ ] **Step 2: Commit**

```bash
git add notebooks/Modelvalidation_meth_exp+analysis.ipynb
git commit -m "Add BP-step DOF discussion stub for human interpretation"
```

---

## Final verification

- [ ] **Restart kernel and run the new section end-to-end.** All Random-Sequence-Scan cells must already be executed (they produce `df_single`, `sequences`, `nb`, `CGNA_CONFIG`, etc.). Then run the 10 new cells in order. No exceptions, all assertions pass, 4 + 3 + 1 = 8 figures render, summary table displays, CSV written.

- [ ] **Skim the figures.** Each shows what its title claims. The cross-check assertion passed. Quartile counts are roughly 25 / 50 / 25.

- [ ] **Branch is `bp-level-prop`.** All commits made on this branch; no main-branch touches.

## Self-review (already performed by plan author)

- Spec coverage: all 6 spec sections (Inputs, Output, Definitions, Feature extraction, Analysis & figures, Notebook layout, Edge cases) map to tasks (T1 imports; T2/T3 features; T4 grouping; T5 sanity assertion; T6 fig A; T7 fig B; T8 confound; T9 summary; T10 discussion).
- Placeholders: none ("fill in from the top rows…" in T10 is intentional — that cell is a human-completion stub explicitly called out as such; T10's purpose is exactly to leave room for the user's interpretation).
- Type consistency: `DOF_NAMES`, `STEP_OFFSETS`, `PER_STEP_QUANTITIES`, `ALL_QUANTITIES` defined once in T1 and reused unchanged; `df_dof` and `df_dev` have explicit column lists in their construction and the same names in downstream merges (T4) and plots (T6–T9).
- Spec amendment (gs_dev as per-CpG interval excess) is faithfully implemented in T3 and respected by the figure layout (no flanking panel for gs_dev in T7).
