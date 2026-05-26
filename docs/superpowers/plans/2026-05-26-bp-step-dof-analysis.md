# BP-step DOF change on single-site methylation — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a notebook section that plots the mean ΔDOF (six panels — tilt, roll, twist, shift, slide, rise) at the methylated CG step as a function of position along the 147 bp window.

**Architecture:** Two cells appended to `notebooks/Modelvalidation_meth_exp+analysis.ipynb`: a feature-extraction cell that loops over `df_single` and builds `df_dof_change`, and a plotting cell.

**Tech Stack:** `cgnaplus_bps_params`, `_revert_terminal_mn`, `apply_methylation_to_sequence`, pandas, matplotlib.

---

## Conventions

- **DOF order from CGNA+** (with `euler_definition=True`, `translations_in_nm=True`): `[tilt, roll, twist, shift, slide, rise]`. Rotations in radians, translations in nm.
- **Step index for CpG at position p**: the CG step is `gs[p]`. (`gs` has shape `(L−1, 6)` for sequence length `L=147`.)
- **Sequence canonicalisation**: pass sequences through `_revert_terminal_mn` before `cgnaplus_bps_params`, matching the existing scan.
- **Branch**: `bp-level-prop`. All commits on this branch.

---

## Task 1: Append section heading + feature extraction

**Files:**
- Modify: `notebooks/Modelvalidation_meth_exp+analysis.ipynb` (append after the Random-Sequence Scan section)

- [ ] **Step 1: Append a markdown cell**

```markdown
## Methylation-induced change in BP-step DOFs vs position

For each single-site methylation scanned above, we compute the change in
the six CGNA+ base-pair-step DOFs (tilt, roll, twist, shift, slide, rise)
at the methylated CG step, then average across the random-sequence ensemble
as a function of position along the 147 bp window.
```

- [ ] **Step 2: Append a code cell that builds `df_dof_change`**

```python
# ── Δ(BP-step DOFs) at the methylated CG step, per CpG-methylation event ──
DOF_NAMES = ["tilt", "roll", "twist", "shift", "slide", "rise"]

def _gs(seq):
    seq_can = _revert_terminal_mn(seq)
    gs, _ = cgnaplus_bps_params(
        sequence=seq_can,
        group_split=CGNA_CONFIG.group_split,
        parameter_set_name=CGNA_CONFIG.parameter_set_name,
    )
    return gs

# Cache gs_un per sequence (50 calls).
gs_un_cache = {s["label"]: _gs(s["sequence"]) for s in sequences}
seq_by_label = {s["label"]: s["sequence"] for s in sequences}

rows = []
for r in df_single.itertuples(index=False):
    label, p = r.seq_label, int(r.cpg_pos)
    gs_un = gs_un_cache[label]
    if p >= gs_un.shape[0]:        # CpG at the very last bp has no step p
        continue
    seq_m = apply_methylation_to_sequence(seq_by_label[label], {p})
    gs_m = _gs(seq_m)
    delta = gs_m[p] - gs_un[p]
    rows.append({"seq_label": label, "cpg_pos": p, **dict(zip(DOF_NAMES, delta))})

df_dof_change = pd.DataFrame(rows)
print(f"df_dof_change: {len(df_dof_change)} rows  (CpGs at last bp skipped: "
      f"{len(df_single) - len(df_dof_change)})")
df_dof_change.head(5)
```

- [ ] **Step 3: Run the cell**

Expected: `~784` rows (minus at most a handful of CpGs that sit at the last bp of a sequence, where no step `p` exists). DataFrame head shows six numeric DOF columns.

- [ ] **Step 4: Commit**

```bash
git add notebooks/Modelvalidation_meth_exp+analysis.ipynb
git commit -m "Extract ΔDOF at methylated CG step into df_dof_change"
```

---

## Task 2: Plot the six-panel positional figure

**Files:**
- Modify: `notebooks/Modelvalidation_meth_exp+analysis.ipynb` (append after Task 1)

- [ ] **Step 1: Append a plotting cell**

```python
# ── Six-panel plot: mean ± SD of ΔDOF vs methylated CpG position ──────────
import matplotlib.pyplot as plt

fig, axes = plt.subplots(6, 1, figsize=(10, 12), sharex=True)
positions = np.arange(WINDOW_SIZE)

# Aggregate per position
grouped = df_dof_change.groupby("cpg_pos")
mean_by_pos = grouped[DOF_NAMES].mean()
std_by_pos  = grouped[DOF_NAMES].std(ddof=1)

for ax, dof in zip(axes, DOF_NAMES):
    # individual events
    ax.scatter(df_dof_change["cpg_pos"], df_dof_change[dof],
               s=8, alpha=0.25, color="steelblue", label="individual CpG sites")
    # mean line + ±1 SD band
    p = mean_by_pos.index.to_numpy()
    m = mean_by_pos[dof].to_numpy()
    sd = std_by_pos[dof].to_numpy()
    ax.plot(p, m, color="crimson", lw=1.6, label=f"mean Δ{dof}")
    ax.fill_between(p, m - sd, m + sd, color="crimson", alpha=0.15, label="±1 SD")
    ax.axhline(0, color="k", lw=0.5, ls="--")
    ax.axvline(DYAD, color="k", lw=0.5, ls=":", label="dyad" if dof == DOF_NAMES[0] else None)
    ax.set_ylabel(f"Δ{dof}")

axes[0].legend(loc="upper right", fontsize=8)
axes[-1].set_xlabel("methylated CpG position (bp along 147 bp window)")
fig.suptitle("Methylation-induced change in BP-step DOFs at the methylated CG step")
plt.tight_layout()
plt.show()
```

- [ ] **Step 2: Run the cell**

Expected: one figure with six stacked panels sharing the x-axis. Each panel shows the scatter of per-event ΔDOF values, a red mean line, a ±1 SD band, the zero line, and the dyad position.

- [ ] **Step 3: Commit**

```bash
git add notebooks/Modelvalidation_meth_exp+analysis.ipynb
git commit -m "Plot ΔDOF vs methylated CpG position (6-panel figure)"
```

---

## Final verification

- [ ] Restart kernel, run the Random-Sequence Scan section to populate `df_single`, `sequences`, `CGNA_CONFIG`, `WINDOW_SIZE`, `DYAD`, then run the two new cells. No errors; figure renders.
- [ ] Working tree is clean; both commits are on `bp-level-prop`.

## Self-review

- Spec coverage: feature-extraction (T1) + six-panel plot (T2) — that's the whole spec.
- Placeholders: none.
- Type consistency: `DOF_NAMES` order matches the CGNA+ convention used everywhere else; `gs_un_cache` and `seq_by_label` defined and used in the same cell.
- Out-of-scope items are explicitly omitted (no quartile bins, no significance tests, no `gs_dev`, no flanking steps).
