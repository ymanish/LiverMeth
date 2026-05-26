# Base-pair-step DOF analysis at stabilizing vs destabilizing CpGs

**Branch:** `bp-level-prop`
**Date:** 2026-05-26
**Status:** Design

## Motivation

The Random-Sequence Scan in `notebooks/Modelvalidation_meth_exp+analysis.ipynb`
shows that the single-site methylation effect ΔΔG varies strongly with CpG
position along the 147 bp window — some positions are systematically
stabilizing (ΔΔG < 0), others destabilizing (ΔΔG > 0). We want to know
**what, at the base-pair-step level, distinguishes the most-stabilizing from
the most-destabilizing CpG sites**.

Concretely: are the six rigid-body DOFs of the CG step (and its immediate
neighbours) — shift, slide, rise, tilt, roll, twist — systematically
different between the two groups? Does methylation deform those DOFs more in
one group than the other? And is the deviation from the nucleosomal target
larger at destabilizing positions?

## Inputs

Reused from the existing notebook section (no recomputation of the energy
scan):

- `df_single` — DataFrame of single-site methylation results
  (cols: `seq_label, cpg_pos, dG, dG_un, ddG, dist_from_dyad, phase10`, …).
  784 rows across 50 sequences.
- `sequences` — list of 50 dicts with `label`, `sequence`, `cpg_positions`.
- `CGNA_CONFIG` — the active CGNA+ parameter set.
- `cgnaplus_bps_params` — generates `(gs, stiff)` for a sequence;
  `gs[i]` is the 6-DOF groundstate for step `i`.
- `nb.nuc_mu0` — the nucleosomal midstep-triad target for the bound state;
  used to derive a per-step target in DOF space.

## Output

A new top-level section appended to
`notebooks/Modelvalidation_meth_exp+analysis.ipynb`:

> `## Base-pair-step DOF analysis at stabilizing vs destabilizing CpGs`

placed after the Random-Sequence Scan section. No new module file is created
unless the feature-extraction cell exceeds ~80 lines, in which case it is
lifted to `femodules/bp_step_features.py`.

## Definitions

**Step quantities** — for each CpG at position `p` and each step
`s ∈ {p-1, p, p+1}` we compute four 6-vectors:

| Name        | Definition                       | Interpretation                              |
|-------------|----------------------------------|---------------------------------------------|
| `gs_un`     | `cgnaplus(seq_un)[s]`            | Unmethylated equilibrium geometry           |
| `gs_meth`   | `cgnaplus(seq_meth)[s]`          | Methylated equilibrium geometry             |
| `delta_gs`  | `gs_meth − gs_un`                | Methylation-induced shift                   |
| `gs_dev`    | `gs_un − μ₀[s]`                  | Mismatch with nucleosomal target            |

`μ₀[s]` is derived from `nb.nuc_mu0` (the midstep triads at the nucleosomal
phosphate-binding positions) by converting consecutive triads into the
6 DOFs at each step `s`. Step indices in DNA coordinates align directly with
the nucleosomal step indices since the bound window is fully placed at
`left=0, right=13`.

**Step scope** — 3-step window centred on the CG step:
`step_offset ∈ {−1, 0, +1}`. CpGs where the window would fall outside
`[0, 145]` are dropped.

**Strand symmetrization** — not applied in v1. We report the forward CG step
only. (Easy to add later by averaging `gs(CG)` with reverse-complemented
`gs(GC)` in flipped coordinates.)

**Grouping** — ΔΔG quartile labels on `df_single`:
- `bottom_25` — ΔΔG ≤ 25th percentile (most stabilizing)
- `mid_50` — middle 50% (reference)
- `top_25` — ΔΔG ≥ 75th percentile (most destabilizing)

## Feature extraction

Build a long-form DataFrame `df_dof` with one row per
(CpG site × step_offset × quantity):

```
seq_label, cpg_pos, step_offset, quantity, shift, slide, rise, tilt, roll, twist, ddG, ddG_bin
```

Procedure:

1. For each unique `seq_label`, compute `gs_un` once and cache it (50 calls).
2. For each CpG site `p` in that sequence:
   - Build `seq_meth` by methylating only position `p`.
   - Compute `gs_meth`.
   - For each `step_offset ∈ {−1, 0, +1}` with valid index:
     - Emit four rows (one per quantity) with the 6-vector unpacked into
       columns.
3. Merge `ddG` and `ddG_bin` from `df_single` on `(seq_label, cpg_pos)`.

Approx. cost: 50 + 784 calls to `cgnaplus_bps_params`. Cheap; runs in the
foreground notebook cell, no parallelism needed.

## Analysis & figures

### Primary figures (approach A) — one per quantity, 4 total

For each quantity in `{gs_un, gs_meth, delta_gs, gs_dev}`:

- 1×6 panel grid (one per DOF).
- Each panel: bar of mean ± SEM for `bottom_25` vs `top_25` at
  `step_offset = 0`. Thin horizontal line for `mid_50` mean.
- Annotation in each panel: Spearman ρ of that DOF vs ΔΔG across all CpGs
  (continuous correlation; ignores binning), with p-value.

### Supplementary figures (approach B) — flanking-step profiles, 4 total

For each quantity, 2×3 grid (6 DOFs). Each panel: line plot of mean ± SEM
across `step_offset ∈ {−1, 0, +1}` for `bottom_25` vs `top_25`. Shows
whether the signal is local to CG or extends to flanks.

### Sanity-check panel — 1 figure

Two side-by-side histograms over `bottom_25` and `top_25`:
- CpG position along 147 bp window (`cpg_pos`).
- Rotational phase (`phase10`).

If the bins segregate strongly in position, the DOF differences may be
positional rather than sequence-driven; this confound is flagged in the
notebook discussion.

### Summary table — 1 cell

One row per `(quantity, DOF, step_offset)`:

```
quantity, DOF, step_offset, mean_bottom, mean_top, delta_mean, spearman_rho, mannwhitney_p
```

Rendered inline. Sorted by `|spearman_rho|` descending for the
`step_offset = 0` rows so the most-informative DOFs are at the top.

## Notebook layout

Approximately 7 new cells appended after the Random-Sequence Scan section:

1. Feature extraction → `df_dof`.
2. Quartile binning + merge with `df_single`.
3. Primary figures A (loop over 4 quantities).
4. Supplementary figures B (loop over 4 quantities).
5. Sanity-check histogram figure.
6. Summary table.
7. Short markdown cell summarising what the bottom-vs-top contrast reveals.

## Edge cases & decisions

- **Edge CpGs**: drop sites where `p-1 < 0` or `p+1 ≥ 146`.
- **Caching**: `gs_un` cached per `seq_label`; not per CpG.
- **μ₀ derivation**: convert `nb.nuc_mu0` triads to 6-DOF steps using the
  same triad-to-step convention as PolyCG. If the existing module exposes a
  helper for this, reuse it; otherwise implement a small inverse-triad
  routine in the notebook and cross-check by reconstructing one step.
- **Test / sanity assertion**: one assertion cell that recomputes ΔΔG for a
  single row from the cached `gs_meth` (running it through
  `nb.calculate_free_energy` with the methylated sequence) and confirms it
  matches `df_single.ddG` to within numerical tolerance. Catches any
  sequence-or-indexing mismatch between feature extraction and the original
  scan.

## Out of scope (v1)

- Strand symmetrization (forward CG only).
- Multivariate predictors of ΔΔG from DOFs (a regression model is a natural
  follow-up but not part of this analysis).
- Stiffness-matrix features — we only look at groundstate DOFs here.
- Applying the same analysis to non-random sequences (e.g., the TNF set).
