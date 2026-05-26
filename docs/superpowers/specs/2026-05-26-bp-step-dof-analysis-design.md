# Base-pair-step DOF change on single-site methylation

**Branch:** `bp-level-prop`
**Date:** 2026-05-26
**Status:** Design (revised — minimal scope)

## Goal

For each single-site methylation already scanned in the Random-Sequence Scan
section, compute the change in each of the six CGNA+ base-pair-step DOFs at
the methylated CG step, and plot the ensemble-averaged ΔDOF as a function of
position along the 147 bp window — six panels, one per DOF.

Same plot style as the existing ΔΔG figure (individual points, mean line,
±1 SD band), but with six DOF panels stacked or laid out instead of one
energy panel.

## Inputs

- `sequences` and `df_single` from the existing Random-Sequence Scan
  section.
- `cgnaplus_bps_params` via the same call pattern as
  `nucleosome_breath_modular.py`:

  ```python
  gs, stiff = cgnaplus_bps_params(
      sequence=sequence,
      group_split=CGNA_CONFIG.group_split,
      parameter_set_name=CGNA_CONFIG.parameter_set_name,
  )
  ```

- Existing helpers: `_revert_terminal_mn`, `apply_methylation_to_sequence`.

## Output

A new section appended to
`notebooks/Modelvalidation_meth_exp+analysis.ipynb`:

> `## Methylation-induced change in BP-step DOFs vs position`

containing:
1. A feature-extraction cell that builds a small DataFrame `df_dof_change`.
2. One figure with six panels (one per DOF) showing mean ± SD of ΔDOF at
   the methylated CG step vs CpG position along the window.

## Procedure

For each row of `df_single` (one methylated CpG at position `p` in some
sequence `s`):

1. Compute `gs_un` for sequence `s` (cached per sequence; 50 calls total).
2. Methylate only position `p` and compute `gs_meth`.
3. Take `delta_gs = gs_meth[p] − gs_un[p]` — a 6-vector at the methylated
   CG step.
4. Emit one row: `seq_label, cpg_pos, tilt, roll, twist, shift, slide, rise`.

Aggregate across all CpG-methylation events: for each position `p` along
0..145, compute mean and standard deviation of each DOF over the events
that methylated position `p`.

Plot: six-panel figure, `x = p`, `y = ΔDOF(p)`. Individual events as
scatter, mean as line, ±1 SD as shaded band — matching the visual style of
the existing ΔΔG plot.

## Out of scope (easy to add later if useful)

- Flanking-step changes (Δgs at steps `p−1` and `p+1`).
- Stabilizing-vs-destabilizing quartile splits.
- Statistical significance tests.
- Deviation from a nucleosomal target.

These are deliberately omitted to keep the first pass simple.
