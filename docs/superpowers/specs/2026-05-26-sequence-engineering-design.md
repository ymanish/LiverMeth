# Engineering DNA Sequences with Targeted Methylation Response

**Branch:** `bp-level-prop`
**Date:** 2026-05-26
**Status:** Design

## Goal

Build a library of 147 bp DNA sequences in three classes, evaluated against
full-CpG methylation:

- `stable_on_meth`   — ΔΔG = G(full-meth) − G(unmeth) ≪ 0 (methylation stabilizes binding)
- `unstable_on_meth` — ΔΔG ≫ 0 (methylation destabilizes binding)
- `inert_to_meth`    — ΔΔG ≈ 0 (no CpGs, no methylation effect)

The design is rule-based: a Widom-style CpG-free backbone gives nucleosome
positioning; CpG injection at empirically-stabilizing or -destabilizing
positions drives the methylation response.

## Inputs

- `df` from the existing Random-Sequence Scan section of
  `notebooks/Modelvalidation_meth_exp+analysis.ipynb` (50 sequences,
  single-site methylation, mean ΔΔG per CpG position).
- Existing helpers: `cgnaplus_bps_params`, `_revert_terminal_mn`,
  `apply_methylation_to_sequence`, `NucleosomeBreath.calculate_free_energy`
  from `femodules/nucleosome_breath_modular.py`.
- `WINDOW_SIZE = 147`, `DYAD = 73`, `CGNA_CONFIG`, `N_WORKERS` already
  available from that notebook.

## Literature rules used

From `docs/pre19.pdf` (Zuiddam & Schiessel, PRE 99, 012422, 2019):
- 147 bp wrapping energy is a sum of *local* (≈trinucleotide) contributions.
- Lowest-energy nucleosomal sequences are ~80% GC, with CC/GG/GC/CG steps
  concentrated where the major groove faces the histone (integer SHL
  positions).
- Highest-energy sequences are A/T-rich with A-tracts of 5–6 bp.

From the methylation literature
(NAR 50, 1864, 2022; Sci Reports 3, 2121, 2013):
- Unmethylated CpGs prefer minor-groove-out positions; methylated CpGs
  prefer minor-groove-in positions.
- Methylation increases roll and decreases twist at CG steps.
- The net effect of full methylation on nucleosome stability depends on
  *where* the CpGs sit relative to SHLs.

These rules motivate the choice of backbone (positioning-favoring,
CpG-free) and the use of an empirical positional ΔΔG map (which already
encodes the SHL-dependent methylation response specific to our CGNA+
parameterization).

## Architecture

Five-stage pipeline, implemented as a small module
`src/pipeline/sequence_design.py` and exercised from a new notebook
`notebooks/Sequence_engineering.ipynb`:

```
  1. position_map     from existing df → P_stab, P_destab
       ↓
  2. widom_backbone   periodic, CpG-free, 147 bp
       ↓
  3. inject_cpgs      place CGs at chosen positions (per class)
       ↓
  4. score_full_meth  full-meth ΔΔG via the existing CGNA+ free-energy call
       ↓
  5. plot             ΔΔG distributions per class + top-10 picks
```

## Components

```python
def position_map(
    df: pd.DataFrame,
    stab_quantile: float = 0.25,
    destab_quantile: float = 0.75,
) -> dict:
    """{'stab':       np.ndarray of positions in bottom quantile of mean ΔΔG,
        'destab':     np.ndarray of positions in top quantile of mean ΔΔG,
        'mean_ddg':   np.ndarray of length WINDOW_SIZE-1, per-position mean.}"""

def widom_backbone(
    length: int = 147,
    dyad: int = 73,
    motif: str = "GGCCCAATTT",
    seed: int | None = None,
) -> str:
    """Periodic CpG-free nucleosome-favoring backbone. `motif` is a
       10-bp tile (GC-rich half at major-groove-in, AT-rich half at
       minor-groove-in). Phase is set so the GC-rich half is centered
       at the dyad (major groove faces the histone at SHL 0) and the
       AT-rich half is centered at SHL ±0.5. `seed`, if given,
       performs in-half A↔T and G↔C swaps for backbone diversity
       (rotational frame preserved, no CGs introduced)."""

def inject_cpgs(backbone: str, positions: Iterable[int]) -> str:
    """Return backbone with seq[p:p+2] = 'CG' at each p (non-overlapping).
       Raise on overlapping/adjacent positions."""

def score_full_meth(
    seq: str,
    nb: NucleosomeBreath,
) -> tuple[float, float, float]:
    """Compute (G_unmeth, G_meth_full, ΔΔG). The methylation pattern is
       every CpG in seq."""

def design_library(
    class_name: Literal["stable", "unstable", "inert"],
    pos_map: dict,
    nb: NucleosomeBreath,
    n_candidates: int = 100,
    cpg_count_range: tuple[int, int] = (4, 12),
    n_backbones: int = 20,
    seed: int = 0,
) -> pd.DataFrame:
    """Generate n_candidates sequences for one class. Returns a
       DataFrame with columns
       (class, backbone_seed, cpg_positions, n_cpg, G_un, G_meth, ddG, seq).
       For 'inert', n_cpg is forced to 0 (no injection)."""
```

## Per-class procedure

**Setup (once):**
1. Build `pos_map` from `df`.
2. Pre-build a pool of `n_backbones=20` backbones (different seeds).
3. Construct one `NucleosomeBreath` instance for scoring.

**stable_on_meth:**
```python
for i in range(n_candidates):
    backbone  = random.choice(backbone_pool)
    n_cpg     = rng.integers(*cpg_count_range)
    positions = rng.choice(pos_map['stab'], size=n_cpg, replace=False)
    seq       = inject_cpgs(backbone, positions)
    G_un, G_m, dd = score_full_meth(seq, nb)
```

**unstable_on_meth:** same loop with `pos_map['destab']`.

**inert_to_meth:** same loop with `n_cpg = 0` (just score the bare
backbones).

Combine into one `df_lib` (~300 rows).

## Compute

Each candidate requires 2 free-energy calls (unmeth + full-meth).
300 candidates × 2 = 600 calls, parallelized via `ProcessPoolExecutor`
with `N_WORKERS` — same order as the existing 50-sequence × ~16-CpG scan.

## Output / plots

1. **Class-separation histogram:** three overlaid ΔΔG distributions
   (one per class), plus a fourth overlay computed from the original 50
   random sequences (full-methylated) as a "did we beat random?" baseline.
   Vertical line at ΔΔG = 0.

2. **Top-10 picks per class:** sorted `df_lib.head(10)` per class.
   `stable` ascending by ΔΔG, `unstable` descending, `inert` ascending by
   `|ΔΔG|`. Columns: `class, n_cpg, G_un, G_meth, ddG, seq`.

## Out of scope (easy to add later)

- Graph-based shortest-path optimization (the `pre19.pdf` method) to
  guarantee extremes rather than rely on random sampling within the rule
  constraints.
- Iterative refinement (mutate-and-rescore loops, GA, simulated
  annealing).
- Partial methylation patterns (some CpGs methylated, others not).
- Positional CpG-density vs mean_ddg overlay (sanity-check panel).
- Hamming diversity / clustering metrics within each class.
- Cross-class transitions (single-mutation paths between stable and
  unstable libraries).
