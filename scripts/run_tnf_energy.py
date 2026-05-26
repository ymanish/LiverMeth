#!/usr/bin/env python
"""
Tnf Region — Nucleosome Free Energy Landscape Pipeline
=======================================================

Reads FASTA, ATAC peak BED, and methylation .cov files, then computes the
sliding-window nucleosome free energy across selected portions of the Tnf
region.

The ``--region`` flag selects which genomic regions to analyse:
  - ``nucleosome``  : complement of ATAC peaks (nucleosome-occupied)
  - ``atac``        : the ATAC peak regions themselves

Results are saved into separate subdirectories under the output root:
  ``<out-dir>/nucleosome_occupied/thresh_X.XX/``
  ``<out-dir>/atac_peaks/thresh_X.XX/``

Uses ``multiprocessing.Pool`` to distribute batches of subsequences.

Usage
-----
    python scripts/run_tnf_energy.py                                # defaults (nucleosome)
    python scripts/run_tnf_energy.py --region atac                  # ATAC peaks
    python scripts/run_tnf_energy.py --region nucleosome --batch-size 200 -j 4
    python scripts/run_tnf_energy.py --help
"""

from __future__ import annotations
import os
if os.environ.get("IMPORT_ENV_SETTINGS", "1") == "1":
    from src.config.env_settings import *  # Triggers env_settings import
import argparse
import sys
import time
from pathlib import Path

import pandas as pd

# ── project root on sys.path so ``src.*`` and ``femodules.*`` resolve ───────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.pipeline import (
    load_fasta,
    load_atac_peaks,
    nucleosome_occupied_segments,
    atac_peak_segments,
    write_nucleosome_fasta,
    write_atac_fasta,
    find_cpg_positions,
    load_methylated_positions,
    apply_methylation_to_sequence,
    nuc_window_generator,
    batcher,
    run_pool_energy,
    results_to_dataframe,
    fill_from_reference,
    filter_unchanged_windows,
    copy_energies_from_reference,
    get_states,
)


# ── Defaults ─────────────────────────────────────────────────────────────────

DATA_DIR = PROJECT_ROOT / "data" / "260121_data_manish"
DEFAULT_FASTA = DATA_DIR / "Tnf_Region.fasta"
DEFAULT_ATAC = DATA_DIR / "Tnf_ATAC_peaks.bed"
DEFAULT_COV_D0 = DATA_DIR / "Tnf_met_d0.cov"
DEFAULT_COV_D38 = DATA_DIR / "Tnf_met_4d38.cov"
DEFAULT_OUT_DIR = PROJECT_ROOT / "output" / "tnf_energy"


# ── CLI ──────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Nucleosome free energy landscape for the Tnf region.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    io = p.add_argument_group("I/O paths")
    io.add_argument("--fasta", type=Path, default=DEFAULT_FASTA, help="FASTA file for the region.")
    io.add_argument("--atac-bed", type=Path, default=DEFAULT_ATAC, help="ATAC peak BED file.")
    io.add_argument("--cov-d0", type=Path, default=DEFAULT_COV_D0, help="Bismark .cov for d0.")
    io.add_argument("--cov-d38", type=Path, default=DEFAULT_COV_D38, help="Bismark .cov for d38.")
    io.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="Output directory.")
    io.add_argument("--out-prefix", type=str, default="tnf_nuc_energy", help="CSV filename prefix.")

    calc = p.add_argument_group("Calculation parameters")
    calc.add_argument("--window-size", type=int, default=147, help="Sliding window size (bp).")
    calc.add_argument("--step-size", type=int, default=1, help="Step size (bp).")
    calc.add_argument("--batch-size", type=int, default=100, help="Windows per batch for Pool.")
    calc.add_argument("--style", choices=["b_index", "ph_index", "open_sites"],
                      default="b_index", help="Breathing style.")
    calc.add_argument("--only-fullbound", action="store_true",
                      help="Only evaluate the fully-bound state.")
    calc.add_argument("--parameterization", choices=["cgna", "rbp"], default="cgna",
                      help="Free energy parameterization.")
    calc.add_argument("--cgna-param-set", type=str,
                      default="Di_hmethyl_methylated-hemi_combine",
                      help="CGNA+ parameter set name.")
    calc.add_argument("--meth-threshold", type=float, default=0.5,
                      help="Methylation fraction threshold (0-1).")

    reg = p.add_argument_group("Region selection")
    reg.add_argument("--region", choices=["nucleosome", "atac"],
                     default="nucleosome",
                     help="Which genomic regions to analyse: "
                          "'nucleosome' (complement of ATAC peaks) or "
                          "'atac' (the ATAC peak regions themselves).")

    par = p.add_argument_group("Parallelism")
    par.add_argument("-j", "--n-workers", type=int, default=None,
                     help="Pool workers (default: cpu_count).")

    return p.parse_args()


# ── Pipeline ─────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()
    t0 = time.time()

    # ── Incorporate region and threshold into output directory ────────────────
    region_tag = "nucleosome_occupied" if args.region == "nucleosome" else "atac_peaks"
    thresh_tag = f"thresh_{args.meth_threshold:.2f}"
    args.out_dir = args.out_dir / region_tag / thresh_tag
    print(f"Region               : {args.region}")
    print(f"Methylation threshold : {args.meth_threshold}  →  {args.out_dir}")

    # ── 1. Load FASTA ────────────────────────────────────────────────────────
    print(f"[1/7] Loading FASTA : {args.fasta}")
    fasta = load_fasta(args.fasta)
    print(f"       {fasta.seq_id}  ({len(fasta.sequence)} bp)")

    # ── 2. Load ATAC peaks → segments ────────────────────────────────────────
    print(f"[2/7] Loading ATAC peaks : {args.atac_bed}")
    peaks = load_atac_peaks(args.atac_bed, chrom_num=fasta.chrom_num)

    if args.region == "nucleosome":
        segments = nucleosome_occupied_segments(
            peaks,
            seq_start=fasta.seq_start,
            seq_len=len(fasta.sequence),
            min_length=args.window_size,
        )
        total_bp = sum(e - s for s, e in segments)
        print(f"       {len(segments)} nucleosome-occupied segments  ({total_bp} bp total)")
    else:
        segments = atac_peak_segments(
            peaks,
            seq_start=fasta.seq_start,
            seq_len=len(fasta.sequence),
            min_length=args.window_size,
        )
        total_bp = sum(e - s for s, e in segments)
        print(f"       {len(segments)} ATAC peak segments  ({total_bp} bp total)")

    # ── 3. Write region FASTA (for reference) ────────────────────────────────
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.region == "nucleosome":
        ref_fasta_path = args.out_dir / "Tnf_nucleosome_occupied.fasta"
        write_nucleosome_fasta(
            segments, fasta.sequence, fasta.chrom, fasta.seq_start, ref_fasta_path
        )
    else:
        ref_fasta_path = args.out_dir / "Tnf_atac_peaks.fasta"
        write_atac_fasta(
            segments, fasta.sequence, fasta.chrom, fasta.seq_start, ref_fasta_path
        )
    print(f"[3/7] Region FASTA     : {ref_fasta_path}")

    # ── 4. Load methylation data ─────────────────────────────────────────────
    print(f"[4/7] Loading methylation data ...")
    meth_d0 = load_methylated_positions(
        args.cov_d0, fasta.sequence, fasta.chrom_num,
        fasta.seq_start, fasta.seq_end, threshold=args.meth_threshold,
    )
    meth_d38 = load_methylated_positions(
        args.cov_d38, fasta.sequence, fasta.chrom_num,
        fasta.seq_start, fasta.seq_end, threshold=args.meth_threshold,
    )
    print(f"       d0  methylated CpGs : {len(meth_d0)}")
    print(f"       d38 methylated CpGs : {len(meth_d38)}")

    # ── 5. Build the three sequence variants ─────────────────────────────────
    seq_orig = fasta.sequence
    seq_d0 = apply_methylation_to_sequence(fasta.sequence, meth_d0)
    seq_d38 = apply_methylation_to_sequence(fasta.sequence, meth_d38)

    # breathing states
    style_states = get_states(args.style, only_fullbound=args.only_fullbound)
    print(f"[5b]  Breathing style : {args.style}  "
          f"({'fully-bound only' if args.only_fullbound else f'{len(style_states)} states'}"
          f"  |  {args.n_workers or 'auto'} workers)")

    # ── 6. Per-condition energy calculation ──────────────────────────────────

    def _gen_windows(sequence, label):
        return list(nuc_window_generator(
            sequence, segments,
            seq_id=f"{fasta.seq_id}_{label}",
            window_size=args.window_size,
            step_size=args.step_size,
        ))

    def _make_batches(windows):
        if not windows:
            return []
        return list(batcher(iter(windows), size=args.batch_size))

    pool_kw = dict(
        style=args.style,
        style_states=style_states,
        parameterization=args.parameterization,
        cgna_parameter_set=args.cgna_param_set,
        n_workers=args.n_workers,
    )

    # ── 6a. Unmethylated (reference) ─────────────────────────────────────────
    print(f"\n[6a/Un] Computing unmethylated reference ...")
    un_windows = _gen_windows(seq_orig, "Un")
    t_un = time.time()
    results_un = run_pool_energy(
        _make_batches(un_windows), **pool_kw, label="Un",
    )
    df_un = results_to_dataframe(results_un)
    out_un = args.out_dir / f"{args.out_prefix}_Un.csv"
    df_un.to_csv(out_un, index=False)
    print(f"       {len(un_windows)} windows  |  {len(df_un)} rows  |  "
          f"{time.time() - t_un:.1f}s  |  {out_un.name}")

    # ── 6b. d0 methylated ────────────────────────────────────────────────────
    print(f"\n[6b/d0] Computing d0 methylated ...")
    d0_windows = _gen_windows(seq_d0, "d0")
    t_d0 = time.time()
    results_d0 = run_pool_energy(
        _make_batches(d0_windows), **pool_kw, skip_if_no_meth=True, label="d0",
    )
    df_d0 = results_to_dataframe(results_d0)
    n_nan = int(df_d0["F"].isna().sum())
    df_d0 = fill_from_reference(df_d0, df_un)
    n_filled = n_nan - int(df_d0["F"].isna().sum())
    out_d0 = args.out_dir / f"{args.out_prefix}_d0.csv"
    df_d0.to_csv(out_d0, index=False)
    print(f"       {len(d0_windows)} windows  |  skipped {n_nan}, filled {n_filled} from Un  |  "
          f"{time.time() - t_d0:.1f}s  |  {out_d0.name}")

    # ── 6c. d38 methylated (chained back-fill: d0 → Un) ─────────────────────
    print(f"\n[6c/d38] Computing d38 methylated ...")
    d38_windows = _gen_windows(seq_d38, "d38")
    changed, unchanged = filter_unchanged_windows(d38_windows, d0_windows)
    n_same_d0 = len(unchanged)
    print(f"        {len(d38_windows)} windows  |  {n_same_d0} identical to d0  |  "
          f"{len(changed)} to compute")

    t_d38 = time.time()
    # Compute only changed windows (worker still skips no-M/N ones)
    results_d38 = run_pool_energy(
        _make_batches(changed), **pool_kw, skip_if_no_meth=True, label="d38",
    )
    df_computed = results_to_dataframe(results_d38)

    # Fill worker-skipped (no M/N) windows from Un
    n_nan_un = int(df_computed["F"].isna().sum()) if not df_computed.empty else 0
    if n_nan_un:
        df_computed = fill_from_reference(df_computed, df_un)
    n_filled_un = n_nan_un - (int(df_computed["F"].isna().sum()) if not df_computed.empty else 0)

    # Copy energies from d0 for unchanged windows
    df_from_d0 = copy_energies_from_reference(unchanged, df_d0, style_states)

    # Combine and sort
    df_d38 = pd.concat([df_computed, df_from_d0], ignore_index=True)
    df_d38 = df_d38.sort_values(["id", "subid"]).reset_index(drop=True)

    out_d38 = args.out_dir / f"{args.out_prefix}_d38.csv"
    df_d38.to_csv(out_d38, index=False)
    print(f"        filled {n_filled_un} from Un, {n_same_d0} from d0  |  "
          f"{time.time() - t_d38:.1f}s  |  {out_d38.name}")

    # ── 7. Summary ───────────────────────────────────────────────────────────
    dt_total = time.time() - t0
    print(f"\n[7/7] Pipeline complete in {dt_total:.1f}s")
    print(f"       Output dir : {args.out_dir}")
    for f in sorted(args.out_dir.glob(f"{args.out_prefix}_*.csv")):
        print(f"         {f.name}")


if __name__ == "__main__":
    main()
