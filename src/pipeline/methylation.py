"""
Methylation data loading and MN encoding for CpG dinucleotides.
"""

from __future__ import annotations

import pandas as pd
from pathlib import Path


def find_cpg_positions(sequence: str) -> list[int]:
    """
    Return 0-based positions of all CG dinucleotides in *sequence*
    (position = index of C).
    """
    return [i for i in range(len(sequence) - 1) if sequence[i : i + 2] == "CG"]


def load_methylated_positions(
    cov_path: str | Path,
    sequence: str,
    chrom_num: str,
    seq_start: int,
    seq_end: int,
    threshold: float = 0.5,
) -> set[int]:
    """
    Load a Bismark .cov file and return the set of 0-based positions (in
    FASTA-index space) where CpGs are methylated above *threshold*.

    Parameters
    ----------
    cov_path : str | Path
        Path to the Bismark ``.cov`` file (6-column, tab-separated, no header).
    sequence : str
        Reference FASTA sequence (uppercase).  Used to confirm CG context.
    chrom_num : str
        Chromosome number **without** ``chr`` prefix (e.g. ``"17"``).
    seq_start : int
        1-based genomic start of the FASTA region.
    seq_end : int
        1-based genomic end of the FASTA region.
    threshold : float
        Minimum methylation fraction (0-1) to call a site methylated.

    Returns
    -------
    set[int]
        0-based FASTA positions of methylated C in confirmed CpG context.
    """
    cov_path = Path(cov_path)
    df = pd.read_csv(
        cov_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "meth_pct", "count_M", "count_U"],
        dtype={"chrom": str},
    )
    df["meth"] = df["meth_pct"] / 100.0

    region = df[
        (df["chrom"] == chrom_num)
        & (df["start"] >= seq_start)
        & (df["start"] < seq_end)
    ].copy()
    region["seq_pos"] = region["start"] - seq_start

    cpg_set = set(find_cpg_positions(sequence))
    methylated: set[int] = set()
    for _, row in region.iterrows():
        pos = int(row["seq_pos"])
        if pos in cpg_set and row["meth"] >= threshold:
            methylated.add(pos)
    return methylated


def apply_methylation_to_sequence(
    sequence: str,
    methylated_positions: set[int],
) -> str:
    """
    Replace methylated CpG dinucleotides with **MN** encoding.

    * ``C → M``  (methylated cytosine)
    * ``G → N``  (guanine paired with methylated C)

    Only positions where ``sequence[pos:pos+2] == "CG"`` **and** *pos* is in
    *methylated_positions* are modified.
    """
    seq = list(sequence)
    for pos in sorted(methylated_positions):
        if pos + 1 < len(seq) and seq[pos] == "C" and seq[pos + 1] == "G":
            seq[pos] = "M"
            seq[pos + 1] = "N"
    return "".join(seq)
