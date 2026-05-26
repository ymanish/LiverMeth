"""
FASTA and BED file I/O, nucleosome-occupied region extraction.
"""

from __future__ import annotations

import pandas as pd
from Bio import SeqIO
from pathlib import Path
from typing import NamedTuple


# ── FASTA loading ────────────────────────────────────────────────────────────

class FastaRecord(NamedTuple):
    """Minimal container for a loaded FASTA record."""
    seq_id: str        # e.g. "chr17:35417357-35421983"
    sequence: str      # uppercase DNA string
    chrom: str         # e.g. "chr17"
    chrom_num: str     # e.g. "17" (without "chr")
    seq_start: int     # 1-based genomic start
    seq_end: int       # 1-based genomic end


def load_fasta(fasta_path: str | Path) -> FastaRecord:
    """
    Load the first record from a FASTA file and parse genomic coordinates
    from the header (expected format: ``>chr17:35417357-35421983``).

    Parameters
    ----------
    fasta_path : str | Path
        Path to the FASTA file.

    Returns
    -------
    FastaRecord
    """
    fasta_path = Path(fasta_path)
    record = next(SeqIO.parse(fasta_path, "fasta"))
    seq_id = record.id
    sequence = str(record.seq).upper()

    chrom, coords = seq_id.split(":")
    seq_start, seq_end = map(int, coords.split("-"))
    chrom_num = chrom.replace("chr", "")

    return FastaRecord(
        seq_id=seq_id,
        sequence=sequence,
        chrom=chrom,
        chrom_num=chrom_num,
        seq_start=seq_start,
        seq_end=seq_end,
    )


# ── ATAC peak loading ───────────────────────────────────────────────────────

def load_atac_peaks(
    bed_path: str | Path,
    chrom_num: str,
) -> pd.DataFrame:
    """
    Load an ATAC-seq BED file and filter to a single chromosome.

    BED format expected: chrom  start  end  name  (0-based, half-open).

    Parameters
    ----------
    bed_path : str | Path
        Path to the BED file.
    chrom_num : str
        Chromosome number without "chr" prefix (e.g. "17").

    Returns
    -------
    pd.DataFrame
        Sorted DataFrame with columns: chrom, start, end, name.
    """
    bed_path = Path(bed_path)
    df = pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name"],
        dtype={"chrom": str},
    )
    peaks = (
        df[df["chrom"] == chrom_num]
        .copy()
        .sort_values("start")
        .reset_index(drop=True)
    )
    return peaks


# ── Nucleosome-occupied regions ──────────────────────────────────────────────

def nucleosome_occupied_segments(
    peaks: pd.DataFrame,
    seq_start: int,
    seq_len: int,
    min_length: int = 147,
) -> list[tuple[int, int]]:
    """
    Compute the complement of ATAC peaks within the FASTA region — i.e. the
    nucleosome-occupied segments.

    Parameters
    ----------
    peaks : pd.DataFrame
        ATAC peaks DataFrame (BED 0-based coords) from :func:`load_atac_peaks`.
    seq_start : int
        1-based genomic start of the FASTA region.
    seq_len : int
        Length of the FASTA sequence.
    min_length : int
        Minimum segment length to keep (default: 147).

    Returns
    -------
    list[tuple[int, int]]
        List of (fasta_start, fasta_end) half-open intervals in FASTA-index
        coordinates (0-based).
    """
    region_0start = seq_start - 1  # genomic 0-based start

    # Convert BED → FASTA coords and clip
    fasta_starts = (peaks["start"] - region_0start).clip(lower=0)
    fasta_ends = (peaks["end"] - region_0start).clip(upper=seq_len)

    # Merge overlapping peaks
    merged: list[tuple[int, int]] = []
    for fs, fe in zip(fasta_starts, fasta_ends):
        if merged and fs <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], int(fe)))
        else:
            merged.append((int(fs), int(fe)))

    # Complement: gaps between / around peaks
    segments: list[tuple[int, int]] = []
    cursor = 0
    for ps, pe in merged:
        if cursor < ps:
            segments.append((cursor, ps))
        cursor = pe
    if cursor < seq_len:
        segments.append((cursor, seq_len))

    # Filter by minimum length
    segments = [(s, e) for s, e in segments if (e - s) >= min_length]
    return segments


def atac_peak_segments(
    peaks: pd.DataFrame,
    seq_start: int,
    seq_len: int,
    min_length: int = 147,
) -> list[tuple[int, int]]:
    """
    Return the ATAC peak regions themselves (as opposed to the complement).

    Parameters
    ----------
    peaks : pd.DataFrame
        ATAC peaks DataFrame (BED 0-based coords) from :func:`load_atac_peaks`.
    seq_start : int
        1-based genomic start of the FASTA region.
    seq_len : int
        Length of the FASTA sequence.
    min_length : int
        Minimum segment length to keep (default: 147).

    Returns
    -------
    list[tuple[int, int]]
        List of (fasta_start, fasta_end) half-open intervals in FASTA-index
        coordinates (0-based).
    """
    region_0start = seq_start - 1  # genomic 0-based start

    # Convert BED → FASTA coords and clip
    fasta_starts = (peaks["start"] - region_0start).clip(lower=0)
    fasta_ends = (peaks["end"] - region_0start).clip(upper=seq_len)

    # Merge overlapping peaks
    merged: list[tuple[int, int]] = []
    for fs, fe in zip(fasta_starts, fasta_ends):
        if merged and fs <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], int(fe)))
        else:
            merged.append((int(fs), int(fe)))

    # Filter by minimum length
    merged = [(s, e) for s, e in merged if (e - s) >= min_length]
    return merged


def write_atac_fasta(
    segments: list[tuple[int, int]],
    sequence: str,
    chrom: str,
    seq_start: int,
    out_path: str | Path,
    line_width: int = 60,
) -> Path:
    """
    Write a multi-record FASTA with one record per ATAC peak segment.

    Parameters
    ----------
    segments : list[tuple[int, int]]
        FASTA-coordinate segments from :func:`atac_peak_segments`.
    sequence : str
        Full FASTA sequence string.
    chrom : str
        Chromosome name (e.g. "chr17").
    seq_start : int
        1-based genomic start of the FASTA region.
    out_path : str | Path
        Destination path for the output FASTA file.
    line_width : int
        Characters per FASTA line (default: 60).

    Returns
    -------
    Path
        The output path (resolved).
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []
    for i, (s, e) in enumerate(segments):
        seg_seq = sequence[s:e]
        genomic_start = seq_start + s
        genomic_end = seq_start + e - 1
        header = f">{chrom}:{genomic_start}-{genomic_end}_atac{i + 1}"
        lines.append(header)
        lines.extend(seg_seq[j : j + line_width] for j in range(0, len(seg_seq), line_width))

    out_path.write_text("\n".join(lines) + "\n")
    return out_path


def write_nucleosome_fasta(
    segments: list[tuple[int, int]],
    sequence: str,
    chrom: str,
    seq_start: int,
    out_path: str | Path,
    line_width: int = 60,
) -> Path:
    """
    Write a multi-record FASTA with one record per nucleosome-occupied segment.

    Parameters
    ----------
    segments : list[tuple[int, int]]
        FASTA-coordinate segments from :func:`nucleosome_occupied_segments`.
    sequence : str
        Full FASTA sequence string.
    chrom : str
        Chromosome name (e.g. "chr17").
    seq_start : int
        1-based genomic start of the FASTA region.
    out_path : str | Path
        Destination path for the output FASTA file.
    line_width : int
        Characters per FASTA line (default: 60).

    Returns
    -------
    Path
        The output path (resolved).
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []
    for i, (s, e) in enumerate(segments):
        seg_seq = sequence[s:e]
        genomic_start = seq_start + s
        genomic_end = seq_start + e - 1
        header = f">{chrom}:{genomic_start}-{genomic_end}_nuc{i + 1}"
        lines.append(header)
        lines.extend(seg_seq[j : j + line_width] for j in range(0, len(seg_seq), line_width))

    out_path.write_text("\n".join(lines) + "\n")
    return out_path
