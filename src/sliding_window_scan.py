"""
Sliding window scan for nucleosome free energy calculation.

This module provides functionality to scan long DNA sequences using a sliding window
approach, calculating nucleosome free energy at each position using the NucleosomeBreath
class.
"""

from typing import List, Optional, Dict
import pandas as pd
import numpy as np
from pathlib import Path

from femodules.CgnaNucFreeEnergy import NucleosomeBreath
from src.config.custom_types import FreeEnergyResult


WINDOW_SIZE = 147  # Standard nucleosome DNA length


class SlidingWindowScanner:
    """
    Scanner for calculating nucleosome free energies across long DNA sequences
    using a sliding window approach.
    """

    def __init__(
        self,
        window_size: int = WINDOW_SIZE,
        step_size: int = 1,
    ):
        """
        Initialize the sliding window scanner.

        Parameters
        ----------
        window_size : int
            Size of the sliding window in base pairs (default: 147)
        step_size : int
            Step size for sliding the window (default: 1 bp)
        """
        self.window_size = window_size
        self.step_size = step_size
        
        # Initialize NucleosomeBreath instance
        self.nuc_breath = NucleosomeBreath()

    def scan_sequence(
        self,
        sequence: str,
        seq_id: Optional[str] = None,
        left: int = 0,
        right: int = 13,
        kresc_factor: float = 1.0,
        style: str = "b_index",
        bound_ends: str = "exclude",
        verbose: bool = True,
    ) -> pd.DataFrame:
        """
        Scan a DNA sequence with a sliding window and calculate free energies.

        Parameters
        ----------
        sequence : str
            DNA sequence to scan (any length >= window_size)
        seq_id : Optional[str]
            Identifier for the sequence (for tracking)
        left : int
            Left bound index for nucleosome binding (default: 0)
        right : int
            Right bound index for nucleosome binding (default: 13)
        kresc_factor : float
            Rescaling factor for K matrix (default: 1.0)
        style : str
            Style for left/right specification (default: "b_index")
        bound_ends : str
            How to treat bound ends: "include" or "exclude" (default: "exclude")
        verbose : bool
            Print progress information (default: True)

        Returns
        -------
        pd.DataFrame
            DataFrame with columns:
            - seq_id: sequence identifier
            - window_index: starting position of window (0-based)
            - start: start position in sequence
            - end: end position in sequence
            - window_seq: the 147bp window sequence
            - dyad_position: position of dyad in original sequence (center of window)
            - F: total free energy
            - F_entropy: entropic contribution
            - F_enthalpy: enthalpic contribution
            - F_freedna: free DNA energy
            - dF: differential free energy (F - F_freedna)
        """
        seq_length = len(sequence)
        
        if seq_length < self.window_size:
            raise ValueError(
                f"Sequence length ({seq_length}) is shorter than "
                f"window size ({self.window_size})"
            )

        if seq_id is None:
            seq_id = "sequence"

        # Calculate number of windows
        num_windows = (seq_length - self.window_size) // self.step_size + 1
        
        if verbose:
            print(f"[SlidingWindowScanner] Scanning sequence: {seq_id}")
            print(f"  Sequence length: {seq_length} bp")
            print(f"  Window size: {self.window_size} bp")
            print(f"  Step size: {self.step_size} bp")
            print(f"  Number of windows: {num_windows}")

        records: List[Dict] = []

        for i, start in enumerate(range(0, seq_length - self.window_size + 1, self.step_size)):
            end = start + self.window_size
            window_seq = sequence[start:end]
            
            # Calculate dyad position (center of window)
            dyad_position = start + self.window_size // 2

            # Calculate free energy for this window
            result = self.nuc_breath.calculate_free_energy_soft(
                seq601=window_seq,
                left=left,
                right=right,
                id=seq_id,
                subid=f"win_{start}_{end}",
                kresc_factor=kresc_factor,
                style=style,
                bound_ends=bound_ends,
            )

            records.append({
                "seq_id": seq_id,
                "window_index": i,
                "start": start,
                "end": end,
                "dyad_position": dyad_position,
                "window_seq": window_seq,
                "F": result.F,
                "F_entropy": result.F_entropy,
                "F_enthalpy": result.F_enthalpy,
                "F_freedna": result.F_freedna,
                "dF": result.F - result.F_freedna,
            })

            if verbose and (i + 1) % 10 == 0:
                print(f"  Progress: {i + 1}/{num_windows} windows completed")

        if verbose:
            print(f"[SlidingWindowScanner] Scan complete: {len(records)} windows analyzed")

        df = pd.DataFrame(records)
        return df

    def scan_multiple_sequences(
        self,
        sequences: List[str],
        seq_ids: Optional[List[str]] = None,
        left: int = 0,
        right: int = 13,
        kresc_factor: float = 1.0,
        style: str = "b_index",
        bound_ends: str = "exclude",
        verbose: bool = True,
    ) -> pd.DataFrame:
        """
        Scan multiple DNA sequences with sliding windows.

        Parameters
        ----------
        sequences : List[str]
            List of DNA sequences to scan
        seq_ids : Optional[List[str]]
            List of sequence identifiers (if None, will be auto-generated)
        left, right, kresc_factor, style, bound_ends : 
            Same as scan_sequence
        verbose : bool
            Print progress information

        Returns
        -------
        pd.DataFrame
            Combined DataFrame with results from all sequences
        """
        if seq_ids is None:
            seq_ids = [f"seq_{i}" for i in range(len(sequences))]

        if len(sequences) != len(seq_ids):
            raise ValueError("Number of sequences must match number of seq_ids")

        all_results = []

        for seq, seq_id in zip(sequences, seq_ids):
            if verbose:
                print(f"\n{'='*60}")
            
            df = self.scan_sequence(
                sequence=seq,
                seq_id=seq_id,
                left=left,
                right=right,
                kresc_factor=kresc_factor,
                style=style,
                bound_ends=bound_ends,
                verbose=verbose,
            )
            all_results.append(df)

        combined_df = pd.concat(all_results, ignore_index=True)
        
        if verbose:
            print(f"\n{'='*60}")
            print(f"[SlidingWindowScanner] All sequences scanned")
            print(f"  Total sequences: {len(sequences)}")
            print(f"  Total windows: {len(combined_df)}")

        return combined_df

    def save_results(
        self,
        df: pd.DataFrame,
        output_path: Path,
        include_sequences: bool = False,
    ):
        """
        Save scan results to CSV file.

        Parameters
        ----------
        df : pd.DataFrame
            Results DataFrame from scan_sequence or scan_multiple_sequences
        output_path : Path
            Path to output CSV file
        include_sequences : bool
            Whether to include window sequences in output (default: False)
            Set to False to reduce file size for large scans
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if not include_sequences and "window_seq" in df.columns:
            df_to_save = df.drop(columns=["window_seq"])
        else:
            df_to_save = df

        df_to_save.to_csv(output_path, index=False)
        print(f"[SlidingWindowScanner] Results saved to: {output_path}")


def scan_sequence_windows(
    sequence: str,
    seq_id: Optional[str] = None,
    left: int = 0,
    right: int = 13,
    kresc_factor: float = 1.0,
    style: str = "b_index",
    bound_ends: str = "exclude",
    nuc_method: str = "crystal",
    free_dna_method: Optional[str] = None,
    window_size: int = WINDOW_SIZE,
    step_size: int = 1,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Convenience function to scan a single sequence with a sliding window.

    This is a wrapper around SlidingWindowScanner for simple use cases.

    Parameters
    ----------
    sequence : str
        DNA sequence to scan
    seq_id : Optional[str]
        Sequence identifier
    left, right, kresc_factor, style, bound_ends :
        Parameters for calculate_free_energy_soft
    nuc_method, free_dna_method :
        Parameters for NucleosomeBreath initialization
    window_size : int
        Window size in bp (default: 147)
    step_size : int
        Step size for sliding (default: 1)
    verbose : bool
        Print progress

    Returns
    -------
    pd.DataFrame
        Scan results
    """
    scanner = SlidingWindowScanner(
        nuc_method=nuc_method,
        free_dna_method=free_dna_method,
        window_size=window_size,
        step_size=step_size,
    )

    return scanner.scan_sequence(
        sequence=sequence,
        seq_id=seq_id,
        left=left,
        right=right,
        kresc_factor=kresc_factor,
        style=style,
        bound_ends=bound_ends,
        verbose=verbose,
    )


if __name__ == "__main__":
    # Example usage
    import time
    from src.config.path import OUT_DIR

    # Example sequence (longer than 147 bp)
    test_seq = "CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT" + "ATCGATCGATCGATCG" * 10

    print(f"Test sequence length: {len(test_seq)} bp")

    # Scan with sliding window
    start_time = time.perf_counter()
    
    df_results = scan_sequence_windows(
        sequence=test_seq,
        seq_id="test_sequence",
        step_size=10,  # Use larger step for faster example
        verbose=True,
    )

    end_time = time.perf_counter()

    print(f"\nTime taken: {end_time - start_time:.2f} seconds")
    print(f"\nResults preview:")
    print(df_results.head())

    # Save results
    output_file = OUT_DIR / "sliding_window_example.csv"
    scanner = SlidingWindowScanner()
    scanner.save_results(df_results, output_file, include_sequences=False)
