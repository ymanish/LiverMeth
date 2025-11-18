"""
Free Energy Landscape Analysis Module

This module provides functionality to calculate nucleosome free energy landscapes
across long DNA sequences using concurrent processing for efficient computation.

Features:
- Generate 147bp subsequences from sequences of any length
- Concurrent processing using multiprocessing
- Structured data types for clear identification and handling
- Export results to CSV and pickle formats
- Progress tracking and error handling
"""

import os
import time
from typing import List, Optional, Dict, Tuple
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
import pickle

import pandas as pd
import numpy as np
from tqdm import tqdm

from src.config.custom_types import (
    SubsequenceTask,
    LandscapeResult,
    SequenceMetadata,
    FreeEnergyResult,
)


# Worker function for multiprocessing - must be at module level
def _calculate_subsequence_energy(
    task: SubsequenceTask,
    left: int,
    right: int,
    kresc_factor: float,
    style: str,
    bound_ends: str,
    nuc_method: str,
    free_dna_method: Optional[str],
) -> LandscapeResult:
    """
    Worker function to calculate free energy for a single subsequence.
    
    This function is called in separate processes, so it must import
    NucleosomeBreath within the function.
    """
    from femodules.CgnaNucFreeEnergy import NucleosomeBreath
    
    # Initialize NucleosomeBreath (done per process)
    nb = NucleosomeBreath(
        nuc_method=nuc_method,
        free_dna_method=free_dna_method
    )
    
    # Calculate free energy
    result: FreeEnergyResult = nb.calculate_free_energy_soft(
        seq601=task.sequence,
        left=left,
        right=right,
        id=task.id,
        subid=task.subid,
        kresc_factor=kresc_factor,
        style=style,
        bound_ends=bound_ends,
    )
    
    # Calculate dyad position (center of subsequence)
    dyad_position = task.start_pos + len(task.sequence) // 2
    
    # Create LandscapeResult
    landscape_result = LandscapeResult(
        id=task.id,
        subid=task.subid,
        sequence=task.sequence,
        start_pos=task.start_pos,
        end_pos=task.end_pos,
        index=task.index,
        dyad_position=dyad_position,
        F=result.F,
        F_entropy=result.F_entropy,
        F_enthalpy=result.F_enthalpy,
        F_freedna=result.F_freedna,
        dF=result.F - result.F_freedna,
        left=left,
        right=right,
    )
    
    return landscape_result


class LandscapeAnalyzer:
    """
    Analyzer for calculating nucleosome free energy landscapes across DNA sequences.
    
    This class handles:
    - Subsequence generation from long sequences
    - Concurrent processing for efficient calculation
    - Result storage and export
    """
    
    def __init__(
        self,
        nuc_method: str = "crystal",
        free_dna_method: Optional[str] = None,
        window_size: int = 147,
        step_size: int = 1,
        n_workers: Optional[int] = None,
    ):
        """
        Initialize the LandscapeAnalyzer.
        
        Parameters
        ----------
        nuc_method : str
            Nucleosome method for energy calculation (default: 'crystal')
        free_dna_method : Optional[str]
            Free DNA method for multiharmonic parameterization (default: None)
        window_size : int
            Size of subsequence window in bp (default: 147)
        step_size : int
            Step size between subsequences in bp (default: 1)
        n_workers : Optional[int]
            Number of worker processes. If None, uses os.cpu_count()
        """
        self.nuc_method = nuc_method
        self.free_dna_method = free_dna_method
        self.window_size = window_size
        self.step_size = step_size
        self.n_workers = n_workers or os.cpu_count()
        
        print(f"[LandscapeAnalyzer] Initialized with {self.n_workers} workers")
    
    def generate_subsequences(
        self,
        sequence: str,
        seq_id: str,
    ) -> Tuple[List[SubsequenceTask], SequenceMetadata]:
        """
        Generate 147bp subsequences from a long DNA sequence.
        
        Parameters
        ----------
        sequence : str
            DNA sequence of any length (must be >= window_size)
        seq_id : str
            Identifier for the main sequence
        
        Returns
        -------
        Tuple[List[SubsequenceTask], SequenceMetadata]
            List of subsequence tasks and metadata about the sequence
        """
        seq_length = len(sequence)
        
        if seq_length < self.window_size:
            raise ValueError(
                f"Sequence length ({seq_length}) is shorter than "
                f"window size ({self.window_size})"
            )
        
        tasks: List[SubsequenceTask] = []
        index = 0
        
        for start_pos in range(0, seq_length - self.window_size + 1, self.step_size):
            end_pos = start_pos + self.window_size
            subseq = sequence[start_pos:end_pos]
            subid = f"{seq_id}_subseq_{index:04d}_pos_{start_pos}_{end_pos}"
            
            task = SubsequenceTask(
                id=seq_id,
                subid=subid,
                sequence=subseq,
                start_pos=start_pos,
                end_pos=end_pos,
                index=index,
            )
            tasks.append(task)
            index += 1
        
        metadata = SequenceMetadata(
            id=seq_id,
            length=seq_length,
            num_subsequences=len(tasks),
            window_size=self.window_size,
            step_size=self.step_size,
        )
        
        return tasks, metadata
    
    def calculate_landscape(
        self,
        sequence: str,
        seq_id: str,
        left: int = 0,
        right: int = 13,
        kresc_factor: float = 1.0,
        style: str = "b_index",
        bound_ends: str = "exclude",
        verbose: bool = True,
    ) -> Tuple[List[LandscapeResult], SequenceMetadata]:
        """
        Calculate free energy landscape for a DNA sequence using concurrent processing.
        
        Parameters
        ----------
        sequence : str
            DNA sequence to analyze
        seq_id : str
            Identifier for the sequence
        left : int
            Left binding index (default: 0)
        right : int
            Right binding index (default: 13)
        kresc_factor : float
            Rescaling factor for K matrix (default: 1.0)
        style : str
            Style for left/right specification (default: "b_index")
        bound_ends : str
            How to treat bound ends: "include" or "exclude" (default: "exclude")
        verbose : bool
            Show progress bar and information (default: True)
        
        Returns
        -------
        Tuple[List[LandscapeResult], SequenceMetadata]
            List of landscape results and sequence metadata
        """
        # Generate subsequences
        if verbose:
            print(f"\n[LandscapeAnalyzer] Processing sequence: {seq_id}")
            print(f"  Sequence length: {len(sequence)} bp")
        
        tasks, metadata = self.generate_subsequences(sequence, seq_id)
        
        if verbose:
            print(f"  Generated {len(tasks)} subsequences")
            print(f"  Window size: {self.window_size} bp")
            print(f"  Step size: {self.step_size} bp")
            print(f"  Using {self.n_workers} worker processes")
        
        # Prepare worker function with fixed parameters
        worker_func = partial(
            _calculate_subsequence_energy,
            left=left,
            right=right,
            kresc_factor=kresc_factor,
            style=style,
            bound_ends=bound_ends,
            nuc_method=self.nuc_method,
            free_dna_method=self.free_dna_method,
        )
        
        # Process tasks concurrently
        results: List[LandscapeResult] = []
        start_time = time.time()
        
        with ProcessPoolExecutor(max_workers=self.n_workers) as executor:
            # Submit all tasks
            futures = {executor.submit(worker_func, task): task for task in tasks}
            
            # Collect results with progress bar
            if verbose:
                pbar = tqdm(total=len(tasks), desc="  Calculating energies")
            
            for future in as_completed(futures):
                try:
                    result = future.result()
                    results.append(result)
                    if verbose:
                        pbar.update(1)
                except Exception as e:
                    task = futures[future]
                    print(f"\n[ERROR] Failed to process {task.subid}: {e}")
            
            if verbose:
                pbar.close()
        
        # Sort results by index to maintain order
        results.sort(key=lambda x: x.index)
        
        elapsed_time = time.time() - start_time
        
        if verbose:
            print(f"  Completed in {elapsed_time:.2f} seconds")
            print(f"  Average time per subsequence: {elapsed_time/len(tasks):.3f} seconds")
        
        return results, metadata
    
    def calculate_multiple_landscapes(
        self,
        sequences: List[str],
        seq_ids: List[str],
        left: int = 0,
        right: int = 13,
        kresc_factor: float = 1.0,
        style: str = "b_index",
        bound_ends: str = "exclude",
        verbose: bool = True,
    ) -> Tuple[List[LandscapeResult], List[SequenceMetadata]]:
        """
        Calculate free energy landscapes for multiple sequences.
        
        Parameters
        ----------
        sequences : List[str]
            List of DNA sequences
        seq_ids : List[str]
            List of sequence identifiers
        left, right, kresc_factor, style, bound_ends :
            Same as calculate_landscape
        verbose : bool
            Show progress information
        
        Returns
        -------
        Tuple[List[LandscapeResult], List[SequenceMetadata]]
            Combined list of all results and list of metadata for each sequence
        """
        if len(sequences) != len(seq_ids):
            raise ValueError("Number of sequences must match number of seq_ids")
        
        all_results: List[LandscapeResult] = []
        all_metadata: List[SequenceMetadata] = []
        
        for seq, seq_id in zip(sequences, seq_ids):
            results, metadata = self.calculate_landscape(
                sequence=seq,
                seq_id=seq_id,
                left=left,
                right=right,
                kresc_factor=kresc_factor,
                style=style,
                bound_ends=bound_ends,
                verbose=verbose,
            )
            all_results.extend(results)
            all_metadata.append(metadata)
        
        if verbose:
            print(f"\n[LandscapeAnalyzer] All sequences processed")
            print(f"  Total sequences: {len(sequences)}")
            print(f"  Total subsequences: {len(all_results)}")
        
        return all_results, all_metadata
    
    def results_to_dataframe(
        self,
        results: List[LandscapeResult],
        include_sequences: bool = True,
    ) -> pd.DataFrame:
        """
        Convert landscape results to a pandas DataFrame.
        
        Parameters
        ----------
        results : List[LandscapeResult]
            List of landscape results
        include_sequences : bool
            Whether to include subsequence sequences (default: True)
        
        Returns
        -------
        pd.DataFrame
            DataFrame with all results
        """
        data = []
        for r in results:
            row = {
                'id': r.id,
                'subid': r.subid,
                'index': r.index,
                'start_pos': r.start_pos,
                'end_pos': r.end_pos,
                'dyad_position': r.dyad_position,
                'F': r.F,
                'F_entropy': r.F_entropy,
                'F_enthalpy': r.F_enthalpy,
                'F_freedna': r.F_freedna,
                'dF': r.dF,
                'left': r.left,
                'right': r.right,
            }
            if include_sequences:
                row['sequence'] = r.sequence
            data.append(row)
        
        df = pd.DataFrame(data)
        return df
    
    def save_results(
        self,
        results: List[LandscapeResult],
        metadata: List[SequenceMetadata],
        output_dir: Path,
        prefix: str = "landscape",
        save_pickle: bool = True,
        save_csv: bool = True,
        include_sequences_in_csv: bool = False,
    ):
        """
        Save landscape results to disk.
        
        Parameters
        ----------
        results : List[LandscapeResult]
            Landscape results to save
        metadata : List[SequenceMetadata]
            Sequence metadata
        output_dir : Path
            Output directory
        prefix : str
            Prefix for output files (default: "landscape")
        save_pickle : bool
            Save results as pickle file (default: True)
        save_csv : bool
            Save results as CSV file (default: True)
        include_sequences_in_csv : bool
            Include sequences in CSV (default: False, to reduce file size)
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save as pickle (preserves full data structure)
        if save_pickle:
            pickle_path = output_dir / f"{prefix}_results.pkl"
            with open(pickle_path, 'wb') as f:
                pickle.dump({'results': results, 'metadata': metadata}, f)
            print(f"[LandscapeAnalyzer] Saved pickle: {pickle_path}")
        
        # Save as CSV
        if save_csv:
            df = self.results_to_dataframe(results, include_sequences=include_sequences_in_csv)
            csv_path = output_dir / f"{prefix}_results.csv"
            df.to_csv(csv_path, index=False)
            print(f"[LandscapeAnalyzer] Saved CSV: {csv_path}")
            
            # Save metadata
            metadata_df = pd.DataFrame([m._asdict() for m in metadata])
            metadata_path = output_dir / f"{prefix}_metadata.csv"
            metadata_df.to_csv(metadata_path, index=False)
            print(f"[LandscapeAnalyzer] Saved metadata: {metadata_path}")
    
    @staticmethod
    def load_results(pickle_path: Path) -> Tuple[List[LandscapeResult], List[SequenceMetadata]]:
        """
        Load results from a pickle file.
        
        Parameters
        ----------
        pickle_path : Path
            Path to pickle file
        
        Returns
        -------
        Tuple[List[LandscapeResult], List[SequenceMetadata]]
            Loaded results and metadata
        """
        with open(pickle_path, 'rb') as f:
            data = pickle.load(f)
        return data['results'], data['metadata']


def calculate_landscape(
    sequence: str,
    seq_id: str,
    output_dir: Optional[Path] = None,
    left: int = 0,
    right: int = 13,
    kresc_factor: float = 1.0,
    style: str = "b_index",
    bound_ends: str = "exclude",
    nuc_method: str = "crystal",
    free_dna_method: Optional[str] = None,
    window_size: int = 147,
    step_size: int = 1,
    n_workers: Optional[int] = None,
    save_results: bool = True,
    verbose: bool = True,
) -> Tuple[List[LandscapeResult], SequenceMetadata, pd.DataFrame]:
    """
    Convenience function to calculate free energy landscape for a single sequence.
    
    Parameters
    ----------
    sequence : str
        DNA sequence to analyze
    seq_id : str
        Sequence identifier
    output_dir : Optional[Path]
        Directory to save results (if save_results=True)
    left, right, kresc_factor, style, bound_ends :
        Parameters for free energy calculation
    nuc_method, free_dna_method :
        Parameters for NucleosomeBreath
    window_size, step_size :
        Parameters for subsequence generation
    n_workers : Optional[int]
        Number of worker processes
    save_results : bool
        Whether to save results to disk
    verbose : bool
        Show progress information
    
    Returns
    -------
    Tuple[List[LandscapeResult], SequenceMetadata, pd.DataFrame]
        Results list, metadata, and DataFrame
    """
    analyzer = LandscapeAnalyzer(
        nuc_method=nuc_method,
        free_dna_method=free_dna_method,
        window_size=window_size,
        step_size=step_size,
        n_workers=n_workers,
    )
    
    results, metadata = analyzer.calculate_landscape(
        sequence=sequence,
        seq_id=seq_id,
        left=left,
        right=right,
        kresc_factor=kresc_factor,
        style=style,
        bound_ends=bound_ends,
        verbose=verbose,
    )
    
    df = analyzer.results_to_dataframe(results)
    
    if save_results and output_dir is not None:
        analyzer.save_results(
            results=results,
            metadata=[metadata],
            output_dir=output_dir,
            prefix=seq_id,
        )
    
    return results, metadata, df


if __name__ == "__main__":
    # Example usage
    from src.config.path import OUT_DIR
    
    # Test sequence
    test_seq = "CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT" + "ATCGATCG" * 30
    
    print(f"Test sequence length: {len(test_seq)} bp")
    
    # Calculate landscape
    results, metadata, df = calculate_landscape(
        sequence=test_seq,
        seq_id="test_sequence",
        output_dir=OUT_DIR / "landscape_test",
        step_size=10,  # Use larger step for faster example
        n_workers=4,
        verbose=True,
    )
    
    print("\nResults preview:")
    print(df.head())
    print(f"\nMetadata: {metadata}")
