"""
Main analyzer class for landscape analysis.

This module provides the high-level LandscapeAnalyzer class that orchestrates
the entire analysis pipeline.
"""

import os
import time
from typing import List, Optional, Tuple
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial

from tqdm import tqdm

from .types import LandscapeResult, SequenceMetadata
from .tasks import generate_subsequences
from .workers import calculate_subsequence_energy
from .io_utils import (
    results_to_dataframe,
    save_results,
    load_results,
    save_multiple_results,
)


class LandscapeAnalyzer:
    """
    Analyzer for calculating nucleosome free energy landscapes.
    
    This class provides a high-level interface for analyzing DNA sequences,
    handling subsequence generation, parallel computation, and result storage.
    
    Parameters
    ----------
    nuc_method : str, optional
        Nucleosome method for energy calculation (default: 'crystal')
    free_dna_method : Optional[str], optional
        Free DNA method for multiharmonic parameterization (default: None)
    window_size : int, optional
        Size of subsequence window in bp (default: 147)
    step_size : int, optional
        Step size between subsequences in bp (default: 1)
    n_workers : Optional[int], optional
        Number of worker processes. If None, uses os.cpu_count()
    
    Examples
    --------
    >>> analyzer = LandscapeAnalyzer(window_size=147, step_size=10)
    >>> results, meta = analyzer.calculate_landscape(sequence, "my_seq")
    >>> analyzer.save_results(results, meta, Path("output"))
    """
    
    def __init__(
        self,
        nuc_method: str = "crystal",
        free_dna_method: Optional[str] = None,
        window_size: int = 147,
        step_size: int = 1,
        n_workers: Optional[int] = None,
    ):
        """Initialize the LandscapeAnalyzer."""
        self.nuc_method = nuc_method
        self.free_dna_method = free_dna_method
        self.window_size = window_size
        self.step_size = step_size
        self.n_workers = n_workers or os.cpu_count()
        
        print(f"[LandscapeAnalyzer] Initialized with {self.n_workers} workers")
    
    def calculate_landscape(
        self,
        sequence: str,
        seq_id: str,
        left: int = 0,
        right: int = 13,
        kresc_factor: float = 1.0,
        style: str = "b_index",
        bound_ends: str = "exclude",
        store_sequences: bool = True,
        verbose: bool = True,
    ) -> Tuple[List[LandscapeResult], SequenceMetadata]:
        """
        Calculate free energy landscape for a DNA sequence.
        
        Uses concurrent processing to efficiently analyze a long DNA sequence
        by breaking it into overlapping subsequences.
        
        Parameters
        ----------
        sequence : str
            DNA sequence to analyze
        seq_id : str
            Identifier for the sequence
        left : int, optional
            Left binding index (default: 0)
        right : int, optional
            Right binding index (default: 13)
        kresc_factor : float, optional
            Rescaling factor for K matrix (default: 1.0)
        style : str, optional
            Style for left/right specification (default: "b_index")
            Options: 'b_index', 'ph_index', 'open_sites'
        bound_ends : str, optional
            How to treat bound ends: "include" or "exclude" (default: "exclude")
        store_sequences : bool, optional
            Store sequences in results (default: True)
            Set to False to save memory for large datasets
        verbose : bool, optional
            Show progress bar and information (default: True)
        
        Returns
        -------
        Tuple[List[LandscapeResult], SequenceMetadata]
            List of landscape results and sequence metadata
        
        Examples
        --------
        >>> analyzer = LandscapeAnalyzer()
        >>> results, meta = analyzer.calculate_landscape(
        ...     sequence="ACGT"*200,
        ...     seq_id="test_seq",
        ...     store_sequences=False
        ... )
        """
        # Generate subsequences
        if verbose:
            print(f"\n[LandscapeAnalyzer] Processing sequence: {seq_id}")
            print(f"  Sequence length: {len(sequence)} bp")
        
        tasks, metadata = generate_subsequences(
            sequence, seq_id, self.window_size, self.step_size
        )
        
        if verbose:
            print(f"  Generated {len(tasks)} subsequences")
            print(f"  Window size: {self.window_size} bp")
            print(f"  Step size: {self.step_size} bp")
            print(f"  Storing sequences: {store_sequences}")
            print(f"  Using {self.n_workers} worker processes")
        
        # Prepare worker function with fixed parameters
        worker_func = partial(
            calculate_subsequence_energy,
            left=left,
            right=right,
            kresc_factor=kresc_factor,
            style=style,
            bound_ends=bound_ends,
            nuc_method=self.nuc_method,
            free_dna_method=self.free_dna_method,
            store_sequence=store_sequences,
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
                    print(f"\n[ERROR] Failed to process subid={task.subid}: {e}")
            
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
        store_sequences: bool = True,
        verbose: bool = True,
    ) -> List[Tuple[List[LandscapeResult], SequenceMetadata]]:
        """
        Calculate landscapes for multiple sequences.
        
        Parameters
        ----------
        sequences : List[str]
            List of DNA sequences
        seq_ids : List[str]
            List of identifiers for sequences
        left : int, optional
            Left binding index (default: 0)
        right : int, optional
            Right binding index (default: 13)
        kresc_factor : float, optional
            Rescaling factor (default: 1.0)
        style : str, optional
            Style specification (default: "b_index")
        bound_ends : str, optional
            Bound ends treatment (default: "exclude")
        store_sequences : bool, optional
            Store sequences in results (default: True)
        verbose : bool, optional
            Show progress (default: True)
        
        Returns
        -------
        List[Tuple[List[LandscapeResult], SequenceMetadata]]
            List of (results, metadata) tuples for each sequence
        """
        if len(sequences) != len(seq_ids):
            raise ValueError("Number of sequences must match number of IDs")
        
        all_results = []
        
        for sequence, seq_id in zip(sequences, seq_ids):
            results, metadata = self.calculate_landscape(
                sequence=sequence,
                seq_id=seq_id,
                left=left,
                right=right,
                kresc_factor=kresc_factor,
                style=style,
                bound_ends=bound_ends,
                store_sequences=store_sequences,
                verbose=verbose,
            )
            all_results.append((results, metadata))
        
        return all_results
    
    # Expose utility functions as methods
    def results_to_dataframe(self, results, metadata=None):
        """Convert results to DataFrame. See io_utils.results_to_dataframe."""
        return results_to_dataframe(results, metadata)
    
    def save_results(self, results, metadata, output_dir, **kwargs):
        """Save results to disk. See io_utils.save_results."""
        return save_results(results, metadata, output_dir, **kwargs)
    
    def load_results(self, pickle_path):
        """Load results from disk. See io_utils.load_results."""
        return load_results(pickle_path)
    
    def save_multiple_results(self, all_results, output_dir, **kwargs):
        """Save multiple results. See io_utils.save_multiple_results."""
        return save_multiple_results(all_results, output_dir, **kwargs)
