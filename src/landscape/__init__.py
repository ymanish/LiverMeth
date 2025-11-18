"""
Nucleosome Free Energy Landscape Analysis

A modular package for calculating nucleosome free energy landscapes across
long DNA sequences using concurrent processing.

Main Components
---------------
- LandscapeAnalyzer: High-level interface for analysis
- Types: Data structures (SubsequenceTask, LandscapeResult, SequenceMetadata)
- Tasks: Subsequence generation
- Workers: Parallel energy calculations
- I/O: Save/load results in various formats

Quick Start
-----------
>>> from src.landscape import LandscapeAnalyzer
>>> analyzer = LandscapeAnalyzer(window_size=147, step_size=10)
>>> results, meta = analyzer.calculate_landscape(sequence, "my_seq")
>>> analyzer.save_results(results, meta, "output")

For convenience, use the quick function:
>>> from src.landscape import calculate_landscape
>>> results, meta = calculate_landscape(sequence, "my_seq")
"""

from .analyzer import LandscapeAnalyzer
from .types import SubsequenceTask, LandscapeResult, SequenceMetadata
from .tasks import generate_subsequences
from .workers import calculate_subsequence_energy
from .io_utils import (
    results_to_dataframe,
    save_results,
    load_results,
    save_multiple_results,
)


# Convenience function for quick analysis
def calculate_landscape(
    sequence: str,
    seq_id: str,
    window_size: int = 147,
    step_size: int = 1,
    left: int = 0,
    right: int = 13,
    store_sequences: bool = True,
    verbose: bool = True,
    **kwargs
):
    """
    Quick function to calculate free energy landscape.
    
    This is a convenience wrapper around LandscapeAnalyzer for simple use cases.
    
    Parameters
    ----------
    sequence : str
        DNA sequence to analyze
    seq_id : str
        Identifier for the sequence
    window_size : int, optional
        Subsequence window size (default: 147)
    step_size : int, optional
        Step between windows (default: 1)
    left : int, optional
        Left binding index (default: 0)
    right : int, optional
        Right binding index (default: 13)
    store_sequences : bool, optional
        Store sequences in results (default: True)
    verbose : bool, optional
        Show progress (default: True)
    **kwargs
        Additional arguments passed to LandscapeAnalyzer
    
    Returns
    -------
    Tuple[List[LandscapeResult], SequenceMetadata]
        Landscape results and metadata
    
    Examples
    --------
    >>> from src.landscape import calculate_landscape
    >>> results, meta = calculate_landscape("ACGT"*200, "test")
    >>> print(f"Analyzed {len(results)} positions")
    """
    analyzer = LandscapeAnalyzer(
        window_size=window_size,
        step_size=step_size,
        **kwargs
    )
    
    return analyzer.calculate_landscape(
        sequence=sequence,
        seq_id=seq_id,
        left=left,
        right=right,
        store_sequences=store_sequences,
        verbose=verbose,
    )


__all__ = [
    'LandscapeAnalyzer',
    'SubsequenceTask',
    'LandscapeResult',
    'SequenceMetadata',
    'generate_subsequences',
    'calculate_subsequence_energy',
    'results_to_dataframe',
    'save_results',
    'load_results',
    'save_multiple_results',
    'calculate_landscape',
]
