"""
I/O utilities for landscape analysis results.

This module provides functions to save and load landscape analysis results
in various formats (pickle, CSV, DataFrame).
"""

import pickle
from pathlib import Path
from typing import List, Optional, Dict, Any
import pandas as pd
from .types import LandscapeResult, SequenceMetadata


def results_to_dataframe(
    results: List[LandscapeResult],
    metadata: Optional[SequenceMetadata] = None,
) -> pd.DataFrame:
    """
    Convert landscape results to a pandas DataFrame.
    
    Parameters
    ----------
    results : List[LandscapeResult]
        List of landscape results
    metadata : Optional[SequenceMetadata]
        Optional sequence metadata to include as columns
    
    Returns
    -------
    pd.DataFrame
        DataFrame with all results and optional metadata
    
    Examples
    --------
    >>> df = results_to_dataframe(results, metadata)
    >>> df.columns
    Index(['id', 'subid', 'sequence', 'start_pos', 'end_pos', ...])
    """
    # Convert results to list of dicts
    data = [result._asdict() for result in results]
    df = pd.DataFrame(data)
    
    # Add metadata columns if provided
    if metadata is not None:
        df['seq_length'] = metadata.length
        df['window_size'] = metadata.window_size
        df['step_size'] = metadata.step_size
        if metadata.description:
            df['description'] = metadata.description
    
    return df


def save_results(
    results: List[LandscapeResult],
    metadata: SequenceMetadata,
    output_dir: Path,
    save_pickle: bool = True,
    save_csv: bool = True,
) -> Dict[str, Path]:
    """
    Save landscape results to disk.
    
    Parameters
    ----------
    results : List[LandscapeResult]
        List of landscape results to save
    metadata : SequenceMetadata
        Sequence metadata
    output_dir : Path
        Directory to save results
    save_pickle : bool, optional
        Save results as pickle (default: True)
    save_csv : bool, optional
        Save results as CSV (default: True)
    
    Returns
    -------
    Dict[str, Path]
        Dictionary mapping format to saved file path
    
    Examples
    --------
    >>> paths = save_results(results, metadata, Path("output"))
    >>> paths['pickle']
    PosixPath('output/my_seq_landscape.pkl')
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    saved_paths = {}
    
    # Save as pickle
    if save_pickle:
        pickle_path = output_dir / f"{metadata.id}_landscape.pkl"
        data = {
            'results': results,
            'metadata': metadata,
        }
        with open(pickle_path, 'wb') as f:
            pickle.dump(data, f)
        saved_paths['pickle'] = pickle_path
    
    # Save as CSV
    if save_csv:
        csv_path = output_dir / f"{metadata.id}_landscape.csv"
        df = results_to_dataframe(results, metadata)
        df.to_csv(csv_path, index=False)
        saved_paths['csv'] = csv_path
    
    return saved_paths


def load_results(pickle_path: Path) -> tuple:
    """
    Load landscape results from pickle file.
    
    Parameters
    ----------
    pickle_path : Path
        Path to pickle file
    
    Returns
    -------
    tuple
        (results, metadata) tuple
    
    Examples
    --------
    >>> results, metadata = load_results(Path("output/my_seq_landscape.pkl"))
    >>> len(results)
    100
    """
    with open(pickle_path, 'rb') as f:
        data = pickle.load(f)
    
    return data['results'], data['metadata']


def save_multiple_results(
    all_results: List[tuple],  # List of (results, metadata) tuples
    output_dir: Path,
    combined_csv: bool = True,
) -> Dict[str, Any]:
    """
    Save results from multiple sequences.
    
    Parameters
    ----------
    all_results : List[tuple]
        List of (results, metadata) tuples for multiple sequences
    output_dir : Path
        Directory to save results
    combined_csv : bool, optional
        Create a combined CSV with all sequences (default: True)
    
    Returns
    -------
    Dict[str, Any]
        Dictionary with saved paths and statistics
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    saved_info = {
        'individual_files': [],
        'num_sequences': len(all_results),
        'total_subsequences': 0,
    }
    
    # Save individual files
    for results, metadata in all_results:
        paths = save_results(results, metadata, output_dir)
        saved_info['individual_files'].append(paths)
        saved_info['total_subsequences'] += len(results)
    
    # Create combined CSV
    if combined_csv and len(all_results) > 0:
        all_dfs = []
        for results, metadata in all_results:
            df = results_to_dataframe(results, metadata)
            all_dfs.append(df)
        
        combined_df = pd.concat(all_dfs, ignore_index=True)
        combined_path = output_dir / "combined_landscape.csv"
        combined_df.to_csv(combined_path, index=False)
        saved_info['combined_csv'] = combined_path
    
    return saved_info
