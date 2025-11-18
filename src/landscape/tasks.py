"""
Task generation for landscape analysis.

This module handles generating subsequences from DNA sequences.
"""

from typing import List, Tuple
from .types import SubsequenceTask, SequenceMetadata


def generate_subsequences(
    sequence: str,
    seq_id: str,
    window_size: int = 147,
    step_size: int = 1,
) -> Tuple[List[SubsequenceTask], SequenceMetadata]:
    """
    Generate subsequences from a DNA sequence with sliding window.
    
    Parameters
    ----------
    sequence : str
        DNA sequence to process (must be >= window_size)
    seq_id : str
        Identifier for the main sequence
    window_size : int, optional
        Size of each subsequence window in bp (default: 147)
    step_size : int, optional
        Step size between windows in bp (default: 1)
    
    Returns
    -------
    Tuple[List[SubsequenceTask], SequenceMetadata]
        List of subsequence tasks and sequence metadata
    
    Raises
    ------
    ValueError
        If sequence is shorter than window_size
    
    Examples
    --------
    >>> tasks, meta = generate_subsequences("ACGT"*50, "test_seq")
    >>> len(tasks)
    54
    >>> tasks[0].subid  # Start position
    0
    >>> tasks[1].subid
    1
    """
    seq_length = len(sequence)
    
    if seq_length < window_size:
        raise ValueError(
            f"Sequence length ({seq_length}) is shorter than "
            f"window size ({window_size})"
        )
    
    tasks: List[SubsequenceTask] = []
    index = 0
    
    for start_pos in range(0, seq_length - window_size + 1, step_size):
        end_pos = start_pos + window_size
        subseq = sequence[start_pos:end_pos]
        
        task = SubsequenceTask(
            id=seq_id,
            subid=start_pos,  # Simplified: just the start position
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
        window_size=window_size,
        step_size=step_size,
    )
    
    return tasks, metadata
