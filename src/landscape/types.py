"""
Data types for landscape analysis.

This module defines the core data structures used throughout the landscape
analysis pipeline.
"""

from typing import NamedTuple, Optional


class SubsequenceTask(NamedTuple):
    """
    Task for processing a single subsequence.
    
    Attributes
    ----------
    id : str
        Main sequence identifier
    subid : int
        Start position of subsequence (0-based index)
    sequence : str
        The DNA subsequence (typically 147bp)
    start_pos : int
        Start position in original sequence
    end_pos : int
        End position in original sequence (exclusive)
    index : int
        Sequential index of this subsequence
    """
    id: str
    subid: int
    sequence: str
    start_pos: int
    end_pos: int
    index: int


class LandscapeResult(NamedTuple):
    """
    Result of free energy calculation for a subsequence.
    
    Attributes
    ----------
    id : str
        Main sequence identifier
    subid : int
        Start position of subsequence
    sequence : Optional[str]
        The subsequence (None if store_sequences=False)
    start_pos : int
        Start position in original sequence
    end_pos : int
        End position in original sequence
    index : int
        Sequential index
    dyad_position : int
        Center position (dyad) in original sequence
    F : float
        Total free energy
    F_entropy : float
        Entropic contribution
    F_enthalpy : float
        Enthalpic contribution
    F_freedna : float
        Free DNA energy
    dF : float
        Differential free energy (F - F_freedna)
    left : int
        Left binding index used
    right : int
        Right binding index used
    """
    id: str
    subid: int
    sequence: Optional[str]
    start_pos: int
    end_pos: int
    index: int
    dyad_position: int
    F: float
    F_entropy: float
    F_enthalpy: float
    F_freedna: float
    dF: float
    left: int
    right: int


class SequenceMetadata(NamedTuple):
    """
    Metadata about a sequence being analyzed.
    
    Attributes
    ----------
    id : str
        Sequence identifier
    length : int
        Length of the full sequence
    num_subsequences : int
        Number of subsequences generated
    window_size : int
        Size of subsequence window
    step_size : int
        Step size between subsequences
    description : Optional[str]
        Optional description
    """
    id: str
    length: int
    num_subsequences: int
    window_size: int
    step_size: int
    description: Optional[str] = None
