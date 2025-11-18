# src/config/var.py
# Created on 2025-04-06
# Author: Manish Yadav

"""
Custom data types for free energy calculations.

Note: For new landscape analysis code, see src/landscape/types.py which provides
improved versions with simplified subid (int instead of str) and optional sequence storage.

These types are maintained for backward compatibility with existing code.
"""

from typing import NamedTuple, Optional, List

class FreeEnergyResult(NamedTuple):
    """
    Container for energy calculation results.
    
    Used by CgnaNucFreeEnergy.calculate_free_energy_soft()
    """
    F: float
    F_entropy: float
    F_enthalpy: float
    F_freedna: float
    id: Optional[str]
    subid: Optional[str]


class SubsequenceTask(NamedTuple):
    """
    Represents a subsequence task for parallel processing.
    
    Attributes:
        id: Main sequence identifier
        subid: Subsequence identifier (e.g., 'subseq_0_147')
        sequence: The 147bp subsequence
        start_pos: Start position in the original sequence (0-based)
        end_pos: End position in the original sequence (exclusive)
        index: Subsequence index number
    """
    id: str
    subid: str
    sequence: str
    start_pos: int
    end_pos: int
    index: int


class LandscapeResult(NamedTuple):
    """
    Container for free energy landscape calculation results.
    
    Attributes:
        id: Main sequence identifier
        subid: Subsequence identifier
        sequence: The 147bp subsequence
        start_pos: Start position in original sequence
        end_pos: End position in original sequence
        index: Subsequence index
        dyad_position: Center position (dyad) in original sequence
        F: Total free energy
        F_entropy: Entropic contribution
        F_enthalpy: Enthalpic contribution
        F_freedna: Free DNA energy
        dF: Differential free energy (F - F_freedna)
        left: Left binding index used
        right: Right binding index used
    """
    id: str
    subid: str
    sequence: str
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
    Metadata for a main sequence being analyzed.
    
    Attributes:
        id: Sequence identifier
        length: Length of the full sequence
        num_subsequences: Number of 147bp subsequences generated
        window_size: Size of subsequence window (default 147)
        step_size: Step size between subsequences
        description: Optional description
    """
    id: str
    length: int
    num_subsequences: int
    window_size: int
    step_size: int
    description: Optional[str] = None

# class FreeEnergyResultSimple(NamedTuple):
#     """Container for simple energy calculation results."""
#     F: float
#     id: Optional[str]
#     subid: Optional[str]= None



# class ProcessedSequence(NamedTuple):
#     id: str
#     subid: str
#     sequence: str
#     start_site: int
#     end_site: int


# class NuclBreathingResult(NamedTuple):
#     """Container for nucleosome breathing results."""
#     id: str
#     subid: str
#     sequence: str
#     leftbind_indx: int
#     rightbind_indx: int
#     F_vals: FreeEnergyResult
#     Adsorp_F: float = 0.0


# class FREEDNAResult(NamedTuple):
#     """Container for free DNA energy results."""
#     id: str
#     subid: str
#     leftbind_indx: int
#     rightbind_indx: int
#     Ffree_bound: float
#     Ffree_unbound: float

# class SequenceTask(NamedTuple):
#     """
#     Represents a sequence task with id, subid, and sequence.
#     """
#     id: str
#     subid: str
#     sequence: str