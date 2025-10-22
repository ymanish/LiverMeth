# src/config/var.py
# Created on 2025-04-06
# Author: Manish Yadav

from typing import NamedTuple, Optional

class FreeEnergyResult(NamedTuple):
    """Container for energy calculation results."""
    F: float
    F_entropy: float
    F_enthalpy: float
    F_freedna: float
    id: Optional[str]
    subid: Optional[str]

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