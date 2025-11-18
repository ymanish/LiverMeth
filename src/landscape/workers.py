"""
Worker functions for parallel energy calculations.

This module contains the worker function that runs in separate processes
to calculate free energies for subsequences.
"""

from typing import Optional
from .types import SubsequenceTask, LandscapeResult
from src.config.custom_types import FreeEnergyResult


def calculate_subsequence_energy(
    task: SubsequenceTask,
    left: int,
    right: int,
    kresc_factor: float,
    style: str,
    bound_ends: str,
    nuc_method: str,
    free_dna_method: Optional[str],
    store_sequence: bool = True,
) -> LandscapeResult:
    """
    Calculate free energy for a single subsequence.
    
    This function is executed in separate processes for parallel computation.
    It imports NucleosomeBreath inside the function to avoid pickling issues.
    
    Parameters
    ----------
    task : SubsequenceTask
        The subsequence task to process
    left : int
        Left binding index
    right : int
        Right binding index
    kresc_factor : float
        Rescaling factor for K matrix
    style : str
        Style for left/right specification ('b_index', 'ph_index', 'open_sites')
    bound_ends : str
        How to treat bound ends ('include' or 'exclude')
    nuc_method : str
        Nucleosome method for calculation
    free_dna_method : Optional[str]
        Free DNA method (for multiharmonic parameterization)
    store_sequence : bool, optional
        Whether to store the sequence in results (default: True)
    
    Returns
    -------
    LandscapeResult
        Free energy calculation results for this subsequence
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
        subid=str(task.subid),  # Convert int to str for compatibility
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
        sequence=task.sequence if store_sequence else None,
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
