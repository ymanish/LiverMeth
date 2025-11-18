"""
Free energy calculation utilities for multiharmonic parameterization.

This module handles the calculation of free energies for bound and unbound
DNA regions when using multiharmonic parameterization.
"""

import numpy as np
from typing import Tuple
from .binding_sites import select_phosphate_sites


def calculate_log_determinant_energy(stiffness_matrix: np.ndarray) -> float:
    """
    Calculate free energy from stiffness matrix log determinant.
    
    Uses the formula: F = -0.5 * N * log(2π) + 0.5 * log|K|
    where K is the stiffness matrix.
    
    Parameters
    ----------
    stiffness_matrix : np.ndarray
        Stiffness matrix (2D array)
    
    Returns
    -------
    float
        Free energy contribution
    
    Notes
    -----
    This calculation is based on the partition function for a harmonic
    approximation to the DNA elastic energy.
    """
    logdet_sign, logdet = np.linalg.slogdet(stiffness_matrix)
    n_dof = len(stiffness_matrix)
    return -0.5 * n_dof * np.log(2 * np.pi) + 0.5 * logdet


def split_bound_unbound_regions(
    full_stiffness_bound: np.ndarray,
    full_stiffness_unbound: np.ndarray,
    left_open: int,
    right_open: int,
    phosphate_sites: list
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Split stiffness matrices into bound and unbound regions.
    
    Parameters
    ----------
    full_stiffness_bound : np.ndarray
        Full stiffness matrix with bound parameterization
    full_stiffness_unbound : np.ndarray
        Full stiffness matrix with unbound parameterization
    left_open : int
        Number of open phosphates on left
    right_open : int
        Number of open phosphates on right
    phosphate_sites : list
        List of phosphate binding site positions
    
    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (bound_block_stiffness, unbound_stiffness)
        - bound_block_stiffness: Stiffness for bound region
        - unbound_stiffness: Stiffness for unbound regions
    
    Notes
    -----
    Each phosphate site corresponds to 6 DOF (degrees of freedom).
    The bound region is extracted as a contiguous block, while unbound
    regions (left and right) are concatenated.
    """
    # Get bound site locations
    bound_locs = phosphate_sites[left_open:len(phosphate_sites) - right_open]
    
    # Calculate DOF indices (6 DOF per site)
    start_idx = 6 * bound_locs[0]
    end_idx = 6 * (bound_locs[-1] + 1)
    
    # Extract bound region
    bound_block_stiffness = full_stiffness_bound[start_idx:end_idx, start_idx:end_idx]
    
    # Extract unbound regions (left + right)
    unbound_indices = np.r_[0:start_idx, end_idx:full_stiffness_bound.shape[0]]
    unbound_stiffness = full_stiffness_unbound[unbound_indices][:, unbound_indices]
    
    return bound_block_stiffness, unbound_stiffness


def calculate_multiharmonic_energy(
    seq: str,
    left_open: int,
    right_open: int,
    bound_ends: str,
    genstiff_bound,
    genstiff_unbound
) -> Tuple[float, float, float]:
    """
    Calculate free energies for bound/unbound DNA with multiharmonic parameters.
    
    This function handles the case where different elastic parameterizations
    are used for bound vs unbound DNA regions.
    
    Parameters
    ----------
    seq : str
        DNA sequence
    left_open : int
        Number of open phosphates on left side
    right_open : int
        Number of open phosphates on right side
    bound_ends : str
        Treatment of outermost base pair steps ('include' or 'exclude')
    genstiff_bound : GenStiffness
        Generator for bound DNA stiffness
    genstiff_unbound : GenStiffness
        Generator for unbound DNA stiffness
    
    Returns
    -------
    Tuple[float, float, float]
        (Fe_bound_dna, Fe_unbound_dna, Fe_total_freedna)
        - Fe_bound_dna: Free energy of bound DNA region
        - Fe_unbound_dna: Free energy of unbound DNA regions
        - Fe_total_freedna: Total free DNA energy (reference)
    
    Raises
    ------
    ValueError
        If bound_ends is not 'include' or 'exclude'
        If stiffness matrices have incompatible sizes
    
    Notes
    -----
    When bound_ends='include' and there are no open sites (fully bound),
    all DNA is treated as bound, and unbound energy is zero.
    """
    if bound_ends not in {"include", "exclude"}:
        raise ValueError("bound_ends must be either 'include' or 'exclude'.")
    
    # Generate stiffness matrices
    full_stiff_unbound, gs_dna = genstiff_unbound.gen_params(seq, use_group=True)
    full_stiff_bound, gs_dna = genstiff_bound.gen_params(seq, use_group=True)
    
    # Validate matrix sizes
    if full_stiff_unbound.shape[0] != full_stiff_bound.shape[0]:
        raise ValueError(
            "Stiffness matrices for unbound and bound DNA must have the same size."
        )
    
    # Calculate total free DNA energy (reference)
    Fe_total_freedna = calculate_log_determinant_energy(full_stiff_unbound)
    
    # Special case: fully bound with included ends
    if bound_ends == "include" and left_open == 0 and right_open == 0:
        Fe_bound_dna = calculate_log_determinant_energy(full_stiff_bound)
        return Fe_bound_dna, 0.0, Fe_total_freedna
    
    # Get phosphate binding sites
    phosphate_sites = select_phosphate_sites()
    
    # Split into bound and unbound regions
    bound_block_stiff, unbound_stiffness = split_bound_unbound_regions(
        full_stiff_bound,
        full_stiff_unbound,
        left_open,
        right_open,
        phosphate_sites
    )
    
    # Calculate energies for each region
    Fe_bound_dna = calculate_log_determinant_energy(bound_block_stiff)
    Fe_unbound_dna = calculate_log_determinant_energy(unbound_stiffness)
    
    return Fe_bound_dna, Fe_unbound_dna, Fe_total_freedna
