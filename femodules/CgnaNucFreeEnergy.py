"""
Nucleosome Free Energy Calculator using CGNA+ Parameterization

This module calculates nucleosome binding free energies using the CGNA+
parameterization with 6 degrees of freedom per base pair step (similar to RBP model).

The module supports:
- Soft binding model with variable contact points
- Multiharmonic parameterization for free DNA
- Different binding site specifications (bound index, phosphate index, open sites)

Key Classes
-----------
NucleosomeBreath : Main calculator class
    Handles nucleosome free energy calculations with flexible binding configurations

References
----------
CGNA+: A sequence-dependent coarse-grain model of double-stranded DNA
Di Stefano, M. et al. (2019)
"""

import os
import numpy as np
from typing import Optional, Tuple

# Import environment settings if flag is set (default: import)
if os.environ.get("IMPORT_ENV_SETTINGS", "1") == "1":
    from src.config.env_settings import *

from src.config.custom_types import FreeEnergyResult
from femodules import NUC_STATE_PATH, K_POSRESC_PATH

from methods import (
    nucleosome_free_energy,
    read_nucleosome_triads,
    GenStiffness,
    calculate_midstep_triads
)
from binding_model import binding_model_free_energy
from methods.PolyCG.polycg import cgnaplus_bps_params


class NucleosomeBreath:
    """
    Calculator for nucleosome binding free energies.
    
    This class handles the calculation of nucleosome free energies using CGNA+
    parameterization with a soft binding model that allows variable contact points.
    
    Parameters
    ----------
    nuc_method : str, optional
        Method for nucleosome parameterization (default: 'crystal')
        Options: 'crystal', 'hybrid', etc.
    free_dna_method : Optional[str], optional
        Method for free DNA parameterization (default: None)
        Used for multiharmonic parameterization. If None, uses single harmonic.
    
    Attributes
    ----------
    nuc_method : str
        Nucleosome calculation method
    free_dna_method : Optional[str]
        Free DNA calculation method
    nuctriads : array
        Nucleosome triad conformations from crystal structure
    K_resc : array
        Rescaled stiffness matrix for binding sites
    nuc_mu0 : array
        Midstep triads for phosphate binding sites
    
    Examples
    --------
    >>> nb = NucleosomeBreath(nuc_method='crystal')
    >>> result = nb.calculate_free_energy_soft(
    ...     seq601=sequence,
    ...     left=0,
    ...     right=13
    ... )
    >>> print(f"Free energy: {result.F:.2f}")
    """

    def __init__(
        self,  
        nuc_method: str = 'crystal',
        free_dna_method: Optional[str] = None
    ):
        """Initialize NucleosomeBreath calculator."""
        self.nuc_method = nuc_method
        self.free_dna_method = free_dna_method

        # Load nucleosome structure data
        self.triadfn = NUC_STATE_PATH
        self.nuctriads = read_nucleosome_triads(self.triadfn)

        # Load rescaled stiffness matrix
        self.fn = K_POSRESC_PATH
        self.K_resc = np.load(self.fn)

        # Calculate midstep triads for all binding sites
        self.nuc_mu0 = calculate_midstep_triads(
            triad_ids=self._select_phosphate_bind_sites(), 
            nucleosome_triads=self.nuctriads
        )



    def _select_phosphate_bind_sites(self, left: int = 0, right: int = 13) -> list:
        """
        Select phosphate binding sites on the nucleosome.
        
        The nucleosome has 28 phosphate binding sites (14 on each strand).
        This method selects a subset based on bound-site indices.
        
        Parameters
        ----------
        left : int, optional
            Left bound-site index (0-13, default: 0)
        right : int, optional
            Right bound-site index (0-13, default: 13)
        
        Returns
        -------
        list
            Indices of selected phosphate binding sites
        
        Notes
        -----
        Each bound-site index corresponds to 2 phosphate positions.
        Total 28 phosphate sites = 14 bound sites × 2 phosphates per site.
        """
        # Crystal structure phosphate binding site positions
        phosphate_bind_sites = [
            2, 6, 14, 17, 24, 29, 34, 38, 
            45, 49, 55, 59, 65, 69, 76, 
            80, 86, 90, 96, 100, 107, 111, 
            116, 121, 128, 131, 139, 143
        ]
        
        return phosphate_bind_sites[left*2:(right*2)+2]

    def _get_left_right_open(
        self,
        left: int,
        right: int,
        style: str = "b_index"
    ) -> Tuple[int, int]:
        """
        Convert binding indices to number of open phosphate sites.
        
        This method converts different representations of binding configurations
        into the number of open (unbound) phosphate sites on each end.
        
        Parameters
        ----------
        left : int
            Left specification (interpretation depends on style)
        right : int
            Right specification (interpretation depends on style)
        style : str, optional
            Specification style (default: 'b_index')
            Options:
            - 'b_index': Bound-site indices (0-13)
            - 'ph_index': Phosphate-site indices (0-27)  
            - 'open_sites': Direct count of open phosphates
        
        Returns
        -------
        Tuple[int, int]
            (l_open, r_open) - number of open phosphates on left and right
        
        Raises
        ------
        ValueError
            If style is not one of the valid options
        
        Examples
        --------
        >>> nb = NucleosomeBreath()
        >>> l_open, r_open = nb._get_left_right_open(1, 11, 'b_index')
        >>> print(l_open, r_open)  # 2, 4
        
        Notes
        -----
        Total phosphate sites: 28 (14 per strand)
        Bound-site i covers phosphates (2*i, 2*i+1)
        """
        if style == "b_index":
            # Bound-site indices (0-13)
            # Each bound site covers 2 phosphates
            # Example: left=1, right=11 → l_open=2, r_open=4
            l_open = 2 * left
            r_open = 28 - (2 * right) - 2

        elif style == "ph_index":
            # Phosphate-site indices (0-27)
            # Direct phosphate specification
            # Example: left=1, right=26 → l_open=1, r_open=1
            l_open = left
            r_open = 28 - right - 1
           
        elif style == "open_sites":
            # Direct specification of open phosphates
            # No conversion needed
            # Example: left=2, right=4 → l_open=2, r_open=4
            l_open = left
            r_open = right

        else:
            raise ValueError(
                f"Invalid style '{style}'. "
                "Use 'b_index', 'ph_index', or 'open_sites'."
            )
        
        return l_open, r_open

    def calculate_free_energy_soft(
        self,
        seq601: str,
        left: int,
        right: int, 
        id: Optional[str] = None,
        subid: Optional[str] = None,
        kresc_factor: float = 1.0,
        style: str = "b_index",
        bound_ends: str = "exclude"
    ) -> FreeEnergyResult:
        """
        Calculate nucleosome binding free energy with soft binding model.
        
        This is the main calculation method using a soft binding model where
        DNA can have variable contact points with the histone core.
        
        Parameters
        ----------
        seq601 : str
            DNA sequence (typically 147bp for nucleosome)
        left : int
            Left binding specification (interpretation depends on style)
        right : int
            Right binding specification (interpretation depends on style)
        id : Optional[str], optional
            Sequence identifier
        subid : Optional[str], optional
            Subsequence identifier
        kresc_factor : float, optional
            Rescaling factor for K matrix (default: 1.0)
            Use to adjust binding strength
        style : str, optional
            Binding specification style (default: 'b_index')
            See _get_left_right_open() for options
        bound_ends : str, optional
            Treatment of outermost base pair steps (default: 'exclude')
            Options:
            - 'include': Include as bound DNA
            - 'exclude': Treat as separate from histone core
        
        Returns
        -------
        FreeEnergyResult
            Named tuple containing:
            - F: Total free energy
            - F_entropy: Entropic contribution
            - F_enthalpy: Enthalpic contribution
            - F_freedna: Free DNA energy
            - id, subid: Identifiers
        
        Examples
        --------
        >>> nb = NucleosomeBreath()
        >>> result = nb.calculate_free_energy_soft(
        ...     seq601="ACGT"*37,
        ...     left=0,
        ...     right=13
        ... )
        >>> print(f"ΔF = {result.F - result.F_freedna:.2f}")
        
        Notes
        -----
        Uses CGNA+ parameters with 6 DOF per base pair step.
        The binding model accounts for DNA bending and histone-DNA contacts.
        """
        
        # Generate CGNA+ parameters for the sequence
        gs, stiff = cgnaplus_bps_params(
            sequence=seq601, 
            group_split=True,
            parameter_set_name='Di_hmethyl_methylated-hemi_combine',
        )

        # Convert binding indices to open phosphates
        l_open, r_open = self._get_left_right_open(left, right, style)

        # Calculate binding free energy
        F_dict = binding_model_free_energy(
            gs,
            stiff,    
            self.nuc_mu0,
            self.K_resc * kresc_factor,
            left_open=l_open,
            right_open=r_open,
            use_correction=True,
        )

        # Handle multiharmonic parameterization if specified
        if self.free_dna_method is not None:
            # Calculate separate energies for bound and unbound DNA regions
            bound_freedna_fe, unbound_freedna_fe, new_freedna_fe = \
                self.calculate_bound_nonbound_dna_energy(
                    seq=seq601, 
                    left_open=l_open, 
                    right_open=r_open,
                    bound_ends=bound_ends
                )
            
            # Recalculate total free energy with multiharmonic correction
            F_601_new = (F_dict['F'] - F_dict['F_freedna']) + bound_freedna_fe + unbound_freedna_fe
            F_free_new = new_freedna_fe
            F_entropy_new = F_601_new - F_dict['F_enthalpy'] 
            F_enthalpy = F_dict['F_enthalpy']
            
            return FreeEnergyResult(F_601_new, F_entropy_new, F_enthalpy, F_free_new, id, subid)

        # Standard single harmonic calculation
        F601 = F_dict['F']
        F_entropy = F_dict['F_entropy']
        F_enthalpy = F_dict['F_enthalpy']
        F_freedna = F_dict['F_freedna']

        return FreeEnergyResult(F601, F_entropy, F_enthalpy, F_freedna, id, subid)
    
    def calculate_bound_nonbound_dna_energy(
        self,
        seq: str,
        left_open: int,
        right_open: int,
        bound_ends: str = "include"
    ) -> Tuple[float, float, float]:
        """
        Calculate free energy for bound and unbound DNA regions separately.
        
        Used for multiharmonic parameterization where different regions of DNA
        have different elastic properties.
        
        Parameters
        ----------
        seq : str
            DNA sequence
        left_open : int
            Number of open phosphates on left side
        right_open : int
            Number of open phosphates on right side
        bound_ends : str, optional
            Treatment of outermost base pair steps (default: 'include')
            - 'include': Include in bound DNA calculation
            - 'exclude': Treat separately from histone core
        
        Returns
        -------
        Tuple[float, float, float]
            (Fe_bound_dna, Fe_unbound_dna, new_free_dna_fe)
            - Fe_bound_dna: Free energy of bound DNA region
            - Fe_unbound_dna: Free energy of unbound DNA region
            - new_free_dna_fe: Total free DNA energy
        
        Raises
        ------
        ValueError
            If bound_ends is not 'include' or 'exclude'
            If stiffness matrices have incompatible sizes
        
        Notes
        -----
        This method is only called when free_dna_method is not None.
        It uses different parameterizations for bound vs unbound DNA.
        """

        if bound_ends not in {"include", "exclude"}:
            raise ValueError("bound_ends must be either 'include' or 'exclude'.")

        full_stiff_dna_unbound, gs_dna = self.genstiff_freedna.gen_params(seq, use_group=True)
        full_stiff_dna_bound, gs_dna = self.genstiff_nuc.gen_params(seq, use_group=True)

        if full_stiff_dna_unbound.shape[0] != full_stiff_dna_bound.shape[0]:
            raise ValueError("Stiffness matrices for unbound and bound DNA must have the same size.")

        phosphate_bind_sites = self._select_phosphate_bind_sites()

        logdet_sign, logdet = np.linalg.slogdet(full_stiff_dna_unbound)
        new_free_dna_fe = -0.5*len(full_stiff_dna_unbound)*np.log(2*np.pi) + 0.5*logdet



        if bound_ends == "include" and left_open == 0 and right_open == 0:
            logdet_sign_b, logdet_b = np.linalg.slogdet(full_stiff_dna_bound)
            Fe_bound_dna = -0.5*len(full_stiff_dna_bound)*np.log(2*np.pi) + 0.5*logdet_b
            return Fe_bound_dna, 0.0, new_free_dna_fe
        
        bound_locs = phosphate_bind_sites[left_open:len(phosphate_bind_sites)-right_open]
        
        start = 6 * bound_locs[0]
        end = 6 * (bound_locs[-1] + 1)
        bound_block_stiff = full_stiff_dna_bound[start:end, start:end]
        
        unbound_indices = np.r_[0:start, end:full_stiff_dna_bound.shape[0]] 
        unbound_stiffness = full_stiff_dna_unbound[unbound_indices][:, unbound_indices]


        logdet_sign_b, logdet_b = np.linalg.slogdet(bound_block_stiff)
        Fe_bound_dna = -0.5*len(bound_block_stiff)*np.log(2*np.pi) + 0.5*logdet_b
        
        
        logdet_sign_ub, logdet_ub = np.linalg.slogdet(unbound_stiffness)
        Fe_unbound_dna = -0.5*len(unbound_stiffness)*np.log(2*np.pi) + 0.5*logdet_ub

        return Fe_bound_dna, Fe_unbound_dna, new_free_dna_fe

    def calculate_free_energy_hard(self, seq147:str, left:int, right:int, id:Optional[str]=None, subid:Optional[str]=None)-> FreeEnergyResult:
        stiff, gs = self.genstiff_nuc.gen_params(seq147, use_group=True)



        mid_array = self._select_phosphate_bind_sites(left, right)
        F_dict = nucleosome_free_energy(gs, stiff, mid_array, self.nuctriads, use_correction=True)


        F601 = F_dict['F']
        F_entrop = F_dict['F_entropy']
        F_entalap = F_dict['F_enthalpy']
        F_free = F_dict['F_freedna']



        # return F601, F_entrop, F_entalap, F_free, F_Diff
        return FreeEnergyResult(F601, F_entrop, F_entalap, F_free, id, subid)



if __name__ == "__main__":
    import time
    start = time.perf_counter()
    Seq601 = "CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGMTCTAGCACCGCTTAAACGCAMGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT"
    # Seq601 = "CG"*73+"C"
    # Seq601 = "TT"*73+"A"
    # Seq601 = "CG"+"MN"*71+"CGC"
    # Seq601 = "CG"+"CG"*71+"CGC"
    # Seq601 = "CG"+"HK"*71+"CGC"


    # Example usage
    nuc_breath = NucleosomeBreath()
    result = nuc_breath.calculate_free_energy_soft(seq601=Seq601, left=0, right=13, style="b_index")
    print(f"Free energy: {result.F}, Entropy: {result.F_entropy}, Enthalpy: {result.F_enthalpy}, Free DNA energy: {result.F_freedna}, dF: {result.F - result.F_freedna}")


    # free_dna_energy.get_free_dna_energy(fullseq=Seq601, left=0, right=13, style="b_index")
    # print(f"Free

    # for i in range(10):
    #     nuc_breath = NucleosomeBreath(nuc_method='hybrid', free_dna_method=None)
    #     result = nuc_breath.calculate_free_energy_soft(seq601=Seq601, left=0, right=13, style="b_index")
    #     # print(f"Left: {i}, Right: 13, Result: {result}")
    #     print(f"Free energy: {result.F}, Entropy: {result.F_entropy}, Enthalpy: {result.F_enthalpy}, Free DNA energy: {result.F_freedna}, dF: {result.F - result.F_freedna}")

    # # # for i in range(0, 14):
    #     result = nuc_breath.calculate_free_energy_soft(seq601=Seq601, left=0, right=i, style="b_index")
    #     print(f"Left: 0, Right: {i}, Result: {result}")
    

    # for i in range(0, 14):
    #     for j in range(0, 14):
    #         if j >= i:
    #             result = nuc_breath.calculate_free_energy_soft(seq601=Seq601, left=i, right=j, style="b_index")
                # print(f"L:{i}, R:{j} ------ Free energy: {result.F}, Entropy: {result.F_entropy}, Enthalpy: {result.F_enthalpy}, Free DNA energy: {result.F_freedna}, dF: {result.F - result.F_freedna}")
    # result = nuc_breath.calculate_free_energy_soft(seq601=Seq601, left=0, right=13, style="b_index")
    end = time.perf_counter()
    print(f"Time taken: {end - start:.2f} seconds")




    
    # nuc_breath.calculate_nonbound_dna_energy(seq=Seq601, left_open=0, right_open=0, bound_ends="exclude")

