"""
Refactored Nucleosome Free Energy Calculator.

This is a cleaner, modular version of CgnaNucFreeEnergy that separates
concerns into different sub-modules.
"""

import os
# Import environment settings if flag is set (default: import)
if os.environ.get("IMPORT_ENV_SETTINGS", "1") == "1":
    from src.config.env_settings import *


import numpy as np
from typing import Optional, Union
from src.config.custom_types import FreeEnergyResult
from femodules import NUC_STATE_PATH, K_POSRESC_PATH
from femodules.binding_sites import select_phosphate_sites, convert_to_open_sites
from femodules.energy_calc import calculate_multiharmonic_energy
from femodules.config import RBPConfig, CgnaConfig, DEFAULT_CGNA_CONFIG

from methods import (
    nucleosome_free_energy,
    read_nucleosome_triads,
    GenStiffness,
    calculate_midstep_triads
)
from binding_model import binding_model_free_energy
from methods.PolyCG.polycg import cgnaplus_bps_params


class NucleosomeBreathModular:
    """
    Modular calculator for nucleosome binding free energies.
    
    This is a refactored version of NucleosomeBreath with better organization:
    - Binding site logic → binding_sites.py
    - Energy calculations → energy_calc.py
    - Core orchestration → This class
    
    Parameters
    ----------
    config : Union[RBPConfig, CgnaConfig]
        Configuration object for either RBP or CGNA+ parameterization.
        Only one can be provided, not both.
        
        - RBPConfig: For Rigid Base Pair model
        - CgnaConfig: For CGNA+ model (default)
    
    Attributes
    ----------
    config : Union[RBPConfig, CgnaConfig]
        Active configuration
    param_type : str
        Type of parameterization ('rbp' or 'cgna')
    nuctriads : array
        Nucleosome triad conformations from crystal structure
    K_resc : array
        Rescaled stiffness matrix for binding sites
    nuc_mu0 : array
        Midstep triads for phosphate binding sites
    genstiff_nuc : GenStiffness
        Stiffness generator for nucleosome (only for RBP with multiharmonic)
    genstiff_freedna : GenStiffness
        Stiffness generator for free DNA (only for RBP with multiharmonic)
    
    Examples
    --------
    >>> # CGNA+ parameterization (default)
    >>> from femodules.config import CgnaConfig
    >>> config = CgnaConfig()
    >>> nb = NucleosomeBreathModular(config)
    >>> result = nb.calculate_free_energy("ACGT"*37, left=0, right=13)
    
    >>> # RBP parameterization
    >>> from femodules.config import RBPConfig
    >>> config = RBPConfig(nuc_method='crystal')
    >>> nb = NucleosomeBreathModular(config)
    
    >>> # RBP with multiharmonic
    >>> config = RBPConfig(nuc_method='crystal', free_dna_method='free_dna')
    >>> nb = NucleosomeBreathModular(config)
    
    >>> # Custom CGNA+ parameters
    >>> config = CgnaConfig(parameter_set_name='cgnaplus_bsc1')
    >>> nb = NucleosomeBreathModular(config)
    
    See Also
    --------
    femodules.config : Configuration classes (RBPConfig, CgnaConfig)
    femodules.binding_sites : Binding site configuration utilities
    femodules.energy_calc : Energy calculation functions
    """
    
    def __init__(
        self,
        config: Optional[Union[RBPConfig, CgnaConfig]] = None
    ):
        """
        Initialize the nucleosome free energy calculator.
        
        Parameters
        ----------
        config : Union[RBPConfig, CgnaConfig], optional
            Configuration object. If None, uses default CGNA+ configuration.
        
        Raises
        ------
        ValueError
            If config is not RBPConfig or CgnaConfig
        """
        # Use default CGNA config if none provided
        if config is None:
            config = DEFAULT_CGNA_CONFIG
        
        # Validate config type
        if not isinstance(config, (RBPConfig, CgnaConfig)):
            raise ValueError(
                f"config must be RBPConfig or CgnaConfig, got {type(config).__name__}"
            )
        
        self.config = config
        
        # Determine parameter type
        if isinstance(config, RBPConfig):
            self.param_type = 'rbp'
        else:
            self.param_type = 'cgna'
        
        # Load nucleosome structure data
        self._load_nucleosome_data()
        
        # Initialize GenStiffness for RBP parameterization
        if self.param_type == 'rbp':
            config_rbp = self.config  # type: RBPConfig
            self.genstiff_nuc = GenStiffness(method=config_rbp.nuc_method)
            # Initialize free DNA GenStiffness only if multiharmonic
            if config_rbp.is_multiharmonic():
                self.genstiff_freedna = GenStiffness(method=config_rbp.free_dna_method)
    
    def _load_nucleosome_data(self):
        """Load nucleosome structure and stiffness data."""
        # Load nucleosome triads from crystal structure
        self.nuctriads = read_nucleosome_triads(NUC_STATE_PATH)
        
        # Load rescaled stiffness matrix
        self.K_resc = np.load(K_POSRESC_PATH)
        
        # Calculate midstep triads for all binding sites
        phosphate_sites = select_phosphate_sites()
        self.nuc_mu0 = calculate_midstep_triads(
            triad_ids=phosphate_sites,
            nucleosome_triads=self.nuctriads
        )
    

    
    def calculate_free_energy(
        self,
        sequence: str,
        left: int,
        right: int,
        id: Optional[str] = None,
        subid: Optional[str] = None,
        kresc_factor: float = 1.0,
        style: str = "b_index",
        bound_ends: str = "exclude"
    ) -> FreeEnergyResult:
        """
        Calculate nucleosome binding free energy.
        
        This is the main method for calculating free energies using the
        soft binding model with CGNA+ parameterization.
        
        Parameters
        ----------
        sequence : str
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
            Adjusts binding strength (e.g., 0.5 = weaker, 2.0 = stronger)
        style : str, optional
            Binding specification style (default: 'b_index')
            Options: 'b_index', 'ph_index', 'open_sites'
        bound_ends : str, optional
            Treatment of outermost base pairs (default: 'exclude')
            Options: 'include', 'exclude'
        
        Returns
        -------
        FreeEnergyResult
            Named tuple containing:
            - F: Total free energy (kT)
            - F_entropy: Entropic contribution (kT)
            - F_enthalpy: Enthalpic contribution (kT)
            - F_freedna: Free DNA energy (kT)
            - id, subid: Identifiers
        
        Examples
        --------
        >>> # Standard calculation
        >>> result = nb.calculate_free_energy("ACGT"*37, 0, 13)
        >>> dF = result.F - result.F_freedna
        >>> print(f"Binding free energy: {dF:.2f} kT")
        
        >>> # Partial unwrapping
        >>> result = nb.calculate_free_energy(
        ...     "ACGT"*37,
        ...     left=2,
        ...     right=11,
        ...     style='b_index'
        ... )
        
        >>> # Weaker binding
        >>> result = nb.calculate_free_energy(
        ...     "ACGT"*37,
        ...     left=0,
        ...     right=13,
        ...     kresc_factor=0.5
        ... )
        
        Notes
        -----
        Uses either CGNA+ or RBP parameters with 6 DOF per base pair step.
        Energy units are in kT (thermal energy units).
        """
        # Generate stiffness parameters based on config type
        if self.param_type == 'cgna':
            config = self.config  # type: CgnaConfig
            gs, stiff = cgnaplus_bps_params(
                sequence=sequence,
                group_split=config.group_split,
                parameter_set_name=config.parameter_set_name,
            )
        else:  # rbp
            # Use pre-initialized GenStiffness from __init__
            stiff, gs = self.genstiff_nuc.gen_params(
                sequence,
                use_group=True,
                sparse=False
            )

        # Convert binding specification to open sites
        l_open, r_open = convert_to_open_sites(left, right, style)
        
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
        
        # Handle multiharmonic parameterization if specified (RBP only)
        if self.param_type == 'rbp' and self.config.is_multiharmonic():
            return self._calculate_with_multiharmonic(
                sequence, l_open, r_open, bound_ends, F_dict, id, subid
            )
        
        # Standard single harmonic result
        return FreeEnergyResult(
            F=F_dict['F'],
            F_entropy=F_dict['F_entropy'],
            F_enthalpy=F_dict['F_enthalpy'],
            F_freedna=F_dict['F_freedna'],
            id=id,
            subid=subid
        )
    
    def _calculate_with_multiharmonic(
        self,
        sequence: str,
        l_open: int,
        r_open: int,
        bound_ends: str,
        F_dict: dict,
        id: Optional[str],
        subid: Optional[str]
    ) -> FreeEnergyResult:
        """
        Calculate free energy with multiharmonic correction.
        
        Internal method that applies multiharmonic parameterization
        to separate bound and unbound DNA regions.
        """
        # Calculate separate energies for bound/unbound regions
        Fe_bound, Fe_unbound, Fe_total_free = calculate_multiharmonic_energy(
            seq=sequence,
            left_open=l_open,
            right_open=r_open,
            bound_ends=bound_ends,
            genstiff_bound=self.genstiff_nuc,
            genstiff_unbound=self.genstiff_freedna
        )
        
        # Recalculate total free energy with multiharmonic correction
        F_total = (F_dict['F'] - F_dict['F_freedna']) + Fe_bound + Fe_unbound
        F_entropy = F_total - F_dict['F_enthalpy']
        
        return FreeEnergyResult(
            F=F_total,
            F_entropy=F_entropy,
            F_enthalpy=F_dict['F_enthalpy'],
            F_freedna=Fe_total_free,
            id=id,
            subid=subid
        )
    
    def calculate_free_energy_hard(
        self,
        sequence: str,
        left: int,
        right: int,
        id: Optional[str] = None,
        subid: Optional[str] = None
    ) -> FreeEnergyResult:
        """
        Calculate free energy with hard binding model.
        
        This method uses a hard constraint model where specified
        phosphate positions are rigidly fixed to the nucleosome structure.
        
        Parameters
        ----------
        sequence : str
            DNA sequence (147bp)
        left : int
            Left bound-site index
        right : int
            Right bound-site index
        id : Optional[str], optional
            Sequence identifier
        subid : Optional[str], optional
            Subsequence identifier
        
        Returns
        -------
        FreeEnergyResult
            Free energy result with hard binding constraints
        
        Notes
        -----
        This method is less commonly used than calculate_free_energy().
        It enforces rigid binding at specified sites rather than soft constraints.
        Only supported for RBP parameterization.
        """
        if self.param_type != 'rbp':
            raise NotImplementedError(
                "Hard binding model only supported for RBP configuration"
            )
        
        # Generate stiffness parameters using pre-initialized GenStiffness
        stiff, gs = self.genstiff_nuc.gen_params(sequence, use_group=True)
        
        # Get binding sites
        mid_array = select_phosphate_sites(left, right)
        
        # Calculate free energy with hard constraints
        F_dict = nucleosome_free_energy(
            gs, stiff, mid_array, self.nuctriads, use_correction=True
        )
        
        return FreeEnergyResult(
            F=F_dict['F'],
            F_entropy=F_dict['F_entropy'],
            F_enthalpy=F_dict['F_enthalpy'],
            F_freedna=F_dict['F_freedna'],
            id=id,
            subid=subid
        )
