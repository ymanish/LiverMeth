"""
Configuration classes for nucleosome free energy calculations.

This module provides configuration classes for different parameterization methods:
- RBPConfig: Rigid Base Pair model configuration
- CgnaConfig: CGNA+ model configuration
"""

from typing import Optional, Literal
from dataclasses import dataclass


@dataclass
class RBPConfig:
    """
    Configuration for Rigid Base Pair (RBP) model calculations.
    
    Parameters
    ----------
    nuc_method : str
        Method for nucleosome parameterization
        Options: 'crystal', 'hybrid', etc.
    free_dna_method : Optional[str]
        Method for free DNA parameterization (default: None)
        Used for multiharmonic parameterization. If None, uses single harmonic.
    
    Examples
    --------
    >>> # Single harmonic
    >>> config = RBPConfig(nuc_method='crystal')
    
    >>> # Multiharmonic
    >>> config = RBPConfig(
    ...     nuc_method='crystal',
    ...     free_dna_method='free_dna'
    ... )
    """
    nuc_method: str = 'crystal'
    free_dna_method: Optional[str] = None
    
    def __post_init__(self):
        """Validate configuration after initialization."""
        valid_nuc_methods = ['crystal', 'hybrid', 'md', 'other']
        if self.nuc_method not in valid_nuc_methods:
            # Allow any string, but warn if not in common list
            pass
    
    def is_multiharmonic(self) -> bool:
        """Check if multiharmonic parameterization is enabled."""
        return self.free_dna_method is not None
    
    def __repr__(self):
        if self.is_multiharmonic():
            return (f"RBPConfig(nuc_method='{self.nuc_method}', "
                   f"free_dna_method='{self.free_dna_method}')")
        return f"RBPConfig(nuc_method='{self.nuc_method}')"


@dataclass
class CgnaConfig:
    """
    Configuration for CGNA+ model calculations.
    
    Parameters
    ----------
    parameter_set_name : str
        Name of the CGNA+ parameter set to use
        Default: 'Di_hmethyl_methylated-hemi_combine'
        
        Common options:
        - 'Di_hmethyl_methylated-hemi_combine': Methylation-aware parameters
        - 'cgnaplus_bsc1': Standard CGNA+ parameters
        - 'abc': ABC parameter set
    group_split : bool
        Whether to use group splitting in parameterization (default: True)
    
    Examples
    --------
    >>> # Default methylation-aware parameters
    >>> config = CgnaConfig()
    
    >>> # Custom parameter set
    >>> config = CgnaConfig(parameter_set_name='cgnaplus_bsc1')
    
    >>> # Without group splitting
    >>> config = CgnaConfig(group_split=False)
    """
    parameter_set_name: str = 'Di_hmethyl_methylated-hemi_combine'
    group_split: bool = True
    
    def __post_init__(self):
        """Validate configuration after initialization."""
        if not self.parameter_set_name:
            raise ValueError("parameter_set_name cannot be empty")
    
    def __repr__(self):
        return (f"CgnaConfig(parameter_set_name='{self.parameter_set_name}', "
               f"group_split={self.group_split})")


def create_rbp_config(
    nuc_method: str = 'crystal',
    free_dna_method: Optional[str] = None
) -> RBPConfig:
    """
    Convenience function to create RBP configuration.
    
    Parameters
    ----------
    nuc_method : str, optional
        Nucleosome method (default: 'crystal')
    free_dna_method : Optional[str], optional
        Free DNA method for multiharmonic (default: None)
    
    Returns
    -------
    RBPConfig
        Configuration object
    
    Examples
    --------
    >>> config = create_rbp_config('crystal')
    >>> config = create_rbp_config('crystal', 'free_dna')
    """
    return RBPConfig(nuc_method=nuc_method, free_dna_method=free_dna_method)


def create_cgna_config(
    parameter_set_name: str = 'Di_hmethyl_methylated-hemi_combine',
    group_split: bool = True
) -> CgnaConfig:
    """
    Convenience function to create CGNA configuration.
    
    Parameters
    ----------
    parameter_set_name : str, optional
        CGNA+ parameter set name
    group_split : bool, optional
        Use group splitting (default: True)
    
    Returns
    -------
    CgnaConfig
        Configuration object
    
    Examples
    --------
    >>> config = create_cgna_config()
    >>> config = create_cgna_config('cgnaplus_bsc1')
    """
    return CgnaConfig(parameter_set_name=parameter_set_name, group_split=group_split)


# Pre-defined common configurations
DEFAULT_RBP_CONFIG = RBPConfig(nuc_method='crystal')
DEFAULT_CGNA_CONFIG = CgnaConfig(parameter_set_name='Di_hmethyl_methylated-hemi_combine')
MULTIHARMONIC_RBP_CONFIG = RBPConfig(nuc_method='crystal', free_dna_method='free_dna')
