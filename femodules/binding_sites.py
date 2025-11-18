"""
Binding site configuration and conversion utilities.

This module handles the phosphate binding site configurations on the nucleosome
and conversions between different binding specification styles.
"""

from typing import Tuple, List


# Crystal structure phosphate binding site positions (28 total)
PHOSPHATE_BINDING_SITES = [
    2, 6, 14, 17, 24, 29, 34, 38, 
    45, 49, 55, 59, 65, 69, 76, 
    80, 86, 90, 96, 100, 107, 111, 
    116, 121, 128, 131, 139, 143
]

# Total number of phosphate binding sites
N_PHOSPHATE_SITES = 28

# Number of bound sites (each covers 2 phosphates)
N_BOUND_SITES = 14


def select_phosphate_sites(left: int = 0, right: int = 13) -> List[int]:
    """
    Select a subset of phosphate binding sites.
    
    The nucleosome has 28 phosphate binding sites (14 on each strand).
    This function selects a subset based on bound-site indices.
    
    Parameters
    ----------
    left : int, optional
        Left bound-site index (0-13, default: 0)
    right : int, optional
        Right bound-site index (0-13, default: 13)
    
    Returns
    -------
    List[int]
        Indices of selected phosphate binding sites
    
    Examples
    --------
    >>> sites = select_phosphate_sites(0, 13)
    >>> len(sites)
    28
    >>> sites = select_phosphate_sites(1, 12)
    >>> len(sites)
    24
    
    Notes
    -----
    Each bound-site index corresponds to 2 phosphate positions.
    Total 28 phosphate sites = 14 bound sites × 2 phosphates per site.
    """
    return PHOSPHATE_BINDING_SITES[left*2:(right*2)+2]


def convert_to_open_sites(
    left: int,
    right: int,
    style: str = "b_index"
) -> Tuple[int, int]:
    """
    Convert binding indices to number of open phosphate sites.
    
    This function converts different representations of binding configurations
    into the number of open (unbound) phosphate sites on each end.
    
    Parameters
    ----------
    left : int
        Left specification (interpretation depends on style)
    right : int
        Right specification (interpretation depends on style)
    style : str, optional
        Specification style (default: 'b_index')
        
        Available styles:
        
        * 'b_index': Bound-site indices (0-13)
            - Each bound site covers 2 phosphates
            - Example: left=1, right=11 → l_open=2, r_open=4
            
        * 'ph_index': Phosphate-site indices (0-27)
            - Direct phosphate specification
            - Example: left=1, right=26 → l_open=1, r_open=1
            
        * 'open_sites': Direct count of open phosphates
            - No conversion needed
            - Example: left=2, right=4 → l_open=2, r_open=4
    
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
    >>> # Bound site indices
    >>> l, r = convert_to_open_sites(1, 11, 'b_index')
    >>> print(l, r)  # 2, 4
    
    >>> # Phosphate indices
    >>> l, r = convert_to_open_sites(2, 25, 'ph_index')
    >>> print(l, r)  # 2, 2
    
    >>> # Direct open sites
    >>> l, r = convert_to_open_sites(3, 5, 'open_sites')
    >>> print(l, r)  # 3, 5
    
    Notes
    -----
    Total phosphate sites: 28 (14 per strand)
    Bound-site i covers phosphates (2*i, 2*i+1)
    """
    if style == "b_index":
        # Bound-site indices (0-13)
        # Each bound site covers 2 phosphates
        l_open = 2 * left
        r_open = N_PHOSPHATE_SITES - (2 * right) - 2

    elif style == "ph_index":
        # Phosphate-site indices (0-27)
        # Direct phosphate specification
        l_open = left
        r_open = N_PHOSPHATE_SITES - right - 1
       
    elif style == "open_sites":
        # Direct specification of open phosphates
        # No conversion needed
        l_open = left
        r_open = right

    else:
        raise ValueError(
            f"Invalid style '{style}'. "
            "Use 'b_index', 'ph_index', or 'open_sites'."
        )
    
    return l_open, r_open


def get_binding_info(left: int, right: int, style: str = "b_index") -> dict:
    """
    Get complete binding configuration information.
    
    Parameters
    ----------
    left : int
        Left binding specification
    right : int
        Right binding specification
    style : str, optional
        Specification style (default: 'b_index')
    
    Returns
    -------
    dict
        Dictionary with binding configuration details:
        - 'left_open': Number of open sites on left
        - 'right_open': Number of open sites on right
        - 'bound_sites': Number of bound phosphate sites
        - 'unbound_sites': Number of unbound phosphate sites
        - 'binding_sites': List of bound site indices
    
    Examples
    --------
    >>> info = get_binding_info(0, 13, 'b_index')
    >>> info['bound_sites']
    28
    >>> info['left_open']
    0
    """
    l_open, r_open = convert_to_open_sites(left, right, style)
    bound_count = N_PHOSPHATE_SITES - l_open - r_open
    
    # Get the actual phosphate sites that are bound
    if style == "b_index":
        binding_sites = select_phosphate_sites(left, right)
    else:
        # Convert to bound indices first
        start_bound = l_open // 2
        end_bound = (N_PHOSPHATE_SITES - r_open) // 2 - 1
        binding_sites = select_phosphate_sites(start_bound, end_bound)
    
    return {
        'left_open': l_open,
        'right_open': r_open,
        'bound_sites': bound_count,
        'unbound_sites': l_open + r_open,
        'binding_sites': binding_sites,
        'total_sites': N_PHOSPHATE_SITES
    }
