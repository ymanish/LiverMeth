#!/usr/bin/env python3
"""
Example script for Nucleosome Free Energy calculations.

This script demonstrates how to use the modular nucleosome breath calculator
with configuration classes for various parameterization methods.
"""

import time
from femodules.nucleosome_breath_modular import NucleosomeBreathModular
from femodules.binding_sites import get_binding_info
from femodules.config import CgnaConfig, RBPConfig


def example_1_basic_calculation():
    """Example 1: Basic free energy calculation with CGNA+ config."""
    print("\n" + "="*70)
    print("EXAMPLE 1: Basic Free Energy Calculation (CGNA+)")
    print("="*70)
    
    # 601 positioning sequence with methylation
    seq_601 = (
        "CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGMTCTAGCACCGCTTAAACGCAMG"
        "TACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATA"
        "TATACATCCTGT"
    )
    
    print(f"\nSequence: {seq_601[:50]}... (147 bp)")
    print(f"          Contains methylated bases (M = 5mC)")
    
    # Create CGNA+ configuration (default methylation-aware)
    config = CgnaConfig()
    print(f"\n[Configuration]")
    print(f"  Type: CGNA+")
    print(f"  Parameter set: {config.parameter_set_name}")
    print(f"  Group split: {config.group_split}")
    
    # Create calculator
    nb = NucleosomeBreathModular(config)
    
    # Calculate free energy with full binding
    print("\n[Calculating with full binding (left=0, right=13)...]")
    start = time.perf_counter()
    result = nb.calculate_free_energy(
        sequence=seq_601,
        left=0,
        right=13,
        style="b_index"
    )
    elapsed = time.perf_counter() - start
    
    # Display results
    print(f"\n[Results]")
    print(f"  Total free energy (F):     {result.F:10.2f} kT")
    print(f"  Entropic contribution:     {result.F_entropy:10.2f} kT")
    print(f"  Enthalpic contribution:    {result.F_enthalpy:10.2f} kT")
    print(f"  Free DNA energy:           {result.F_freedna:10.2f} kT")
    print(f"  Delta_F:  {result.F - result.F_freedna:10.2f} kT")
    print(f"\n  Calculation time: {elapsed:.4f} seconds")


def example_2_partial_unwrapping():
    """Example 2: Partial unwrapping scenarios with custom CGNA config."""
    print("\n" + "="*70)
    print("EXAMPLE 2: Partial Unwrapping (Custom CGNA Config)")
    print("="*70)
    
    # Simple AT-rich sequence
    sequence = "AT" * 73 + "A"
    print(f"\nSequence: {'AT'*25}... (AT-rich, 147 bp)")
    
    # Create custom CGNA config
    config = CgnaConfig(parameter_set_name='Di_hmethyl_methylated-hemi_combine')
    print(f"\n[Using CGNA+ with {config.parameter_set_name}]")
    
    # Create calculator
    nb = NucleosomeBreathModular(config)
    
    # Test different unwrapping levels
    unwrap_configs = [
        (0, 13, "Fully wrapped"),
        (1, 12, "1 site open on each end"),
        (2, 11, "2 sites open on each end"),
        (3, 10, "3 sites open on each end"),
    ]
    
    print("\n[Testing different unwrapping configurations]")
    print(f"{'Config':<30} {'ΔF (kT)':<15} {'Info'}")
    print("-" * 70)
    
    for left, right, description in unwrap_configs:
        result = nb.calculate_free_energy(
            sequence=sequence,
            left=left,
            right=right,
            style="b_index"
        )
        dF = result.F - result.F_freedna
        
        # Get binding info
        info = get_binding_info(left, right, "b_index")
        
        print(f"{description:<30} {dF:>10.2f}     "
              f"{info['bound_sites']} bound, {info['unbound_sites']} open")




def example_3_rbp_config():
    """Example 3: Basic free energy calculation with RBP config."""
    print("\n" + "="*70)
    print("EXAMPLE 3: Basic Free Energy Calculation (RBP)")
    print("="*70)

    # Standard 601 sequence (A/C/G/T only for RBP workflow)
    seq_601 = (
        "CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGATCTAGCACCGCTTAAACGCAGG"
        "TACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATA"
        "TATACATCCTGT"
    )

    print(f"\nSequence: {seq_601[:50]}... (147 bp)")

    # Create RBP configuration
    config = RBPConfig(nuc_method='crystal')
    print(f"\n[Configuration]")
    print(f"  Type: RBP")
    print(f"  Nucleosome method: {config.nuc_method}")
    print(f"  Multiharmonic: {config.is_multiharmonic()}")

    # Create calculator
    nb = NucleosomeBreathModular(config)

    # Calculate free energy with full binding
    print("\n[Calculating with full binding (left=0, right=13)...]")
    start = time.perf_counter()
    result = nb.calculate_free_energy(
        sequence=seq_601,
        left=0,
        right=13,
        style="b_index"
    )
    elapsed = time.perf_counter() - start

    # Display results
    print(f"\n[Results]")
    print(f"  Total free energy (F):     {result.F:10.2f} kT")
    print(f"  Entropic contribution:     {result.F_entropy:10.2f} kT")
    print(f"  Enthalpic contribution:    {result.F_enthalpy:10.2f} kT")
    print(f"  Free DNA energy:           {result.F_freedna:10.2f} kT")
    print(f"  Delta_F:  {result.F - result.F_freedna:10.2f} kT")
    print(f"\n  Calculation time: {elapsed:.4f} seconds")


def example_6_scanning():
    """Example 6: Scanning left/right binding."""
    print("\n" + "="*70)
    print("EXAMPLE 6: Scanning Binding Configurations")
    print("="*70)
    
    sequence = "ACGT" * 36 + "ACG"
    print(f"\nSequence: {'ACGT'*12}... (147 bp)")
    print("\nScanning right-side binding (left=0, right=0→13)")
    
    nb = NucleosomeBreathModular()
    
    results = []
    for right_idx in range(0, 14):
        result = nb.calculate_free_energy(
            sequence=sequence,
            left=0,
            right=right_idx,
            style="b_index"
        )
        dF = result.F - result.F_freedna
        results.append((right_idx, dF))
    
    # Display scan results
    print("\n[Binding vs right-site index]")
    print(f"{'Right Index':<15} {'ΔF (kT)':<15} {'Bar Chart'}")
    print("-" * 70)
    
    for right_idx, dF in results:
        # Create simple bar chart
        bar_length = int(abs(dF) / 5)  # Scale for visualization
        bar = "█" * bar_length
        print(f"{right_idx:<15} {dF:>10.2f}     {bar}")


def main():
    """Run all examples."""
    print("\n" + "="*70)
    print(" Nucleosome Free Energy Calculator - Examples")
    print(" Modular Version")
    print("="*70)
    
    try:
        # Run all examples
        example_1_basic_calculation()
        example_2_partial_unwrapping()
        example_3_rbp_config()
        example_6_scanning()
        
        print("\n" + "="*70)
        print(" All examples completed successfully!")
        print("="*70)
        print("\nFor more information, see:")
        print("  - femodules/nucleosome_breath_modular.py")
        print("  - femodules/binding_sites.py")
        print("  - femodules/energy_calc.py")
        
    except Exception as e:
        print(f"\n[ERROR] An error occurred: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
