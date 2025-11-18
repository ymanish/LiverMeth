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


def example_3_binding_strength():
    """Example 3: Effect of binding strength."""
    print("\n" + "="*70)
    print("EXAMPLE 3: Varying Binding Strength")
    print("="*70)
    
    # GC-rich sequence
    sequence = "CG" * 73 + "C"
    print(f"\nSequence: {'CG'*25}... (GC-rich, 147 bp)")
    
    # Create calculator
    nb = NucleosomeBreathModular()
    
    # Test different binding strengths
    kresc_factors = [0.5, 0.75, 1.0, 1.25, 1.5]
    
    print("\n[Effect of binding strength (kresc_factor)]")
    print(f"{'Factor':<15} {'ΔF (kT)':<15} {'Description'}")
    print("-" * 70)
    
    for factor in kresc_factors:
        result = nb.calculate_free_energy(
            sequence=sequence,
            left=0,
            right=13,
            kresc_factor=factor
        )
        dF = result.F - result.F_freedna
        
        if factor < 1.0:
            desc = "Weaker binding"
        elif factor > 1.0:
            desc = "Stronger binding"
        else:
            desc = "Normal binding"
        
        print(f"{factor:<15.2f} {dF:>10.2f}     {desc}")


def example_4_binding_styles():
    """Example 4: Different binding specification styles."""
    print("\n" + "="*70)
    print("EXAMPLE 4: Binding Specification Styles")
    print("="*70)
    
    sequence = "ACGT" * 36 + "ACG"
    print(f"\nSequence: {'ACGT'*12}... (147 bp)")
    
    nb = NucleosomeBreathModular()
    
    # Same configuration, different styles
    configs = [
        (0, 13, "b_index", "Bound-site indices (0-13)"),
        (0, 27, "ph_index", "Phosphate indices (0-27)"),
        (0, 0, "open_sites", "Direct open sites"),
    ]
    
    print("\n[Equivalent binding configurations in different styles]")
    print(f"{'Style':<15} {'Left':<8} {'Right':<8} {'ΔF (kT)':<15}")
    print("-" * 70)
    
    for left, right, style, description in configs:
        result = nb.calculate_free_energy(
            sequence=sequence,
            left=left,
            right=right,
            style=style
        )
        dF = result.F - result.F_freedna
        print(f"{style:<15} {left:<8} {right:<8} {dF:>10.2f}")
        print(f"  → {description}")


def example_5_batch_processing():
    """Example 5: Batch processing multiple sequences."""
    print("\n" + "="*70)
    print("EXAMPLE 5: Batch Processing")
    print("="*70)
    
    # Define multiple sequences
    sequences = {
        "AT-rich": "AT" * 73 + "A",
        "GC-rich": "CG" * 73 + "C",
        "Mixed-1": "ACGT" * 36 + "ACG",
        "Mixed-2": "TAGC" * 36 + "TAG",
    }
    
    print(f"\nProcessing {len(sequences)} sequences...")
    
    # Create calculator
    nb = NucleosomeBreathModular()
    
    # Calculate for all sequences
    print("\n[Results for all sequences]")
    print(f"{'Sequence':<15} {'F (kT)':<12} {'ΔF (kT)':<12} {'Time (s)'}")
    print("-" * 70)
    
    for name, seq in sequences.items():
        start = time.perf_counter()
        result = nb.calculate_free_energy(seq, 0, 13)
        elapsed = time.perf_counter() - start
        
        dF = result.F - result.F_freedna
        print(f"{name:<15} {result.F:>10.2f}  {dF:>10.2f}  {elapsed:>8.4f}")


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
        example_3_binding_strength()
        example_4_binding_styles()
        example_5_batch_processing()
        example_6_scanning()
        
        print("\n" + "="*70)
        print(" All examples completed successfully!")
        print("="*70)
        print("\nFor more information, see:")
        print("  - femodules/nucleosome_breath_modular.py")
        print("  - femodules/binding_sites.py")
        print("  - femodules/energy_calc.py")
        print("  - docs/CGNA_MODULE_ARCHITECTURE.md")
        
    except Exception as e:
        print(f"\n[ERROR] An error occurred: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
