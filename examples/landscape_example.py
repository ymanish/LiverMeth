#!/usr/bin/env python3
"""
Simple example script for nucleosome free energy landscape analysis.

This script demonstrates how to use the landscape module to analyze
DNA sequences and calculate nucleosome binding free energies.
"""

from pathlib import Path
from src.landscape import LandscapeAnalyzer, calculate_landscape


def example_simple():
    """Simple example using the convenience function."""
    print("\n" + "="*60)
    print("EXAMPLE 1: Quick analysis with convenience function")
    print("="*60)
    
    # Define a test sequence (601 positioning sequence)
    seq_601 = (
        "CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACG"
        "TACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATA"
        "TATACATCCTGT"
    )
    
    print(f"\nSequence length: {len(seq_601)} bp")
    print(f"First 50 bp: {seq_601[:50]}...")
    
    # Quick analysis with default parameters
    results, metadata = calculate_landscape(
        sequence=seq_601,
        seq_id="seq601_test",
        window_size=147,
        step_size=10,  # Calculate every 10 bp
        store_sequences=False,  # Save memory by not storing sequences
        verbose=True,
    )
    
    print(f"\n[Results]")
    print(f"  Total positions analyzed: {len(results)}")
    print(f"  Sequence stored in results: {results[0].sequence is not None}")
    
    # Show first few results
    print(f"\n[First 3 results]")
    for i in range(min(3, len(results))):
        r = results[i]
        print(f"  Position {r.subid}: F={r.F:.2f}, dF={r.dF:.2f}, dyad={r.dyad_position}")
    
    return results, metadata


def example_detailed():
    """Detailed example using the LandscapeAnalyzer class."""
    print("\n" + "="*60)
    print("EXAMPLE 2: Detailed analysis with LandscapeAnalyzer")
    print("="*60)
    
    # Create analyzer with custom settings
    analyzer = LandscapeAnalyzer(
        window_size=147,
        step_size=1,  # Sliding by 1 bp
        n_workers=4,   # Use 4 parallel processes
    )
    
    # Define a synthetic AT-rich sequence
    test_sequence = "AT" * 100  # 200 bp AT-rich sequence
    
    print(f"\nSequence: {'AT'*25}... (AT-rich, 200 bp)")
    
    # Calculate landscape with specific binding parameters
    results, metadata = analyzer.calculate_landscape(
        sequence=test_sequence,
        seq_id="AT_rich",
        left=0,
        right=13,
        kresc_factor=1.0,
        store_sequences=True,  # Store sequences this time
        verbose=True,
    )
    
    print(f"\n[Metadata]")
    print(f"  ID: {metadata.id}")
    print(f"  Length: {metadata.length} bp")
    print(f"  Subsequences: {metadata.num_subsequences}")
    print(f"  Window: {metadata.window_size} bp")
    print(f"  Step: {metadata.step_size} bp")
    
    # Convert to DataFrame
    df = analyzer.results_to_dataframe(results, metadata)
    print(f"\n[DataFrame shape]: {df.shape}")
    print(f"[DataFrame columns]: {list(df.columns)}")
    print(f"\n[First few rows]:")
    print(df.head())
    
    # Save results
    output_dir = Path("example_output")
    saved_paths = analyzer.save_results(
        results, 
        metadata, 
        output_dir,
        save_pickle=True,
        save_csv=True
    )
    
    print(f"\n[Saved files]:")
    for format_type, path in saved_paths.items():
        print(f"  {format_type}: {path}")
    
    return results, metadata, df


def example_multiple_sequences():
    """Example analyzing multiple sequences."""
    print("\n" + "="*60)
    print("EXAMPLE 3: Multiple sequence analysis")
    print("="*60)
    
    # Create analyzer
    analyzer = LandscapeAnalyzer(window_size=147, step_size=20)
    
    # Define multiple test sequences
    sequences = [
        "CG" * 90,  # GC-rich
        "AT" * 90,  # AT-rich
    ]
    seq_ids = ["GC_rich", "AT_rich"]
    
    print(f"\nAnalyzing {len(sequences)} sequences:")
    for sid, seq in zip(seq_ids, sequences):
        print(f"  - {sid}: {len(seq)} bp")
    
    # Analyze all sequences
    all_results = analyzer.calculate_multiple_landscapes(
        sequences=sequences,
        seq_ids=seq_ids,
        store_sequences=False,
        verbose=True,
    )
    
    print(f"\n[Results Summary]")
    for (results, metadata) in all_results:
        print(f"  {metadata.id}: {len(results)} positions analyzed")
    
    # Save all results
    output_dir = Path("example_output")
    saved_info = analyzer.save_multiple_results(
        all_results,
        output_dir,
        combined_csv=True
    )
    
    print(f"\n[Saved files]")
    print(f"  Sequences: {saved_info['num_sequences']}")
    print(f"  Total positions: {saved_info['total_subsequences']}")
    if 'combined_csv' in saved_info:
        print(f"  Combined CSV: {saved_info['combined_csv']}")
    
    return all_results


def example_loading_results():
    """Example of loading previously saved results."""
    print("\n" + "="*60)
    print("EXAMPLE 4: Loading saved results")
    print("="*60)
    
    # Load results from pickle file
    pickle_path = Path("example_output/AT_rich_landscape.pkl")
    
    if pickle_path.exists():
        results, metadata = LandscapeAnalyzer().load_results(pickle_path)
        
        print(f"\n[Loaded from]: {pickle_path}")
        print(f"  ID: {metadata.id}")
        print(f"  Positions: {len(results)}")
        print(f"  First result: F={results[0].F:.2f}, subid={results[0].subid}")
        
        return results, metadata
    else:
        print(f"\n[File not found]: {pickle_path}")
        print("  Run example_detailed() first to create the file.")
        return None, None


def main():
    """Run all examples."""
    print("\n" + "="*70)
    print(" Nucleosome Free Energy Landscape Analysis - Examples")
    print("="*70)
    
    # Run examples
    try:
        # Example 1: Quick analysis
        results1, meta1 = example_simple()
        
        # Example 2: Detailed analysis
        results2, meta2, df2 = example_detailed()
        
        # Example 3: Multiple sequences
        all_results = example_multiple_sequences()
        
        # Example 4: Load results
        loaded_results, loaded_meta = example_loading_results()
        
        print("\n" + "="*70)
        print(" All examples completed successfully!")
        print("="*70)
        print("\nCheck the 'example_output' directory for saved files.")
        
    except Exception as e:
        print(f"\n[ERROR] An error occurred: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
