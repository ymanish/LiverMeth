# Landscape Module

Modular nucleosome free energy landscape analysis.

## Quick Start

```python
from src.landscape import calculate_landscape

# Quick analysis
results, meta = calculate_landscape(
    sequence="ACGT" * 200,
    seq_id="test",
    step_size=10,
    store_sequences=False
)

# Access results
for r in results[:5]:
    print(f"Position {r.subid}: F={r.F:.2f}, dF={r.dF:.2f}")
```

## Module Structure

- `__init__.py` - Public API and convenience functions
- `types.py` - Data structures (SubsequenceTask, LandscapeResult, SequenceMetadata)
- `tasks.py` - Subsequence generation from DNA sequences
- `workers.py` - Parallel energy calculation workers
- `analyzer.py` - LandscapeAnalyzer orchestration class
- `io_utils.py` - Save/load operations (pickle, CSV, DataFrame)

## Key Features

✨ **Simplified subid** - Uses integer start position instead of complex string  
💾 **Optional sequence storage** - Save memory with `store_sequences=False`  
🔧 **Modular design** - Clean separation of concerns  
📊 **Multiple formats** - Export to pickle, CSV, or DataFrame  
⚡ **Parallel processing** - Multi-core computation support  
📝 **Well documented** - Comprehensive docstrings and examples  

## Example

```python
from src.landscape import LandscapeAnalyzer

# Create analyzer
analyzer = LandscapeAnalyzer(
    window_size=147,
    step_size=10,
    n_workers=4
)

# Analyze sequence
results, meta = analyzer.calculate_landscape(
    sequence="ACGT" * 200,
    seq_id="my_seq",
    store_sequences=False,  # Save memory
    verbose=True
)

# Save results
analyzer.save_results(results, meta, "output")

# Convert to DataFrame
df = analyzer.results_to_dataframe(results, meta)
```

## Run Examples

```bash
python examples/landscape_example.py
```

## Documentation

- **Quick Reference**: `../docs/MODULAR_QUICK_REF.md`
- **Full Guide**: `../docs/MODULAR_LANDSCAPE.md`
- **Comparison**: `../docs/MODULAR_COMPARISON.md`
- **Summary**: `../docs/SUMMARY.md`

## Import

```python
# Main components
from src.landscape import LandscapeAnalyzer, calculate_landscape

# Data structures
from src.landscape import SubsequenceTask, LandscapeResult, SequenceMetadata

# Utilities
from src.landscape import results_to_dataframe, save_results, load_results
```

## Data Structure

```python
LandscapeResult(
    id='my_seq',          # Sequence identifier
    subid=42,             # Start position (int, not string!)
    sequence='ACGT...',   # Or None if store_sequences=False
    start_pos=42,
    end_pos=189,
    index=4,
    dyad_position=115,    # Center position
    F=-123.45,            # Total free energy
    F_entropy=-45.67,
    F_enthalpy=-77.78,
    F_freedna=-100.00,
    dF=-23.45,            # F - F_freedna
    left=0,
    right=13
)
```

## Performance

```python
# Memory efficient (don't store sequences)
results, meta = calculate_landscape(seq, "id", store_sequences=False)

# Faster scan (larger step)
results, meta = calculate_landscape(seq, "id", step_size=20)

# More workers (parallel)
analyzer = LandscapeAnalyzer(n_workers=8)
```
