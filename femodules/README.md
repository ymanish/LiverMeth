# femodules

Modular nucleosome free-energy utilities for the LiverMeth project.

This package provides a refactored workflow for nucleosome binding free-energy calculations, with:

- **CGNA+** parameterization (`CgnaConfig`)
- **RBP** parameterization (`RBPConfig`)
- Binding-site conversion helpers
- Optional multiharmonic energy correction for RBP

> Note: `CgnaNucFreeEnergy.py` is a legacy module. Use `nucleosome_breath_modular.py` for new work.

## Files

- `nucleosome_breath_modular.py` — main calculator class (`NucleosomeBreathModular`)
- `config.py` — config dataclasses (`CgnaConfig`, `RBPConfig`)
- `binding_sites.py` — binding-site selection and style conversion utilities
- `energy_calc.py` — multiharmonic energy utilities
- `CgnaNucFreeEnergy.py` — legacy implementation (kept for compatibility)

## Prerequisites

This module depends on external backend resources (loaded in `femodules/__init__.py`):

- `methods/State/Nucleosome.state`
- `MDParams/nuc_K_pos_resc_sym.npy`

Backend resolution order:

1. `BACKEND_DIR` environment variable (recommended)
2. `/opt/backend`
3. `../backend` relative to this repository root
4. `~/pol/Projects/Codebase/NucFreeEnergy`

If none are found, import will fail with `FileNotFoundError`.

### Environment variables

- `BACKEND_DIR`: absolute path to backend directory with `methods/` and `MDParams/`
- `IMPORT_ENV_SETTINGS`:
  - `"1"` (default): import `src.config.env_settings`
  - `"0"`: skip importing environment settings

## Quick start

### 1) CGNA+ (default)

```python
from femodules.nucleosome_breath_modular import NucleosomeBreathModular
from femodules.config import CgnaConfig

config = CgnaConfig()  # default: Di_hmethyl_methylated-hemi_combine
nb = NucleosomeBreathModular(config)

sequence = "ACGT" * 36 + "ACG"  # 147 bp
result = nb.calculate_free_energy(sequence=sequence, left=0, right=13, style="b_index")

print(result)
print("ΔF =", result.F - result.F_freedna)
```

### 2) RBP

```python
from femodules.nucleosome_breath_modular import NucleosomeBreathModular
from femodules.config import RBPConfig

config = RBPConfig(nuc_method="crystal")
nb = NucleosomeBreathModular(config)

sequence = "ACGT" * 36 + "ACG"  # 147 bp
result = nb.calculate_free_energy(sequence=sequence, left=0, right=13, style="b_index")

print(result)
```

### 3) RBP multiharmonic

```python
from femodules.nucleosome_breath_modular import NucleosomeBreathModular
from femodules.config import RBPConfig

config = RBPConfig(nuc_method="crystal", free_dna_method="free_dna")
nb = NucleosomeBreathModular(config)

result = nb.calculate_free_energy(
    sequence="ACGT" * 36 + "ACG",
    left=1,
    right=12,
    style="b_index",
    bound_ends="exclude",  # or "include"
)
```

## Main API

### `NucleosomeBreathModular(config=None)`

- `config=None` uses default `CgnaConfig`
- Accepts only `RBPConfig` or `CgnaConfig`

### `calculate_free_energy(...)`

```python
calculate_free_energy(
    sequence: str,
    left: int,
    right: int,
    id: str | None = None,
    subid: str | None = None,
    kresc_factor: float = 1.0,
    style: str = "b_index",
    bound_ends: str = "exclude",
)
```

Returns `FreeEnergyResult` with:

- `F`
- `F_entropy`
- `F_enthalpy`
- `F_freedna`
- `id`, `subid`

### `calculate_free_energy_hard(...)`

Hard-binding model (RBP only). Raises `NotImplementedError` for CGNA config.

## Binding style options

The `style` argument supports:

- `"b_index"`: bound-site index space (`0..13`)
- `"ph_index"`: phosphate index space (`0..27`)
- `"open_sites"`: direct counts of open phosphates (left/right)

Helper utilities in `binding_sites.py`:

- `convert_to_open_sites(left, right, style)`
- `get_binding_info(left, right, style)`
- `select_phosphate_sites(left=0, right=13)`

## Running the example script

From project root:

```bash
python examples/femodules_example.py
```

This script demonstrates:

- CGNA+ full binding
- CGNA+ partial unwrapping
- RBP configuration example
- right-index scanning

## Notes

- Use sequences of nucleosome length (typically 147 bp) for expected behavior.
- CGNA+ workflows may use methylation-aware parameter sets (for example base symbol `M` in sequence strings).
- Keep `kresc_factor=1.0` unless you intentionally want to rescale binding strength.
