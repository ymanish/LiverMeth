import os
# Import environment settings if flag is set (default: import)
if os.environ.get("IMPORT_ENV_SETTINGS", "1") == "1":
    from src.config.env_settings import *
import sys
from pathlib import Path

candidates = []
env_backend = os.environ.get("BACKEND_DIR")

if env_backend:                
    candidates.append(Path(env_backend))
candidates.extend([
    Path("/opt/backend"),
    Path(__file__).resolve().parent.parent.parent / "backend",
    Path("~/pol/Projects/Codebase/NucFreeEnergy").expanduser(),
])

backend_dir = next((p for p in candidates if p.is_dir()), None)
print(f"Using backend directory: {backend_dir}")
# Verify base directory exists
if not backend_dir.exists():
    raise FileNotFoundError(f"Base directory not found: {backend_dir}")

if str(backend_dir) not in sys.path:
    sys.path.insert(0, str(backend_dir))

# nucfree_energy = backend_dir / "NucFreeEnergy"
nucfree_energy = backend_dir
sys.path.insert(0, str(nucfree_energy))
print(f"Using NucFreeEnergy directory: {nucfree_energy}")

NUC_STATE_PATH = nucfree_energy / "methods" / "State" / "Nucleosome.state"
if not NUC_STATE_PATH.exists():
    raise FileNotFoundError(f"State file missing: {NUC_STATE_PATH}")

K_POSRESC_PATH = nucfree_energy / "MDParams" / "nuc_K_pos_resc_sym.npy"
if not K_POSRESC_PATH.exists():
    raise FileNotFoundError(f"Parameter file missing: {K_POSRESC_PATH}")

__all__ = ['NUC_STATE_PATH', 'K_POSRESC_PATH']