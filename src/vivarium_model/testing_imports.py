import os
import sys

# Add src to path for process imports
project_root = os.path.join(os.path.dirname(__file__), "../..")
sys.path.insert(0, project_root)

try:
    from src.vivarium_model.binding_kinetics import DRUG_TARGET_PARAMS
    from src.vivarium_model.drug_diffusion import DRUG_PERMEABILITY
    from src.vivarium_model.whole_cell_composite import run_hybrid_simulation
    print("Successfully imported Vivarium model components.")
except ImportError as e:
    print(f"Error importing Vivarium model components: {e}")
