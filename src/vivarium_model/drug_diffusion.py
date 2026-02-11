import numpy as np
from vivarium.core.process import Process


# Drug-specific permeability coefficients (cm/s) — literature estimates
DRUG_PERMEABILITY = {
    "Trimethoprim": 1.2e-5,       # Small, moderately lipophilic
    "Methotrexate": 2.0e-7,       # Charged, poor passive diffusion
    "Fosmidomycin": 5.0e-7,       # Hydrophilic, relies on transporters
    "Cerulenin": 8.0e-6,          # Lipophilic antibiotic
    "Mupirocin": 3.0e-6,          # Moderate lipophilicity
    "Generic Glycolysis Inhibitor": 5.0e-6,
}

# Drug molecular weights (g/mol) — affects membrane crossing
DRUG_MW = {
    "Trimethoprim": 290.3,
    "Methotrexate": 454.4,
    "Fosmidomycin": 183.1,
    "Cerulenin": 223.3,
    "Mupirocin": 500.6,
    "Generic Glycolysis Inhibitor": 200.0,
}


class DrugDiffusion(Process):
    """
    Simulates drug transport across the bacterial membrane.

    For M. genitalium:
      - Cell diameter ≈ 0.3 µm (sphere approximation)
      - Surface area A ≈ 0.283 µm² = 2.83e-9 cm²
      - Volume V ≈ 0.0141 µm³ = 1.41e-14 cm³
      - A/V ≈ 2.0e5 cm⁻¹
    """

    defaults = {
        "drug_name": "Trimethoprim",
        "permeability": None,          # Override or use lookup
        "surface_area_to_volume": 2.0e5,  # 1/cm, for M. genitalium
        "degradation_rate": 1e-4,      # 1/s, intracellular drug degradation
        "time_step": 1.0,              # seconds
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)
        drug = self.parameters["drug_name"]
        if self.parameters["permeability"] is None:
            self.parameters["permeability"] = DRUG_PERMEABILITY.get(drug, 5.0e-6)

    def ports_schema(self):
        return {
            "extracellular": {
                "drug_concentration": {
                    "_default": 0.0,
                    "_updater": "accumulate",
                    "_emit": True,
                    "_properties": {"units": "µM"},
                }
            },
            "cytoplasm": {
                "drug_concentration": {
                    "_default": 0.0,
                    "_updater": "accumulate",
                    "_emit": True,
                    "_properties": {"units": "µM"},
                }
            },
            "membrane": {
                "permeability_factor": {
                    "_default": 1.0,
                    "_updater": "set",
                    "_emit": True,
                }
            },
        }

    def next_update(self, timestep, states):
        P = self.parameters["permeability"]
        av_ratio = self.parameters["surface_area_to_volume"]
        k_deg = self.parameters["degradation_rate"]

        c_out = states["extracellular"]["drug_concentration"]
        c_in = states["cytoplasm"]["drug_concentration"]
        perm_factor = states["membrane"]["permeability_factor"]

        # Effective permeability (can be modulated by efflux pumps, etc.)
        P_eff = P * perm_factor

        # Fick's law diffusion + degradation
        dc_in = (P_eff * av_ratio * (c_out - c_in) - k_deg * c_in) * timestep

        # Drug is consumed from extracellular (conservation)
        # Volume ratio: extracellular >> intracellular, so negligible depletion
        # But we model it for completeness (assume 1000:1 volume ratio)
        volume_ratio = 1000.0
        dc_out = -dc_in / volume_ratio

        return {
            "extracellular": {"drug_concentration": dc_out},
            "cytoplasm": {"drug_concentration": dc_in},
        }