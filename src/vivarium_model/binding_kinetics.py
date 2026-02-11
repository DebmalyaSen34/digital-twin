"""
Enzyme Binding Kinetics Process
Models drug-enzyme interactions using competitive inhibition kinetics.

For each target enzyme:
  - Drug binds with association rate k_on and dissociates with k_off
  - Inhibition constant Ki = k_off / k_on
  - Fractional enzyme activity = 1 / (1 + [Drug_in] / Ki)  (competitive)
  - Or for irreversible: activity decays exponentially with drug exposure

The resulting enzyme activities are used to constrain FBA reaction bounds.
"""

import numpy as np
from vivarium.core.process import Process


# Drug-target binding parameters
# Ki values in µM (inhibition constants from literature)
DRUG_TARGET_PARAMS = {
    "Trimethoprim": {
        "MG228": {  # DHFR
            "Ki": 0.005,             # 5 nM — very tight binding
            "k_on": 1e6,            # 1/M/s — diffusion-limited
            "k_off": 5e-3,          # 1/s
            "inhibition_type": "competitive",
            "hill_coefficient": 1.0,
        },
    },
    "Methotrexate": {
        "MG228": {
            "Ki": 0.001,             # 1 nM — tighter than trimethoprim
            "k_on": 1e7,
            "k_off": 1e-2,
            "inhibition_type": "competitive",
            "hill_coefficient": 1.0,
        },
        "MG006": {
            "Ki": 0.05,
            "k_on": 5e5,
            "k_off": 2.5e-2,
            "inhibition_type": "competitive",
            "hill_coefficient": 1.0,
        },
    },
    "Fosmidomycin": {
        "MG066": {
            "Ki": 0.035,             # 35 nM
            "k_on": 8e5,
            "k_off": 2.8e-2,
            "inhibition_type": "slow_tight",
            "hill_coefficient": 1.0,
        },
    },
    "Cerulenin": {
        "MG212": {
            "Ki": 0.8,              # 0.8 µM
            "k_on": 2e4,
            "k_off": 1.6e-2,
            "inhibition_type": "irreversible",  # Cerulenin is a covalent inhibitor
            "hill_coefficient": 1.0,
            "inactivation_rate": 0.01,  # k_inact (1/s)
        },
        "MG114": {
            "Ki": 1.2,
            "k_on": 1.5e4,
            "k_off": 1.8e-2,
            "inhibition_type": "irreversible",
            "hill_coefficient": 1.0,
            "inactivation_rate": 0.008,
        },
    },
    "Mupirocin": {
        "MG345": {
            "Ki": 0.02,             # 20 nM — very potent
            "k_on": 5e6,
            "k_off": 1e-1,
            "inhibition_type": "competitive",
            "hill_coefficient": 1.0,
        },
    },
    "Generic Glycolysis Inhibitor": {
        "MG041": {
            "Ki": 5.0,
            "k_on": 1e4,
            "k_off": 5e-2,
            "inhibition_type": "competitive",
            "hill_coefficient": 1.5,  # Slight cooperativity
        },
        "MG429": {
            "Ki": 8.0,
            "k_on": 8e3,
            "k_off": 6.4e-2,
            "inhibition_type": "competitive",
            "hill_coefficient": 1.0,
        },
    },
}


class EnzymeKinetics(Process):
    """
    Computes fractional enzyme activity for each drug target based on
    intracellular drug concentration and binding kinetics.

    Supports:
      - Competitive inhibition: activity = 1 / (1 + ([I]/Ki)^n)
      - Irreversible inhibition: dE_active/dt = -k_inact * [I]/(Ki + [I]) * E_active
      - Slow-tight binding: tracks bound/unbound enzyme explicitly
    """

    defaults = {
        "drug_name": "Trimethoprim",
        "time_step": 1.0,
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)
        drug = self.parameters["drug_name"]
        self.target_params = DRUG_TARGET_PARAMS.get(drug, {})

    def ports_schema(self):
        schema = {
            "cytoplasm": {
                "drug_concentration": {
                    "_default": 0.0,
                    "_updater": "accumulate",
                    "_emit": True,
                }
            },
            "enzymes": {},
        }

        # Create a port for each target enzyme's activity
        for gene_id, params in self.target_params.items():
            schema["enzymes"][f"{gene_id}_activity"] = {
                "_default": 1.0,       # Fully active
                "_updater": "set",
                "_emit": True,
            }
            # For irreversible inhibitors, track remaining active enzyme fraction
            if params.get("inhibition_type") == "irreversible":
                schema["enzymes"][f"{gene_id}_active_fraction"] = {
                    "_default": 1.0,
                    "_updater": "set",
                    "_emit": True,
                }

        return schema

    def next_update(self, timestep, states):
        drug_conc = states["cytoplasm"]["drug_concentration"]
        update = {"enzymes": {}}

        for gene_id, params in self.target_params.items():
            Ki = params["Ki"]
            n = params.get("hill_coefficient", 1.0)
            inh_type = params.get("inhibition_type", "competitive")

            if inh_type == "competitive":
                # Standard competitive inhibition (Hill equation variant)
                if drug_conc > 0 and Ki > 0:
                    activity = 1.0 / (1.0 + (drug_conc / Ki) ** n)
                else:
                    activity = 1.0
                update["enzymes"][f"{gene_id}_activity"] = activity

            elif inh_type == "irreversible":
                # Irreversible (covalent) inhibition
                # Active fraction decreases over time with drug exposure
                k_inact = params.get("inactivation_rate", 0.01)
                current_fraction = states["enzymes"].get(f"{gene_id}_active_fraction", 1.0)

                if drug_conc > 0 and Ki > 0:
                    # Rate of inactivation: k_inact * [I] / (Ki + [I])
                    inactivation = k_inact * drug_conc / (Ki + drug_conc) * timestep
                    new_fraction = max(0.0, current_fraction - inactivation)
                else:
                    new_fraction = current_fraction  # No recovery for irreversible

                update["enzymes"][f"{gene_id}_active_fraction"] = new_fraction
                update["enzymes"][f"{gene_id}_activity"] = new_fraction

            elif inh_type == "slow_tight":
                # Slow tight-binding: kinetic approach to equilibrium
                k_on = params.get("k_on", 1e5)
                k_off = params.get("k_off", 1e-2)
                current_activity = states["enzymes"].get(f"{gene_id}_activity", 1.0)

                # dActivity/dt = k_off * (1 - activity) - k_on * [I] * activity
                # At equilibrium: activity_eq = k_off / (k_off + k_on * [I]) = Ki / (Ki + [I])
                d_activity = (k_off * (1 - current_activity) - k_on * drug_conc * 1e-6 * current_activity) * timestep
                new_activity = np.clip(current_activity + d_activity, 0.0, 1.0)
                update["enzymes"][f"{gene_id}_activity"] = new_activity

        return update