"""
Whole-Cell Composite
Composes the three processes (diffusion, kinetics, FBA) into a single
Vivarium simulation that can be stepped forward in time.
"""

from vivarium.core.composer import Composer
from vivarium.core.engine import Engine

from src.vivarium_model.drug_diffusion import DrugDiffusion
from src.vivarium_model.binding_kinetics import EnzymeKinetics, DRUG_TARGET_PARAMS
from src.vivarium_model.fba_process import DynamicFBA

class WholeCellComposite(Composer):
    """
    Composes:
      1. DrugDiffusion  — extracellular → cytoplasm drug transport
      2. EnzymeKinetics — drug-enzyme binding → fractional activity
      3. DynamicFBA     — enzyme activities → constrained FBA → fluxes + growth

    Wiring:
      DrugDiffusion.cytoplasm.drug_concentration
        → EnzymeKinetics.cytoplasm.drug_concentration
      EnzymeKinetics.enzymes.*_activity
        → DynamicFBA.enzymes.*_activity
    """

    defaults = {
        "drug_name": "Trimethoprim",
        "model_path": "src/iPS189.xml",
        "diffusion_timestep": 0.1,   # seconds (fast process)
        "kinetics_timestep": 1.0,    # seconds
        "fba_timestep": 10.0,        # seconds (slow, expensive)
    }

    def generate_processes(self, config):
        return {
            "diffusion": DrugDiffusion({
                "drug_name": config["drug_name"],
                "time_step": config["diffusion_timestep"],
            }),
            "kinetics": EnzymeKinetics({
                "drug_name": config["drug_name"],
                "time_step": config["kinetics_timestep"],
            }),
            "fba": DynamicFBA({
                "model_path": config["model_path"],
                "time_step": config["fba_timestep"],
            }),
        }

    def generate_topology(self, config):
        return {
            "diffusion": {
                "extracellular": ("extracellular",),
                "cytoplasm": ("cytoplasm",),
                "membrane": ("membrane",),
            },
            "kinetics": {
                "cytoplasm": ("cytoplasm",),
                "enzymes": ("enzymes",),
            },
            "fba": {
                "enzymes": ("enzymes",),
                "fluxes": ("fluxes",),
                "global": ("global",),
            },
        }


def run_hybrid_simulation(
    drug_name: str,
    extracellular_concentration: float,  # µM
    total_time: float = 300.0,           # seconds (5 minutes)
    emit_step: float = 10.0,             # seconds between data points
) -> dict:
    """
    Run the full hybrid simulation and return time-series data.

    Returns:
        dict with keys:
          - time: list of time points
          - drug_intracellular: list of intracellular drug conc at each time
          - enzyme_activities: {gene_id: [activity at each time]}
          - growth_rate: list of growth rates at each time
          - fluxes: {rxn_id: [flux at each time]}  (only for affected reactions)
    """
    if drug_name == "Control (No Treatment)":
        # No drug — just run FBA once
        import cobra
        model = cobra.io.read_sbml_model("src/iPS189.xml")
        model.objective = "Biomass"
        sol = model.optimize()
        return {
            "time": [0.0],
            "drug_intracellular": [0.0],
            "enzyme_activities": {},
            "growth_rate": [float(sol.objective_value) if sol.status == "optimal" else 0.0],
            "fluxes": {rxn.id: [float(sol.fluxes.get(rxn.id, 0))] for rxn in model.reactions},
            "wt_growth": float(sol.objective_value) if sol.status == "optimal" else 0.0,
        }

    # Build composite
    composite = WholeCellComposite({
        "drug_name": drug_name,
    })

    # Initial state
    initial_state = {
        "extracellular": {
            "drug_concentration": extracellular_concentration,
        },
        "cytoplasm": {
            "drug_concentration": 0.0,
        },
        "membrane": {
            "permeability_factor": 1.0,
        },
        "enzymes": {},
        "fluxes": {},
        "global": {
            "growth_rate": 0.0,
        },
    }

    # Initialize enzyme activities to 1.0
    target_params = DRUG_TARGET_PARAMS.get(drug_name, {})
    for gene_id in target_params:
        initial_state["enzymes"][f"{gene_id}_activity"] = 1.0
        if target_params[gene_id].get("inhibition_type") == "irreversible":
            initial_state["enzymes"][f"{gene_id}_active_fraction"] = 1.0

    # Also run WT FBA for comparison
    import cobra
    wt_model = cobra.io.read_sbml_model("src/iPS189.xml")
    wt_model.objective = "Biomass"
    wt_sol = wt_model.optimize()
    wt_growth = float(wt_sol.objective_value) if wt_sol.status == "optimal" else 0.0
    wt_fluxes = {rxn.id: float(wt_sol.fluxes.get(rxn.id, 0)) for rxn in wt_model.reactions}

    # Create and run engine
    composite_state = composite.generate()
    engine = Engine(
        processes=composite_state["processes"],
        topology=composite_state["topology"],
        initial_state=initial_state,
        emitter="timeseries",
    )

    engine.update(total_time)
    timeseries = engine.emitter.get_timeseries()

    # Extract results
    times = timeseries.get("time", [0.0])
    drug_intra = []
    enzyme_acts = {gid: [] for gid in target_params}
    growth_rates = []
    flux_series = {}

    for t_idx, t in enumerate(times):
        # Drug concentration
        cyto = timeseries.get(("cytoplasm", "drug_concentration"), [0.0])
        drug_intra.append(cyto[t_idx] if t_idx < len(cyto) else 0.0)

        # Enzyme activities
        for gid in target_params:
            key = ("enzymes", f"{gid}_activity")
            vals = timeseries.get(key, [1.0])
            enzyme_acts[gid].append(vals[t_idx] if t_idx < len(vals) else 1.0)

        # Growth rate
        gr_vals = timeseries.get(("global", "growth_rate"), [0.0])
        growth_rates.append(gr_vals[t_idx] if t_idx < len(gr_vals) else 0.0)

    # Get final fluxes
    final_fluxes = {}
    for rxn_id in wt_fluxes:
        key = ("fluxes", rxn_id)
        vals = timeseries.get(key, [0.0])
        final_fluxes[rxn_id] = vals[-1] if vals else 0.0

    return {
        "time": times,
        "drug_intracellular": drug_intra,
        "enzyme_activities": enzyme_acts,
        "growth_rate": growth_rates,
        "fluxes": final_fluxes,
        "wt_growth": wt_growth,
        "wt_fluxes": wt_fluxes,
    }