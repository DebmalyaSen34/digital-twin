"""
Whole-Cell Composite
Composes the three processes (diffusion, kinetics, FBA) into a single
Vivarium simulation that can be stepped forward in time.
"""

import numpy as np
import cobra

from src.vivarium_model.drug_diffusion import DrugDiffusion, DRUG_PERMEABILITY
from src.vivarium_model.binding_kinetics import EnzymeKinetics, DRUG_TARGET_PARAMS
from src.vivarium_model.fba_process import DynamicFBA

# Try Vivarium import — fall back to manual stepping if unavailable
VIVARIUM_ENGINE_AVAILABLE = False
try:
    from vivarium.core.composer import Composer
    from vivarium.core.engine import Engine
    VIVARIUM_ENGINE_AVAILABLE = True
except ImportError:
    pass


def _run_manual_simulation(
    drug_name: str,
    extracellular_concentration: float,
    total_time: float = 300.0,
    emit_step: float = 10.0,
) -> dict:
    """
    Manual time-stepping simulation that doesn't depend on Vivarium Engine.
    Steps: diffusion (fast) → kinetics → FBA (slow).
    This guarantees correct data flow between processes.
    """

    target_params = DRUG_TARGET_PARAMS.get(drug_name, {})

    # --- Get drug permeability parameters ---
    P = DRUG_PERMEABILITY.get(drug_name, 5.0e-6)
    av_ratio = 2.0e5       # M. genitalium surface-area-to-volume (1/cm)
    k_deg = 1.0e-4         # intracellular degradation rate (1/s)
    volume_ratio = 1000.0  # extracellular / intracellular volume

    # --- Load COBRA model once ---
    model = cobra.io.read_sbml_model("src/iPS189.xml")
    model.objective = "Biomass"

    # WT solution
    wt_sol = model.optimize()
    wt_growth = float(wt_sol.objective_value) if wt_sol.status == "optimal" else 0.0
    wt_fluxes = {}
    for rxn in model.reactions:
        wt_fluxes[rxn.id] = float(wt_sol.fluxes.get(rxn.id, 0.0))

    # Build gene → reaction map
    gene_rxn_map = {}
    for gene in model.genes:
        gene_rxn_map[gene.id] = [rxn.id for rxn in gene.reactions]

    # Store original bounds
    original_bounds = {}
    for rxn in model.reactions:
        original_bounds[rxn.id] = (rxn.lower_bound, rxn.upper_bound)

    # --- State variables ---
    c_out = float(extracellular_concentration)
    c_in = 0.0

    enzyme_activity = {}
    enzyme_active_fraction = {}  # for irreversible inhibitors
    for gid, params in target_params.items():
        enzyme_activity[gid] = 1.0
        if params.get("inhibition_type") == "irreversible":
            enzyme_active_fraction[gid] = 1.0

    # --- Time stepping ---
    dt_diffusion = 0.1      # seconds (fast)
    dt_kinetics = 1.0       # seconds
    dt_fba = max(emit_step, 10.0)  # seconds (slow)

    times = []
    drug_intra_series = []
    enzyme_act_series = {gid: [] for gid in target_params}
    growth_series = []

    current_growth = wt_growth
    current_fluxes = dict(wt_fluxes)

    t = 0.0
    next_emit = 0.0
    next_kinetics = 0.0
    next_fba = 0.0

    while t <= total_time + 1e-9:
        # === EMIT data at regular intervals ===
        if t >= next_emit - 1e-9:
            times.append(t)
            drug_intra_series.append(c_in)
            for gid in target_params:
                enzyme_act_series[gid].append(enzyme_activity[gid])
            growth_series.append(current_growth)
            next_emit += emit_step

        # === DIFFUSION (every dt_diffusion) ===
        dc_in = (P * av_ratio * (c_out - c_in) - k_deg * c_in) * dt_diffusion
        dc_out = -dc_in / volume_ratio
        c_in = max(0.0, c_in + dc_in)
        c_out = max(0.0, c_out + dc_out)

        # === KINETICS (every dt_kinetics) ===
        if t >= next_kinetics - 1e-9:
            for gid, params in target_params.items():
                Ki = params["Ki"]
                n = params.get("hill_coefficient", 1.0)
                inh_type = params.get("inhibition_type", "competitive")

                if inh_type == "competitive":
                    if c_in > 0 and Ki > 0:
                        enzyme_activity[gid] = 1.0 / (1.0 + (c_in / Ki) ** n)
                    else:
                        enzyme_activity[gid] = 1.0

                elif inh_type == "irreversible":
                    k_inact = params.get("inactivation_rate", 0.01)
                    frac = enzyme_active_fraction.get(gid, 1.0)
                    if c_in > 0 and Ki > 0:
                        inactivation = k_inact * c_in / (Ki + c_in) * dt_kinetics
                        frac = max(0.0, frac - inactivation)
                    enzyme_active_fraction[gid] = frac
                    enzyme_activity[gid] = frac

                elif inh_type == "slow_tight":
                    k_on = params.get("k_on", 1e5)
                    k_off = params.get("k_off", 1e-2)
                    curr = enzyme_activity[gid]
                    # dA/dt = k_off*(1-A) - k_on*[I]*A  (with [I] in M → µM*1e-6)
                    d_activity = (k_off * (1 - curr) - k_on * c_in * 1e-6 * curr) * dt_kinetics
                    enzyme_activity[gid] = float(np.clip(curr + d_activity, 0.0, 1.0))

            next_kinetics += dt_kinetics

        # === FBA (every dt_fba) ===
        if t >= next_fba - 1e-9:
            # Reset bounds
            for rxn in model.reactions:
                orig_lb, orig_ub = original_bounds[rxn.id]
                rxn.lower_bound = orig_lb
                rxn.upper_bound = orig_ub

            # Apply enzyme constraints
            for gid, activity in enzyme_activity.items():
                if activity < 1.0 - 1e-9:
                    for rxn_id in gene_rxn_map.get(gid, []):
                        try:
                            rxn = model.reactions.get_by_id(rxn_id)
                            orig_lb, orig_ub = original_bounds[rxn_id]
                            rxn.lower_bound = orig_lb * activity
                            rxn.upper_bound = orig_ub * activity
                        except KeyError:
                            pass

            sol = model.optimize()
            if sol.status == "optimal":
                current_growth = float(sol.objective_value)
                for rxn in model.reactions:
                    current_fluxes[rxn.id] = float(sol.fluxes.get(rxn.id, 0.0))
            else:
                current_growth = 0.0
                for rxn in model.reactions:
                    current_fluxes[rxn.id] = 0.0

            next_fba += dt_fba

        t += dt_diffusion

    return {
        "time": times,
        "drug_intracellular": drug_intra_series,
        "enzyme_activities": enzyme_act_series,
        "growth_rate": growth_series,
        "fluxes": current_fluxes,       # final snapshot
        "wt_growth": wt_growth,
        "wt_fluxes": wt_fluxes,
    }


def run_hybrid_simulation(
    drug_name: str,
    extracellular_concentration: float,
    total_time: float = 300.0,
    emit_step: float = 10.0,
) -> dict:
    """
    Run the full hybrid simulation and return time-series data.
    Uses manual stepping for reliability.
    """
    if drug_name == "Control (No Treatment)" or extracellular_concentration <= 0:
        model = cobra.io.read_sbml_model("src/iPS189.xml")
        model.objective = "Biomass"
        sol = model.optimize()
        gr = float(sol.objective_value) if sol.status == "optimal" else 0.0
        return {
            "time": [0.0],
            "drug_intracellular": [0.0],
            "enzyme_activities": {},
            "growth_rate": [gr],
            "fluxes": {rxn.id: float(sol.fluxes.get(rxn.id, 0)) for rxn in model.reactions},
            "wt_growth": gr,
            "wt_fluxes": {rxn.id: float(sol.fluxes.get(rxn.id, 0)) for rxn in model.reactions},
        }

    print(f"🔬 Running hybrid simulation: {drug_name} @ {extracellular_concentration} µM for {total_time}s")
    result = _run_manual_simulation(
        drug_name=drug_name,
        extracellular_concentration=extracellular_concentration,
        total_time=total_time,
        emit_step=emit_step,
    )

    # Debug summary
    if result["drug_intracellular"]:
        peak_drug = max(result["drug_intracellular"])
        print(f"  → Peak intracellular [Drug]: {peak_drug:.4f} µM")
    for gid, acts in result.get("enzyme_activities", {}).items():
        if acts:
            print(f"  → {gid} final activity: {acts[-1]:.4f}")
    if result["growth_rate"]:
        print(f"  → Final growth rate: {result['growth_rate'][-1]:.4f} (WT: {result['wt_growth']:.4f})")

    return result