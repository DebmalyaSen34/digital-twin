"""
Dynamic FBA Process
Runs COBRApy FBA with reaction bounds constrained by enzyme activities
from the kinetics process.

This bridges the kinetic layer with the constraint-based metabolic model.
"""

import cobra
from vivarium.core.process import Process


class DynamicFBA(Process):
    """
    Constraint-based FBA that reads enzyme activity states and adjusts
    reaction bounds accordingly before optimizing.

    For each gene target with reduced activity:
      - All reactions associated with that gene have bounds scaled by activity
      - reaction.bounds = (lb * activity, ub * activity)
    """

    defaults = {
        "model_path": "src/iPS189.xml",
        "objective": "Biomass",
        "time_step": 1.0,
        "gene_reaction_map": {},   # gene_id -> [reaction_ids]
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)
        # Load model once
        self.model = cobra.io.read_sbml_model(self.parameters["model_path"])
        self.model.objective = self.parameters["objective"]

        # Build gene→reaction map if not provided
        if not self.parameters["gene_reaction_map"]:
            grm = {}
            for gene in self.model.genes:
                grm[gene.id] = [rxn.id for rxn in gene.reactions]
            self.parameters["gene_reaction_map"] = grm

        # Store original bounds for each reaction (so we can restore/scale)
        self._original_bounds = {}
        for rxn in self.model.reactions:
            self._original_bounds[rxn.id] = (rxn.lower_bound, rxn.upper_bound)

    def ports_schema(self):
        schema = {
            "enzymes": {},
            "fluxes": {},
            "global": {
                "growth_rate": {
                    "_default": 0.0,
                    "_updater": "set",
                    "_emit": True,
                }
            },
        }

        # Enzyme activity inputs (one per gene)
        for gene in self.model.genes:
            schema["enzymes"][f"{gene.id}_activity"] = {
                "_default": 1.0,
                "_updater": "set",
                "_emit": True,
            }

        # Flux outputs (one per reaction)
        for rxn in self.model.reactions:
            schema["fluxes"][rxn.id] = {
                "_default": 0.0,
                "_updater": "set",
                "_emit": True,
            }

        return schema

    def next_update(self, timestep, states):
        grm = self.parameters["gene_reaction_map"]

        # Reset all bounds to original
        for rxn in self.model.reactions:
            orig_lb, orig_ub = self._original_bounds[rxn.id]
            rxn.lower_bound = orig_lb
            rxn.upper_bound = orig_ub

        # Apply enzyme activity constraints
        for gene_id, rxn_ids in grm.items():
            activity_key = f"{gene_id}_activity"
            activity = states["enzymes"].get(activity_key, 1.0)

            if activity < 1.0:
                for rxn_id in rxn_ids:
                    try:
                        rxn = self.model.reactions.get_by_id(rxn_id)
                        orig_lb, orig_ub = self._original_bounds[rxn_id]
                        rxn.lower_bound = orig_lb * activity
                        rxn.upper_bound = orig_ub * activity
                    except KeyError:
                        pass

        # Run FBA
        solution = self.model.optimize()

        update = {"fluxes": {}, "global": {}}

        if solution.status == "optimal":
            update["global"]["growth_rate"] = float(solution.objective_value)
            for rxn in self.model.reactions:
                update["fluxes"][rxn.id] = float(solution.fluxes.get(rxn.id, 0.0))
        else:
            update["global"]["growth_rate"] = 0.0
            for rxn in self.model.reactions:
                update["fluxes"][rxn.id] = 0.0

        return update