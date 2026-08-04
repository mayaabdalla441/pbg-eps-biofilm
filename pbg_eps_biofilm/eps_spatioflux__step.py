from typing import Dict, Any
import numpy as np
import cobra.io
from cobra.io import read_sbml_model

from process_bigraph import Composite, allocate_core
from spatio_flux.processes.dfba import DynamicFBA, MODEL_REGISTRY_DFBA

core = allocate_core()
print(core)
print(list(MODEL_REGISTRY_DFBA.keys()))

class EpsFBAStep(DynamicFBA):
    config_schema = {
        "model_file": "string",
        "biomass_reaction_id": {"_type": "string", "_default": "BIOMASS_Ec_iML1515_core_75p37M"},
        "eps_reaction_id": {"_type": "string", "_default": "EX_eps_e"},
        "growth_floor_fraction": {"_type": "float", "_default": 0.9},
        "aa_bounds": "map[float]",
        # Sensitivity-test knob, NOT a validated biological mechanism: scales the FBA-optimal
        # growth rate down by this fraction before the floor-then-maximize-EPS stage, to test
        # how sensitive time-to-depletion/total-EPS are to the (unrealistic) assumption that
        # cells grow at their unconstrained planktonic optimum. Default 1.0 = no-op, exactly
        # reproduces prior validated behavior. See flow_phase_open_issues memory item A.1.
        "growth_rate_cap_fraction": {"_type": "float", "_default": 1.0},
    }

    def initialize(self, config):
        self.model = read_sbml_model(config["model_file"])
        assert self.model.reactions.has_id(config["biomass_reaction_id"]), \
            f'{config["biomass_reaction_id"]} not found in model'
        assert self.model.reactions.has_id(config["eps_reaction_id"]), \
            f'{config["eps_reaction_id"]} not found in model'

    def inputs(self):
        return {
            "substrates": {"_type": "map", "_value": {"_type": "concentration", "_units": "mM"}},
            "biomass": {"_type": "mass", "_units": "gDW"},
        }
    def outputs(self):
        return {
            "substrates": "map[concentration_delta]",
            "biomass": "mass_delta",
            "eps": "mass_delta",
        }
    def update(self, inputs, interval):
        substrates = inputs["substrates"]
        biomass = inputs["biomass"]

        with self.model:
            for aa_id,max_rate in self.config['aa_bounds'].items():
                remaining = substrates.get(aa_id, 0.0)
                cap = min(max_rate, remaining / (biomass * interval)) if biomass > 0 else max_rate
                cap = max(cap, 0.0)
                self.model.reactions.get_by_id(aa_id).lower_bound = -cap

            #Stage #1: maximize growth
            self.model.objective = self.config['biomass_reaction_id']
            sol1 = self.model.optimize()
            if sol1.status != 'optimal':
                # Depletion caps make growth infeasible this step -- flux values from a
                # non-optimal solve aren't meaningful, don't trust them. No growth, no
                # consumption, no EPS this step.
                return {
                    "substrates": {aa_id: 0.0 for aa_id in self.config['aa_bounds']},
                    "biomass": 0.0,
                    "eps": 0.0,
                }
            mu = sol1.fluxes[self.config['biomass_reaction_id']]

            # Sensitivity-test knob (flow_phase_open_issues memory item A.1), NOT a validated
            # biological mechanism: scales down the reference growth rate used to set stage 2's
            # floor, to test sensitivity to the (unrealistic) assumption that cells grow at their
            # unconstrained planktonic optimum. Default 1.0 = no-op.
            mu_reference = mu * self.config['growth_rate_cap_fraction']

            #Stage #2: floor growth at growth_floor_fraction * mu_reference, max EPS[e].
            # Stage 1's mu and aa_flux are NOT used for the reported outputs below -- they were
            # only needed to compute the floor. Reporting them directly would mean claiming both
            # the unconstrained-growth flux solution AND the floored-growth EPS-maximizing flux
            # solution happened in the same step, which isn't mass-balanced: a single finite
            # nutrient pool can only support one flux distribution per step. Stage 2's own solve
            # is the one real, self-consistent answer -- growth, substrate uptake, and EPS must
            # all be read from sol2 (see flow_phase_open_issues memory item E). In practice sol2's
            # realized growth settles at (or essentially at) the floor, since growing beyond it
            # isn't rewarded by the EPS-maximizing objective.
            self.model.reactions.get_by_id(self.config["biomass_reaction_id"]).lower_bound = \
                self.config["growth_floor_fraction"] * mu_reference
            self.model.objective = self.config['eps_reaction_id']
            sol2 = self.model.optimize()
            if sol2.status != 'optimal':
                # Stage 2's constraint is a strict relaxation of stage 1's (a lower floor on an
                # already-feasible reaction), so this should not normally happen -- treat as the
                # same "nothing trustworthy to report" case as an infeasible stage 1.
                return {
                    "substrates": {aa_id: 0.0 for aa_id in self.config['aa_bounds']},
                    "biomass": 0.0,
                    "eps": 0.0,
                }
            mu_realized = sol2.fluxes[self.config['biomass_reaction_id']]
            aa_flux = {aa_id: sol2.fluxes.get(aa_id, 0.0) for aa_id in self.config['aa_bounds']}
            eps_flux = sol2.fluxes[self.config['eps_reaction_id']]

        return{
          "substrates": {aa_id: flux * biomass * interval for aa_id, flux in aa_flux.items()},
          "biomass": mu_realized * biomass * interval,
          "eps": eps_flux * biomass * interval,
        }

