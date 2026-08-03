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
            aa_flux = {aa_id: sol1.fluxes.get(aa_id, 0.0) for aa_id in self.config['aa_bounds']}

            #Stage #2: floor growth, max EPS[e]
            self.model.reactions.get_by_id(self.config["biomass_reaction_id"]).lower_bound = \
                self.config["growth_floor_fraction"] * mu
            self.model.objective = self.config['eps_reaction_id']
            sol2 = self.model.optimize()
            eps_flux = sol2.fluxes[self.config['eps_reaction_id']] if sol2.status == 'optimal' else 0.0

        return{
          "substrates": {aa_id: flux * biomass * interval for aa_id, flux in aa_flux.items()},
          "biomass": mu * biomass * interval,
          "eps": eps_flux * biomass * interval,

        }

