from cdFBA.processes.dfba import dFBA


class EpsCdFBAStep(dFBA):
    """Two-stage growth-then-EPS[e] dFBA Process, extending cdFBA's own dFBA class
    (not just borrowing a utility function from it).

    Reuses dFBA.__init__ as-is -- model loading via model_from_file, medium/bounds
    application, gene_knockout/reaction_knockout support (self.config["changes"]),
    and biomass_identifier auto-detection (self.biomass_identifier, via
    get_objective_reaction) all come from the parent for free. We only override
    inputs()/outputs()/update() for the two-stage EPS[e] logic, which doesn't exist
    in the parent at all -- and keep our own hard-cap amino-acid depletion instead
    of cdFBA's Michaelis-Menten kinetics (the "kinetics" field is inherited in
    config_schema but left unused/empty).
    """

    config_schema = {
        **dFBA.config_schema,
        "eps_reaction_id": {"_type": "string", "_default": "EX_eps_e"},
        "growth_floor_fraction": {"_type": "float", "_default": 0.9},
        "aa_bounds": "map[float]",
        # Sensitivity-test knob, kept in sync with EpsFBAStep -- see that file's comment and
        # flow_phase_open_issues memory item A.1. Default 1.0 = no-op.
        "growth_rate_cap_fraction": {"_type": "float", "_default": 1.0},
    }

    def inputs(self):
        return {
            "shared_environment": "volumetric",
        }

    def outputs(self):
        return {
            "shared_environment": "volumetric",
        }

    def update(self, inputs, interval):
        env = inputs["shared_environment"]
        concentrations = env["concentrations"]
        name = self.config["name"]
        biomass = concentrations.get(name, 0.0)

        with self.model:
            # Hard-cap each amino acid's uptake bound by what's actually left --
            # same depletion mechanism as EpsFBAStep, not cdFBA's own MM kinetics.
            for aa_id, max_rate in self.config["aa_bounds"].items():
                remaining = concentrations.get(aa_id, 0.0)
                cap = min(max_rate, remaining / (biomass * interval)) if biomass > 0 else max_rate
                cap = max(cap, 0.0)
                self.model.reactions.get_by_id(aa_id).lower_bound = -cap

            # Stage 1: maximize growth (reuse the parent's auto-detected biomass_identifier
            # instead of a redundant explicit config field)
            self.model.objective = self.biomass_identifier
            sol1 = self.model.optimize()
            if sol1.status != "optimal":
                # Depletion caps make growth infeasible this step -- don't trust flux
                # values from a non-optimal solve. No growth, no consumption, no EPS.
                zero_update = {aa_id: 0.0 for aa_id in self.config["aa_bounds"]}
                zero_update[name] = 0.0
                zero_update["eps"] = 0.0
                return {"shared_environment": {"counts": zero_update}}

            mu = sol1.fluxes[self.biomass_identifier]
            mu_reference = mu * self.config['growth_rate_cap_fraction']

            # Stage 2: floor growth at growth_floor_fraction * mu_reference, maximize EPS[e].
            # Stage 1's mu and aa_flux are NOT used for the reported outputs below -- only to
            # compute the floor. Reporting them directly would mean claiming both the
            # unconstrained-growth flux solution AND the floored-growth EPS-maximizing flux
            # solution happened in the same step, which isn't mass-balanced. Stage 2's own solve
            # is the one real, self-consistent answer (see flow_phase_open_issues memory item E,
            # and the matching fix/comment in EpsFBAStep).
            self.model.reactions.get_by_id(self.biomass_identifier).lower_bound = \
                self.config["growth_floor_fraction"] * mu_reference
            self.model.objective = self.config["eps_reaction_id"]
            sol2 = self.model.optimize()
            if sol2.status != "optimal":
                zero_update = {aa_id: 0.0 for aa_id in self.config["aa_bounds"]}
                zero_update[name] = 0.0
                zero_update["eps"] = 0.0
                return {"shared_environment": {"counts": zero_update}}

            mu_realized = sol2.fluxes[self.biomass_identifier]
            aa_flux = {aa_id: sol2.fluxes.get(aa_id, 0.0) for aa_id in self.config["aa_bounds"]}
            eps_flux = sol2.fluxes[self.config["eps_reaction_id"]]

        counts_update = {aa_id: flux * biomass * interval for aa_id, flux in aa_flux.items()}
        counts_update[name] = mu_realized * biomass * interval
        counts_update["eps"] = eps_flux * biomass * interval

        return {"shared_environment": {"counts": counts_update}}
