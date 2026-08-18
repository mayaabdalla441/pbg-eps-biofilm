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
        # First-order EPS turnover/shedding rate (1/h), NOT a literature-derived value -- no
        # citable rate constant was found for E. coli matrix turnover or shear-driven erosion
        # despite a real search; treated as an explicit assumption, like growth_floor_fraction.
        # Default 0.0 = no-op, exactly reproduces prior validated behavior (unbounded EPS
        # accumulation). See flow_phase_open_issues memory item C.7.
        "eps_shed_rate": {"_type": "float", "_default": 0.0},
        # Flow-phase Stage 1 scaffold (flow_phase_open_issues memory items D.8/D.9) -- NOT
        # calibrated values. These validate the MECHANICS (depletion-bypass transition, growth-
        # cap numerical safety) independent of picking real numbers, which are still blocked on
        # Nona/literature (see run_flow_scaffold.py for the full rationale). Defaults are no-ops:
        # flow never starts (1e9h) and the carrying capacity is effectively infinite (1e9),
        # exactly reproducing prior validated static-phase behavior.
        "flow_start_time_h": {"_type": "float", "_default": 1e9},
        "biomass_carrying_capacity": {"_type": "float", "_default": 1e9},
        # QUICK, EXPLICITLY-APPROXIMATE stress test of the diffusion-limitation hypothesis
        # (flow_phase_open_issues memory item, 2026-08-04): the well-mixed model has no spatial
        # structure, so this is a phenomenological stand-in, NOT the real mechanism (that needs
        # viva-biofilm's diff_biofilm support -- currently a stub, "unused until Phase B" per
        # world.rs -- ask sent to the PI, response pending). During flow, amino-acid uptake is
        # throttled by a factor that decays from 1.0 (no limitation) toward
        # diffusion_limitation_floor as accumulated EPS grows past diffusion_limitation_scale,
        # approximating reduced nutrient diffusion into a maturing, matrix-dense biofilm.
        # diffusion_limitation_floor=0.25 matches Stewart 2003's organic-solute D_e/D_aq ratio;
        # diffusion_limitation_scale (mg) has NO data/literature anchor -- swept, not fitted.
        # Defaults (floor=1.0, scale=1e9) are a no-op, reproducing prior validated behavior.
        "diffusion_limitation_scale": {"_type": "float", "_default": 1e9},
        "diffusion_limitation_floor": {"_type": "float", "_default": 1.0},
        # Volume-anchoring fix (decline_investigation_summary memory, 2026-08-17): the static/
        # depleting branch's uptake cap (`remaining / (biomass * interval)`) mixes AA
        # concentration (mM = mmol/L) with total biomass (gDW) with NO explicit shared reference
        # volume -- dimensionally this only makes sense if the "pool" is implicitly 1 L, which
        # was never stated or chosen deliberately; it's a units gap that predates this
        # investigation (present since the original static-phase model). Root-caused via the
        # depth-scaffold negative result: because of this gap, real depletion stays negligible
        # regardless of biomass scale, so production never plateaus and no shedding rate can ever
        # catch up to it. Default 1.0 (L) is a DELIBERATE no-op -- it exactly reproduces the prior
        # (buggy-but-already-validated/cited) formula bit-for-bit, so PR #1's static-phase
        # headline numbers and every other existing validated run are UNCHANGED unless a caller
        # opts in explicitly. Set to a real physical volume (e.g. Nona's ~3uL chamber = 3e-6 L,
        # or chamber_volume/num_layers per depth layer) to anchor uptake/depletion to reality.
        "substrate_reference_volume_L": {"_type": "float", "_default": 1.0},
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
            # Read back the current accumulated EPS mass so the shedding term below has
            # something to decay -- eps was previously write-only (accumulate-type output).
            "eps": {"_type": "mass", "_units": "mg"},
            # Framework-maintained simulation clock, read to decide static- vs flow-phase
            # behavior (see flow_start_time_h below).
            "global_time": {"_type": "float", "_default": 0.0},
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
        eps_current = inputs.get("eps", 0.0)
        global_time = inputs.get("global_time", 0.0)

        # Flow-phase scaffold: once flow starts, treat amino-acid concentration as pinned near
        # the inlet value rather than depleting -- justified by Nona's real chamber turnover
        # (~36sec for 5uL/min into a 3uL chamber, ~100 chamber-volumes/hr), fast enough that
        # replenishment swamps consumption. Uptake is still bounded by max_rate (the reaction's
        # own kinetic ceiling), just not by what's "left" in a pool that isn't meaningfully
        # depleting. See flow_phase_open_issues memory item D.8.
        flowing = global_time >= self.config["flow_start_time_h"]

        # Diffusion-limitation stress test (see config_schema comment) -- decays from 1.0 toward
        # diffusion_limitation_floor as eps_current grows past diffusion_limitation_scale.
        floor = self.config["diffusion_limitation_floor"]
        diffusion_factor = floor + (1.0 - floor) * np.exp(-eps_current / self.config["diffusion_limitation_scale"])

        with self.model:
            for aa_id, max_rate in self.config['aa_bounds'].items():
                if flowing:
                    cap = max_rate * diffusion_factor
                else:
                    remaining = substrates.get(aa_id, 0.0)
                    # remaining (mM=mmol/L) * volume_L -> total mmol available in the real pool;
                    # dividing THAT by (biomass * interval) gives mmol/gDW/h, matching max_rate's
                    # units. volume_L=1.0 (default) reproduces the old formula exactly.
                    volume_l = self.config["substrate_reference_volume_L"]
                    remaining_total_mmol = remaining * volume_l
                    cap = min(max_rate, remaining_total_mmol / (biomass * interval)) if biomass > 0 else max_rate
                    cap = max(cap, 0.0)
                self.model.reactions.get_by_id(aa_id).lower_bound = -cap

            #Stage #1: maximize growth
            self.model.objective = self.config['biomass_reaction_id']
            sol1 = self.model.optimize()
            if sol1.status != 'optimal':
                # Depletion caps make growth infeasible this step -- flux values from a
                # non-optimal solve aren't meaningful, don't trust them. No growth, no
                # consumption, no EPS PRODUCTION this step -- but shedding still applies (a
                # mechanical/enzymatic process, not dependent on this step's FBA solve). Found
                # via the eps_shed_rate sweep: hardcoding "eps": 0.0 here silently zeroed out
                # shedding for every step once growth goes infeasible -- exactly the post-
                # plateau regime the shedding mechanism exists to model. See flow_phase_open_
                # issues memory item C.7.
                return {
                    "substrates": {aa_id: 0.0 for aa_id in self.config['aa_bounds']},
                    "biomass": 0.0,
                    "eps": -self.config["eps_shed_rate"] * eps_current * interval,
                }
            mu = sol1.fluxes[self.config['biomass_reaction_id']]

            # Sensitivity-test knob (flow_phase_open_issues memory item A.1), NOT a validated
            # biological mechanism: scales down the reference growth rate used to set stage 2's
            # floor, to test sensitivity to the (unrealistic) assumption that cells grow at their
            # unconstrained planktonic optimum. Default 1.0 = no-op.
            mu_reference = mu * self.config['growth_rate_cap_fraction']

            # Flow-phase scaffold: logistic space/carrying-capacity constraint (flow_phase_open_
            # issues memory items D.8/D.9). Removing the nutrient-depletion limiter above (once
            # flowing) would otherwise let growth blow up unboundedly -- real biofilms plateau at
            # a mature, space/diffusion-limited size instead. Folded into mu_reference (feeding
            # the SAME stage-2 solve below) rather than applied as a post-hoc multiplier on the
            # output, so substrate consumption stays mass-balanced with the reduced growth --
            # exactly the class of bug fixed in flow_phase_open_issues item E. max(0.0, ...)
            # guards against a negative factor if a coarse step ever overshoots the capacity.
            space_factor = max(0.0, 1.0 - biomass / self.config["biomass_carrying_capacity"])
            mu_reference = mu_reference * space_factor

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
                # same "nothing trustworthy to report" case as an infeasible stage 1. Shedding
                # still applies here too (see the matching comment on stage 1's bailout above).
                return {
                    "substrates": {aa_id: 0.0 for aa_id in self.config['aa_bounds']},
                    "biomass": 0.0,
                    "eps": -self.config["eps_shed_rate"] * eps_current * interval,
                }
            mu_realized = sol2.fluxes[self.config['biomass_reaction_id']]
            aa_flux = {aa_id: sol2.fluxes.get(aa_id, 0.0) for aa_id in self.config['aa_bounds']}
            eps_flux = sol2.fluxes[self.config['eps_reaction_id']]

        # First-order shedding: d(EPS)/dt = production - eps_shed_rate * EPS. Applied here (not
        # folded into the FBA solve) since it represents mechanical/enzymatic matrix removal, not
        # a metabolic flux -- eps_shed_rate=0.0 (default) reproduces the prior unbounded-
        # accumulation behavior exactly.
        eps_produced = eps_flux * biomass * interval
        eps_shed = self.config["eps_shed_rate"] * eps_current * interval

        if flowing:
            # Pool is treated as pinned at inlet concentration during flow (see above) -- don't
            # report a depleting delta for it, or the tracked store would drain toward negative
            # values even though physically nothing is running out.
            substrate_deltas = {aa_id: 0.0 for aa_id in self.config['aa_bounds']}
        else:
            # Inverse of the cap-side fix above: flux*biomass*interval is a total mmol change;
            # dividing by volume_L converts it back to a concentration (mM) delta on the
            # `substrates` store, consistent with what the cap now assumes. volume_L=1.0
            # (default) reproduces the old formula exactly.
            volume_l = self.config["substrate_reference_volume_L"]
            substrate_deltas = {aa_id: (flux * biomass * interval) / volume_l
                                 for aa_id, flux in aa_flux.items()}

        return{
          "substrates": substrate_deltas,
          "biomass": mu_realized * biomass * interval,
          "eps": eps_produced - eps_shed,
        }

