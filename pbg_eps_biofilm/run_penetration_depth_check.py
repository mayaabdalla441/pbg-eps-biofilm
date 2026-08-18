"""Step 1 of the depth-resolved spatial model plan (see
composites/eps_biofilm_depth_scaffold.py for the full design rationale): a
cheap, throwaway diagnostic answering ONE question before any spatial code
gets built --

Does a nutrient gradient across the real 150um chamber height plausibly
matter at this model's REALIZED (not enzymatic-ceiling) amino-acid uptake
flux?

Uses the standard zero-order-sink biofilm penetration-depth formula (the same
framing Stewart 2003 itself uses for O2/substrate penetration):

    L_p = sqrt(2 * D_eff * C_bulk / R)

where R = v_uptake (mmol/gDW/h) * B_density (gDW/cm^3) is the volumetric
consumption rate, D_eff = 0.25 * D_aq (Stewart 2003's organic-solute ratio,
already an accepted citation elsewhere in this codebase -- see
eps_spatioflux__step.py's diffusion_limitation_floor comment), and C_bulk is
each amino acid's real inlet concentration (DEFAULT_INITIAL_AA_CONC).

Does NOT read the model's reported `substrates` output (hard-zeroed during
flow by design, see EpsFBAStep.update() -- "Pool is treated as pinned at
inlet concentration during flow... don't report a depleting delta"). Instead
snoops on the underlying cobra Model.optimize() calls via a thin recording
subclass (ZERO changes to EpsFBAStep itself) to capture the REAL stage-2
(EPS-maximizing) solution's flux values, since stage 2's solve is the one
self-consistent, real answer per EpsFBAStep's own design (see its comments).

Guarded by `if __name__ == "__main__":` for the same reason as the other
run_*.py scripts.
"""
import warnings
warnings.filterwarnings('ignore')

import math

from process_bigraph import Composite
from process_bigraph.emitter import gather_emitter_results, emitter_from_wires

from pbg_eps_biofilm.core import build_core
from pbg_eps_biofilm.eps_spatioflux__step import EpsFBAStep
from pbg_eps_biofilm.composites.eps_biofilm_flow_scaffold import (
    eps_biofilm_flow_scaffold, DEFAULT_INITIAL_AA_CONC,
)

CHAMBER_VOLUME_CM3 = 0.2 * 1.0 * 0.015  # 2mm x 10mm x 150um = 0.003 cm^3
D_AQ_CM2_PER_S = 1e-5  # representative small-solute aqueous diffusivity
D_AQ_CM2_PER_H = D_AQ_CM2_PER_S * 3600.0
STEWART_ORGANIC_RATIO = 0.25  # Stewart 2003, D_e/D_aq for organic solutes
D_EFF_CM2_PER_H = STEWART_ORGANIC_RATIO * D_AQ_CM2_PER_H
CHAMBER_HEIGHT_UM = 150.0

CAPTURED = []  # list of {'time': h, 'biomass': gDW, 'aa_flux': {aa_id: mmol/gDW/h}}


class RecordingEpsFBAStep(EpsFBAStep):
    """Snoops on cobra's model.optimize() to capture the real stage-2 solve,
    without touching EpsFBAStep.update()'s logic at all."""

    def initialize(self, config):
        super().initialize(config)
        orig_optimize = self.model.optimize

        def recording_optimize(*a, **kw):
            sol = orig_optimize(*a, **kw)
            self._last_solutions.append(sol)
            return sol

        self._last_solutions = []
        self.model.optimize = recording_optimize

    def update(self, inputs, interval):
        self._last_solutions = []
        result = super().update(inputs, interval)
        if len(self._last_solutions) >= 2:
            sol2 = self._last_solutions[-1]  # stage 2 = EPS-maximizing solve, the real answer
            aa_flux = {aa_id: sol2.fluxes.get(aa_id, 0.0) for aa_id in self.config['aa_bounds']}
            CAPTURED.append({
                'time': inputs.get('global_time', 0.0),
                'biomass': inputs['biomass'],
                'aa_flux': aa_flux,
            })
        return result


def run_recording(total_hours=80.0, dt=0.1, **generator_kwargs):
    core = build_core()
    if 'RecordingEpsFBAStep' not in core.link_registry:
        core.register_link('RecordingEpsFBAStep', RecordingEpsFBAStep)

    state = eps_biofilm_flow_scaffold(dt=dt, **generator_kwargs)
    # Swap the composite's process address to use our recording subclass instead
    # of EpsFBAStep -- everything else (config, wiring) is unchanged.
    state['eps_fba']['address'] = 'local:RecordingEpsFBAStep'
    state['emitter'] = emitter_from_wires({
        "global_time": ["global_time"],
        "fields": ["fields"],
    })

    sim = Composite({"state": state}, core=core)
    sim.run(total_hours)
    return gather_emitter_results(sim)[('emitter',)]


def penetration_depth_um(v_uptake_mmol_gdw_h, biomass_gdw, c_bulk_mM):
    if v_uptake_mmol_gdw_h <= 0 or biomass_gdw <= 0:
        return float('inf')
    b_density_gdw_cm3 = biomass_gdw / CHAMBER_VOLUME_CM3
    r_mmol_cm3_h = v_uptake_mmol_gdw_h * b_density_gdw_cm3  # volumetric consumption rate
    c_bulk_mmol_cm3 = c_bulk_mM * 1e-3
    l_p_cm = math.sqrt(2 * D_EFF_CM2_PER_H * c_bulk_mmol_cm3 / r_mmol_cm3_h)
    return l_p_cm * 1e4  # cm -> um


def main():
    print(f"D_eff = {D_EFF_CM2_PER_H:.5g} cm^2/h (Stewart 2003 organic-solute ratio x D_aq)")
    print(f"Chamber volume = {CHAMBER_VOLUME_CM3:.5g} cm^3, height = {CHAMBER_HEIGHT_UM} um\n")

    print("Running real-calibrated flow scaffold with flux recording (dt=0.1h, coarse -- "
          "this is a sanity check, not the final calibrated run)...")
    run_recording(
        total_hours=80.0, dt=0.1,
        initial_biomass=2.4e-7, biomass_carrying_capacity=4.17e-4,
        growth_rate_cap_fraction=0.095, flow_start_time_h=3.0,
    )

    targets = [5.0, 40.0, 75.0]
    print(f"\n{'t (h)':>8} {'biomass (gDW)':>15} {'AA':>10} {'v_uptake':>12} {'L_p (um)':>10} {'vs 150um':>10}")
    for target_t in targets:
        closest = min(CAPTURED, key=lambda r: abs(r['time'] - target_t))
        biomass = closest['biomass']
        # Report the AA with the largest realized uptake flux (most likely to be
        # diffusion-limiting) at this timepoint.
        aa_flux = closest['aa_flux']
        top_aa = min(aa_flux, key=lambda k: aa_flux[k])  # most negative = largest uptake
        v_uptake = abs(aa_flux[top_aa])
        c_bulk = DEFAULT_INITIAL_AA_CONC[top_aa]
        l_p = penetration_depth_um(v_uptake, biomass, c_bulk)
        verdict = "<< 150um (gradient plausible)" if l_p < CHAMBER_HEIGHT_UM else ">> 150um (no real gradient)"
        print(f"{closest['time']:8.2f} {biomass:15.6g} {top_aa:>10} {v_uptake:12.4f} {l_p:10.2f} {verdict:>10}")


if __name__ == "__main__":
    # NOTE: invoke this via `python -c "import pbg_eps_biofilm.run_penetration_depth_check as m; m.main()"`
    # rather than `python -m ...` -- running as -m makes this file BOTH `__main__` AND get
    # re-imported by its dotted path when core.py's auto-discovery walks the package, creating
    # two separate module objects (two separate CAPTURED lists / RecordingEpsFBAStep classes),
    # so the composite ends up running the OTHER module's class while this one's CAPTURED stays
    # empty. Calling main() only after a normal dotted import avoids that entirely.
    main()
