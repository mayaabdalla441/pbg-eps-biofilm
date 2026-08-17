"""Depth-resolved (vertical) extension of the validated well-mixed flow-phase
model. Full design rationale:

The well-mixed model cannot reproduce Nona's real t=72-80h EPS-coverage
decline, for a structural reason confirmed twice: a production-side throttle
alone can only slow growth (never reverse it), and shedding can't outcompete
production once nutrients are effectively unconstrained under flow -- both
negative results share the same root cause, everything acting on ONE pool.

This composite tests the fix: separate a flow-exposed SURFACE layer (where
erosion is real) from a protected INTERIOR/basal layer (where production
continues undisturbed), so removal and production stop being forced to net
out to zero on the same number. `run_penetration_depth_check.py` confirmed
this is also physically grounded, not just structurally convenient: at this
model's REALIZED (not enzymatic-ceiling) uptake flux, the reaction-diffusion
penetration depth shrinks from ~970um (t=5h, negligible biomass) to ~35-70um
(t=40-75h, biomass near carrying capacity) -- well inside the real 150um
chamber height, so a real nutrient gradient across depth is expected to
matter exactly in the window that matters for the late decline.

N independent EpsFBAStep instances (ZERO code changes to that class), one per
depth layer, each with its OWN local biomass/eps state -- attached matrix
does not diffuse, only amino-acid CONCENTRATIONS diffuse between layers, via
DepthDiffusionStep (see depth_diffusion_step.py -- a small local
finite-difference process; NOT spatio_flux's DiffusionAdvection, which has a
confirmed shape-order bug for any non-square grid like the thin 1xN one this
model needs). Only the SURFACE layer (index 0) is exposed to flow shear --
only it gets a nonzero eps_shed_rate; interior/basal layers are protected.
"""
from viva_superpowers.composite_generator import composite_generator

from pbg_eps_biofilm.composites.eps_biofilm_flow_scaffold import (
    DEFAULT_AA_BOUNDS, DEFAULT_INITIAL_AA_CONC,
)

NUM_LAYERS_DEFAULT = 3
CHAMBER_HEIGHT_CM = 150e-4  # 150um, Nona's real chamber height
D_AQ_CM2_PER_S = 1e-5  # representative small-solute aqueous diffusivity
D_AQ_CM2_PER_H = D_AQ_CM2_PER_S * 3600.0
STEWART_ORGANIC_RATIO = 0.25  # Stewart 2003, D_e/D_aq for organic solutes -- already an
                              # accepted citation elsewhere in this codebase (see
                              # eps_spatioflux__step.py's diffusion_limitation_floor comment).
                              # run_penetration_depth_check.py confirmed this makes a real,
                              # relevant-scale gradient across the 150um chamber height once
                              # the biofilm matures past ~t=40h -- see module docstring above.
D_EFF_CM2_PER_H_DEFAULT = STEWART_ORGANIC_RATIO * D_AQ_CM2_PER_H

BIOMASS_REACTION_ID = "BIOMASS_Ec_iML1515_core_75p37M"
EPS_REACTION_ID = "EX_eps_e"


def _eps_fba_layer(i, *, dt, model_file, growth_floor_fraction, growth_rate_cap_fraction,
                    layer_biomass_carrying_capacity, eps_shed_rate, aa_bounds):
    return {
        "_type": "process",
        "address": "local:EpsFBAStep",
        "config": {
            "model_file": model_file,
            "biomass_reaction_id": BIOMASS_REACTION_ID,
            "eps_reaction_id": EPS_REACTION_ID,
            "growth_floor_fraction": growth_floor_fraction,
            "eps_shed_rate": eps_shed_rate,
            "growth_rate_cap_fraction": growth_rate_cap_fraction,
            # flow_start_time_h intentionally NOT set -- left at EpsFBAStep's own no-op
            # default (1e9, "never flowing"), so every layer ALWAYS uses the static/
            # depletion-cap accounting branch (real local uptake bounded by what's
            # actually present, real substrate deltas reported). This is deliberate,
            # not an oversight: EpsFBAStep's "flowing" branch hard-zeroes reported
            # substrate deltas and ignores local depletion -- a real, well-mixed-bulk
            # approximation (justified by ~100 chamber-volumes/hr turnover) that is
            # correct for the WELL-MIXED model but wrong per-layer here, since it would
            # suppress the exact local-depletion signal the diffusion step needs to
            # respond to, silently flattening any gradient regardless of how deep or
            # dense the layer is. In THIS model, "flow" is represented structurally
            # instead: the dirichlet boundary condition on DepthDiffusionStep pins the
            # surface layer's neighboring concentration at the real bulk/inlet value,
            # and diffusion (not a flowing-flag) carries that resupply down through the
            # stack -- local consumption at each layer is real and can outpace it,
            # which is the actual diffusion-limitation mechanism being tested. First
            # confirmed necessary via mechanics_run() in run_depth_scaffold.py: with
            # flow_start_time_h set to the real 3.0h (as in the well-mixed model), all
            # 3 layers showed byte-identical biomass and a constant, undepleted AA
            # concentration for the full 80h -- the flowing branch was silently
            # preventing any spatial structure from ever developing.
            "biomass_carrying_capacity": layer_biomass_carrying_capacity,
            "aa_bounds": aa_bounds,
            # diffusion_limitation_scale/floor stay at their no-op defaults -- real
            # depth-resolved diffusion limitation is now handled by DepthDiffusionStep
            # directly; using both would double-count the same physics two ways.
        },
        "inputs": {
            "substrates": {aa_id: ["fields", "substrates", aa_id, str(i)] for aa_id in aa_bounds},
            "biomass": ["layers", str(i), "biomass"],
            "eps": ["layers", str(i), "eps"],
            "global_time": ["global_time"],
        },
        "outputs": {
            "substrates": {aa_id: ["fields", "substrates", aa_id, str(i)] for aa_id in aa_bounds},
            "biomass": ["layers", str(i), "biomass"],
            "eps": ["layers", str(i), "eps"],
        },
        "interval": dt,
    }


@composite_generator(
    name="eps_biofilm_depth_scaffold",
    description=(
        "Depth-resolved (N-layer, vertical through the real 150um chamber height) "
        "extension of the validated well-mixed flow-phase model. Tests whether a "
        "spatially-differentiated surface (erosion, via eps_shed_rate on layer 0 only) "
        "vs. protected interior (undisturbed production) structure can reproduce the real "
        "t=72-80h EPS-coverage decline that the well-mixed model cannot -- see "
        "eps_biofilm_depth_scaffold.py's module docstring and "
        "run_penetration_depth_check.py for the full rationale/verification."
    ),
    parameters={
        "model_file": {"type": "string", "default": "eps_ecoli_model_lb_epspool.xml"},
        "num_layers": {"type": "integer", "default": NUM_LAYERS_DEFAULT},
        "growth_floor_fraction": {"type": "float", "default": 0.9},
        "growth_rate_cap_fraction": {"type": "float", "default": 0.095},
        "initial_biomass": {"type": "float", "default": 2.4e-7},
        "biomass_carrying_capacity": {"type": "float", "default": 4.17e-4},
        "dt": {"type": "float", "default": 0.01},
        "surface_eps_shed_rate": {"type": "float", "default": 0.0},
        "diffusion_coeff_cm2_per_h": {"type": "float", "default": D_EFF_CM2_PER_H_DEFAULT},
    },
)
def eps_biofilm_depth_scaffold(core=None, *, model_file="eps_ecoli_model_lb_epspool.xml",
                                num_layers=NUM_LAYERS_DEFAULT, growth_floor_fraction=0.9,
                                growth_rate_cap_fraction=0.095, initial_biomass=2.4e-7,
                                biomass_carrying_capacity=4.17e-4, dt=0.01,
                                surface_eps_shed_rate=0.0,
                                diffusion_coeff_cm2_per_h=D_EFF_CM2_PER_H_DEFAULT):
    aa_bounds = DEFAULT_AA_BOUNDS
    layer_biomass0 = initial_biomass / num_layers
    layer_bmax = biomass_carrying_capacity / num_layers
    depth_per_layer_cm = CHAMBER_HEIGHT_CM / num_layers

    state = {}
    for i in range(num_layers):
        # Only the surface layer (index 0) is exposed to flow shear.
        shed_rate = surface_eps_shed_rate if i == 0 else 0.0
        state[f"eps_fba[{i}]"] = _eps_fba_layer(
            i, dt=dt, model_file=model_file,
            growth_floor_fraction=growth_floor_fraction,
            growth_rate_cap_fraction=growth_rate_cap_fraction,
            layer_biomass_carrying_capacity=layer_bmax,
            eps_shed_rate=shed_rate,
            aa_bounds=aa_bounds,
        )

    state["diffusion"] = {
        "_type": "process",
        "address": "local:DepthDiffusionStep",
        "config": {
            "num_layers": num_layers,
            "aa_ids": list(aa_bounds.keys()),
            "diffusion_coeff_cm2_per_h": diffusion_coeff_cm2_per_h,
            "depth_per_layer_cm": depth_per_layer_cm,
            # Dirichlet boundary at the surface layer: pinned to the same real bulk/
            # inlet concentration already used to justify "pinned near inlet" during
            # flow in the well-mixed model (fast chamber turnover, ~100 volumes/hr).
            "bulk_concentration": dict(DEFAULT_INITIAL_AA_CONC),
        },
        "inputs": {"substrates": ["fields", "substrates"]},
        "outputs": {"substrates": ["fields", "substrates"]},
        "interval": dt,
    }

    # Every layer starts at the same real bulk/inlet concentration -- no biofilm
    # structure exists yet at t=0 to have created a gradient.
    state["fields"] = {
        "substrates": {
            aa_id: {str(i): conc for i in range(num_layers)}
            for aa_id, conc in DEFAULT_INITIAL_AA_CONC.items()
        },
    }
    state["layers"] = {
        str(i): {"biomass": layer_biomass0, "eps": 0.0}
        for i in range(num_layers)
    }
    return state
