from viva_superpowers.composite_generator import composite_generator

DEFAULT_AA_BOUNDS = {
    'EX_ala__L_e': 47.57, 'EX_arg__L_e': 21.36, 'EX_asn__L_e': 31.33, 'EX_asp__L_e': 23.01,
    'EX_cys__L_e': 5.62,  'EX_gln__L_e': 68.50, 'EX_glu__L_e': 86.47, 'EX_gly_e': 24.40,
    'EX_his__L_e': 18.69, 'EX_ile__L_e': 49.04, 'EX_leu__L_e': 90.65, 'EX_lys__L_e': 58.95,
    'EX_met__L_e': 24.74, 'EX_phe__L_e': 34.52, 'EX_pro__L_e': 96.80, 'EX_ser__L_e': 65.67,
    'EX_thr__L_e': 40.47, 'EX_trp__L_e': 6.21,  'EX_tyr__L_e': 32.93, 'EX_val__L_e': 66.67,
}
DEFAULT_INITIAL_AA_CONC = {
    'EX_ala__L_e': 4.418, 'EX_arg__L_e': 1.984, 'EX_asn__L_e': 2.910, 'EX_asp__L_e': 2.137,
    'EX_cys__L_e': 0.522, 'EX_gln__L_e': 6.362, 'EX_glu__L_e': 8.032, 'EX_gly_e': 2.266,
    'EX_his__L_e': 1.736, 'EX_ile__L_e': 4.555, 'EX_leu__L_e': 8.420, 'EX_lys__L_e': 5.475,
    'EX_met__L_e': 2.298, 'EX_phe__L_e': 3.206, 'EX_pro__L_e': 8.991, 'EX_ser__L_e': 6.100,
    'EX_thr__L_e': 3.759, 'EX_trp__L_e': 0.577, 'EX_tyr__L_e': 3.059, 'EX_val__L_e': 6.192,
}

BIOMASS_NAME = "eps_biofilm"


@composite_generator(
    name="eps_biofilm_static_cdfba",
    description=(
        "Static-phase (no-flow) dFBA for E. coli producing pooled EPS[e] on LB medium, "
        "built on cdFBA's dFBA process (shared_environment/Volumetric convention) instead "
        "of spatio-flux's DynamicFBA. Same two-stage growth-then-EPS[e] logic and hard-cap "
        "amino-acid depletion as the spatio-flux version -- for direct comparison between "
        "the two framework implementations of the same science."
    ),
    parameters={
        "model_file": {"type": "string", "default": "eps_ecoli_model_lb_epspool.xml"},
        "growth_floor_fraction": {"type": "float", "default": 0.9},
        "initial_biomass": {"type": "float", "default": 0.01},
        "dt": {"type": "float", "default": 0.01},
        # First-order EPS turnover/shedding rate (1/h), kept in sync with the spatio-flux
        # composite -- see EpsCdFBAStep's comment. 0.0 is a no-op.
        "eps_shed_rate": {"type": "float", "default": 0.0},
    },
)
def eps_biofilm_static_cdfba(core=None, *, model_file="eps_ecoli_model_lb_epspool.xml",
                              growth_floor_fraction=0.9, initial_biomass=0.01, dt=0.01,
                              eps_shed_rate=0.0):
    # dt=0.01h -- convergence-tested on the spatio-flux implementation 2026-08-04
    # (run_dt_convergence.py); not independently re-checked on this cdFBA-based loop, but both
    # implementations have matched exactly at every other checkpoint so far.
    initial_concentrations = dict(DEFAULT_INITIAL_AA_CONC)
    initial_concentrations[BIOMASS_NAME] = initial_biomass
    initial_concentrations["eps"] = 0.0

    return {
        "eps_fba": {
            "_type": "process",
            "address": "local:EpsCdFBAStep",
            "config": {
                "model_file": model_file,
                "name": BIOMASS_NAME,
                "eps_reaction_id": "EX_eps_e",
                "growth_floor_fraction": growth_floor_fraction,
                "eps_shed_rate": eps_shed_rate,
                "aa_bounds": DEFAULT_AA_BOUNDS,
            },
            "inputs": {
                "shared_environment": ["shared_environment"],
            },
            "outputs": {
                "shared_environment": ["shared_environment"],
            },
            "interval": dt,
        },
        "shared_environment": {
            "concentrations": initial_concentrations,
            # MUST be pre-populated with the same keys as concentrations, not left as
            # {} -- confirmed by testing that Volumetric's counts-delta apply silently
            # drops any key not already present in the current "counts" dict, rather
            # than treating a brand-new key as starting from zero. An empty starting
            # dict here causes the entire shared_environment to go empty after the
            # very first step.
            "counts": dict(initial_concentrations),
            # Placeholder volume (L) -- since our own update() reads "concentrations" and
            # writes "counts" (matching cdFBA's own Volumetric apply-convention), volume=1.0
            # keeps counts and concentrations numerically identical for this first test.
            # A real chamber volume (e.g. Nona's 3uL) would plug in here later.
            "volume": 1.0,
        },
    }
