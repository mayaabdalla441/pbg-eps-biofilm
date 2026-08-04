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


@composite_generator(
    name="eps_biofilm_static",
    description=(
        "Static-phase (no-flow) dFBA for E. coli producing pooled EPS[e] on LB medium. "
        "Two-stage optimization per step: maximize growth, floor at growth_floor_fraction, "
        "then maximize EX_eps_e. Amino-acid uptake hard-capped by remaining concentration "
        "each step (not Michaelis-Menten kinetics)."
    ),
    parameters={
        "model_file": {"type": "string", "default": "eps_ecoli_model_lb_epspool.xml"},
        "growth_floor_fraction": {"type": "float", "default": 0.9},
        "initial_biomass": {"type": "float", "default": 0.01},
        "dt": {"type": "float", "default": 0.01},
        # Sensitivity-test knob (see EpsFBAStep) -- 1.0 is a no-op, reproduces validated behavior.
        "growth_rate_cap_fraction": {"type": "float", "default": 1.0},
    },
)
def eps_biofilm_static(core=None, *, model_file="eps_ecoli_model_lb_epspool.xml",
                        growth_floor_fraction=0.9, initial_biomass=0.01, dt=0.01,
                        growth_rate_cap_fraction=1.0):
    # dt=0.01h is a working assumption, NOT yet convergence-tested on the full loop --
    # an earlier single-step spot-check (dt=0.1h vs dt=0.001h) showed ~9% timing drift
    # in when the amino-acid depletion cascade completes.
    return {
        "eps_fba": {
            "_type": "process",
            "address": "local:EpsFBAStep",
            "config": {
                "model_file": model_file,
                "biomass_reaction_id": "BIOMASS_Ec_iML1515_core_75p37M",
                "eps_reaction_id": "EX_eps_e",
                "growth_floor_fraction": growth_floor_fraction,
                "growth_rate_cap_fraction": growth_rate_cap_fraction,
                "aa_bounds": DEFAULT_AA_BOUNDS,
            },
            "inputs": {
                "substrates": ["fields", "substrates"],
                "biomass": ["fields", "biomass"],
            },
            "outputs": {
                "substrates": ["fields", "substrates"],
                "biomass": ["fields", "biomass"],
                "eps": ["fields", "eps"],
            },
            "interval": dt,
        },
        "fields": {
            "substrates": dict(DEFAULT_INITIAL_AA_CONC),
            "biomass": initial_biomass,
            "eps": 0.0,
        },
    }
