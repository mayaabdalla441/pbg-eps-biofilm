"""Validates the Stage 1 flow-phase MECHANICS (flow_phase_open_issues memory items D.8/D.9)
using the eps_biofilm_flow_scaffold composite -- NOT a calibrated model. Two placeholder
mechanisms under test, using illustrative (not literature-derived) numbers:

  1. Depletion-cap bypass at flow_start_time_h=3.0h (matches Nona's real protocol timing --
     static 0-3h, flow from t=3h -- justified by her real chamber turnover math).
  2. A logistic growth cap at biomass_carrying_capacity=2.0 (a PLACEHOLDER chosen only to sit
     clearly above the static model's natural ~1.63 plateau, so the cap -- not nutrient
     depletion -- is visibly what stops growth once flow removes the depletion limiter).

What this run should show, if the mechanics are correct:
  - t=0-3h: same depletion-limited growth as the validated static model (this composite must
    match run_static_composite.py's trajectory exactly up to t=3h, since flow_start_time_h=3.0
    means nothing new is active yet).
  - t=3h: a visible inflection as the depletion cap switches off and growth accelerates again
    (biomass was well below the nutrient-only plateau by t=3h, so there's room to grow).
  - t=3h onward: growth curves smoothly into the carrying capacity (approaches, does not
    overshoot, biomass_carrying_capacity=2.0) rather than blowing up unboundedly.
  - EPS keeps accumulating even once biomass saturates (unchanged design decision -- EPS is not
    capped by biomass_carrying_capacity).

Guarded by `if __name__ == "__main__":` for the same reason as the other run_*.py scripts.
"""
import warnings
warnings.filterwarnings('ignore')

from process_bigraph import Composite
from process_bigraph.emitter import gather_emitter_results, emitter_from_wires

from pbg_eps_biofilm.core import build_core
from pbg_eps_biofilm.composites.eps_biofilm_flow_scaffold import eps_biofilm_flow_scaffold


def run(total_hours=10.0, **generator_kwargs):
    core = build_core()
    state = eps_biofilm_flow_scaffold(**generator_kwargs)
    state['emitter'] = emitter_from_wires({
        "global_time": ["global_time"],
        "fields": ["fields"],
    })

    sim = Composite({"state": state}, core=core)
    sim.run(total_hours)

    return gather_emitter_results(sim)[('emitter',)]


def _aa_total(substrates):
    return sum(substrates.values())


def plot(results, flow_start_time_h=3.0, biomass_carrying_capacity=2.0, outfile="flow_scaffold_trajectory.png"):
    import matplotlib.pyplot as plt

    t = [r['global_time'] for r in results]
    biomass = [r['fields']['biomass'] for r in results]
    eps = [r['fields']['eps'] for r in results]
    aa_total = [_aa_total(r['fields']['substrates']) for r in results]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    series = [
        (axes[0], t, biomass, "Biomass (gDW/L)", "#2a78d6"),
        (axes[1], t, eps, "Accumulated EPS[e] (mg)", "#eb6834"),
        (axes[2], t, aa_total, "Total amino acids remaining (mM, summed)", "#1baf7a"),
    ]
    for ax, xs, ys, title, color in series:
        ax.plot(xs, ys, linewidth=2, color=color)
        ax.axvline(flow_start_time_h, linestyle='--', linewidth=1, color='#94a3b8')
        ax.set_title(title)
        ax.set_xlabel("t (h)")
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].axhline(biomass_carrying_capacity, linestyle=':', linewidth=1, color='#94a3b8')

    plt.tight_layout()
    plt.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")


if __name__ == "__main__":
    results = run(total_hours=10.0)

    print(f"{'t (h)':>8} {'biomass':>12} {'eps (mg)':>12} {'aa_total_remaining':>20}")
    for r in results[::50]:
        t = r['global_time']
        fields = r['fields']
        print(f"{t:8.2f} {fields['biomass']:12.4f} {fields['eps']:12.4f} {_aa_total(fields['substrates']):20.4f}")

    plot(results)
