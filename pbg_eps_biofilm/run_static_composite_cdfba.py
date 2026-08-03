"""Run the cdFBA-based static-phase EPS biofilm composite and plot the trajectory.

Mirrors run_static_composite.py (the spatio-flux version) exactly in structure --
driver script, not a Process/Step definition and not a composite generator. Guarded
by `if __name__ == "__main__":` for the same reason: build_core()'s auto-discovery
imports every top-level module in this package while registering Process/Step
classes, and must never trigger a real simulation run as a side effect of that.
"""
import warnings
warnings.filterwarnings('ignore')

from process_bigraph import Composite
from process_bigraph.emitter import gather_emitter_results, emitter_from_wires

from pbg_eps_biofilm.core import build_core
from pbg_eps_biofilm.composites.eps_biofilm_static_cdfba import eps_biofilm_static_cdfba, BIOMASS_NAME


def run(total_hours=5.0, **generator_kwargs):
    core = build_core()
    state = eps_biofilm_static_cdfba(**generator_kwargs)
    state['emitter'] = emitter_from_wires({
        "global_time": ["global_time"],
        "shared_environment": ["shared_environment"],
    })

    sim = Composite({"state": state}, core=core)
    sim.run(total_hours)

    return gather_emitter_results(sim)[('emitter',)]


def _aa_total(concentrations):
    return sum(v for k, v in concentrations.items() if k not in (BIOMASS_NAME, "eps"))


def plot(results, outfile="static_composite_trajectory_cdfba.png"):
    import matplotlib.pyplot as plt

    t = [r['global_time'] for r in results]
    biomass = [r['shared_environment']['concentrations'][BIOMASS_NAME] for r in results]
    eps = [r['shared_environment']['concentrations']['eps'] for r in results]
    aa_total = [_aa_total(r['shared_environment']['concentrations']) for r in results]

    # Small multiples, one axis each -- same reasoning as the spatio-flux version:
    # biomass, eps, and aa_total_remaining are on very different scales, so a
    # shared/dual axis would be misleading.
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    series = [
        (axes[0], t, biomass, "Biomass (gDW/L)", "#2a78d6"),
        (axes[1], t, eps, "Accumulated EPS[e] (mg)", "#eb6834"),
        (axes[2], t, aa_total, "Total amino acids remaining (mM, summed)", "#1baf7a"),
    ]
    for ax, xs, ys, title, color in series:
        ax.plot(xs, ys, linewidth=2, color=color)
        ax.set_title(title)
        ax.set_xlabel("t (h)")
        ax.spines[['top', 'right']].set_visible(False)

    plt.tight_layout()
    plt.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")


if __name__ == "__main__":
    results = run(total_hours=5.0)

    print(f"{'t (h)':>8} {'biomass':>12} {'eps (mg)':>12} {'aa_total_remaining':>20}")
    for r in results[::20]:
        t = r['global_time']
        concentrations = r['shared_environment']['concentrations']
        biomass = concentrations[BIOMASS_NAME]
        eps = concentrations['eps']
        aa_total = _aa_total(concentrations)
        print(f"{t:8.2f} {biomass:12.4f} {eps:12.4f} {aa_total:20.4f}")

    plot(results)
