"""Design question C.7 (flow_phase_open_issues memory): EPS currently accumulates without any
ceiling once biomass plateaus -- fine over the validated 5h window, not fine over a 72h flow-
phase horizon. Adds a first-order shedding term (d(EPS)/dt = production - eps_shed_rate * EPS)
instead of a hard cap, since eps_shed_rate is a RATE (1/h), not an absolute mass -- it doesn't
need the same real-chamber-scale unit anchoring an absolute EPS_max would (see the memory entry
for why an absolute cap isn't buildable yet: a literature EPS-areal-mass number, applied to
Nona's chamber area, is ~4 orders of magnitude off from this model's current output, the same
unit-anchoring gap already blocking initial_biomass/B_max).

No citable literature rate constant was found for E. coli matrix turnover or shear-driven
erosion despite a real search -- eps_shed_rate is an explicit, flagged ASSUMPTION, like
growth_floor_fraction. Rather than picking one value, this sweeps turnover half-lives spanning
"cap bites hard" to "barely matters within the window": no shedding (control), and half-lives of
12h, 24h, 48h, 96h (half-life T -> rate = ln(2)/T).

Run over 15h (not just the validated 5h) so BOTH phases are visible: the active-growth-phase
effect on accumulation rate, and the post-plateau exponential decay once production stops
(growth and EPS production both go to zero once the binding amino acids are fully depleted --
confirmed by the flat post-plateau tail in the un-shed baseline -- so after ~4.4h this becomes a
clean test of the decay term in isolation, `d(EPS)/dt = -eps_shed_rate * EPS`).

Guarded by `if __name__ == "__main__":` for the same reason as the other run_*.py scripts.
"""
import math
import warnings
warnings.filterwarnings('ignore')

from pbg_eps_biofilm.run_static_composite import run

HALF_LIVES_H = [None, 12, 24, 48, 96]  # None = no shedding (control)
COLORS = {None: "#2a78d6", 12: "#eb6834", 24: "#1baf7a", 48: "#8858d0", 96: "#c0392b"}


def _rate_for(half_life_h):
    return 0.0 if half_life_h is None else math.log(2) / half_life_h


def run_sweep(total_hours=15.0):
    results = {}
    for hl in HALF_LIVES_H:
        rate = _rate_for(hl)
        label = "no shedding" if hl is None else f"T1/2={hl}h"
        print(f"running eps_shed_rate={rate:.5f} 1/h ({label}) ...")
        results[hl] = run(total_hours=total_hours, eps_shed_rate=rate)
    return results


def summarize(results):
    print(f"\n{'half_life':>12} {'rate (1/h)':>12} {'final_eps (mg)':>16} {'eps @ t=4.4h':>14}")
    for hl, rows in results.items():
        eps = [r['fields']['eps'] for r in rows]
        t = [r['global_time'] for r in rows]
        final_eps = eps[-1]
        eps_at_plateau = min(rows, key=lambda r: abs(r['global_time'] - 4.4))['fields']['eps']
        label = "none" if hl is None else f"{hl}h"
        print(f"{label:>12} {_rate_for(hl):12.5f} {final_eps:16.4f} {eps_at_plateau:14.4f}")


def plot(results, outfile="eps_shedding_sensitivity.png"):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for hl, rows in results.items():
        t = [r['global_time'] for r in rows]
        biomass = [r['fields']['biomass'] for r in rows]
        eps = [r['fields']['eps'] for r in rows]
        color = COLORS[hl]
        label = "no shedding" if hl is None else f"T1/2={hl}h"
        axes[0].plot(t, biomass, linewidth=2, color=color, label=label)
        axes[1].plot(t, eps, linewidth=2, color=color, label=label)

    for ax, title in zip(axes, ["Biomass (gDW/L)", "Accumulated EPS[e] (mg)"]):
        ax.set_title(title)
        ax.set_xlabel("t (h)")
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")


if __name__ == "__main__":
    results = run_sweep(total_hours=15.0)
    summarize(results)
    plot(results)
