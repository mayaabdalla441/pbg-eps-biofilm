"""Sensitivity test A.2 (flow_phase_open_issues memory): does the growth_floor_fraction=0.9
choice actually matter, or is the model insensitive to it? Directly informs the still-open
investigation.yaml decisions_needed item ("is the 90% growth-floor assumption... biologically
defensible?").

Sweeps growth_floor_fraction over [0.9, 0.7, 0.5] via the eps_biofilm_static composite (already
exposes this as a parameter -- no process/composite code changes needed for this test, unlike
A.1's growth_rate_cap_fraction which required a new knob).

Guarded by `if __name__ == "__main__":` for the same reason as the other run_*.py scripts --
build_core()'s auto-discovery imports every top-level module in this package while registering
Process/Step classes, and must never trigger a real simulation run as a side effect of that.
"""
import warnings
warnings.filterwarnings('ignore')

from pbg_eps_biofilm.run_static_composite import run

FLOOR_FRACTIONS = [0.9, 0.7, 0.5]
# Fixed categorical color order, consistent with the other sensitivity script.
COLORS = {0.9: "#2a78d6", 0.7: "#eb6834", 0.5: "#1baf7a"}


def run_sweep(total_hours=5.0):
    results = {}
    for floor in FLOOR_FRACTIONS:
        print(f"running growth_floor_fraction={floor} ...")
        results[floor] = run(total_hours=total_hours, growth_floor_fraction=floor)
    return results


def _aa_total(substrates):
    return sum(substrates.values())


def summarize(results):
    print(f"\n{'floor_fraction':>14} {'final_biomass':>14} {'final_eps (mg)':>16} {'t_plateau (h)':>14}")
    for floor, rows in results.items():
        biomass = [r['fields']['biomass'] for r in rows]
        eps = [r['fields']['eps'] for r in rows]
        t = [r['global_time'] for r in rows]
        final_biomass, final_eps = biomass[-1], eps[-1]
        t_plateau = next(
            (t[i] for i in range(1, len(biomass)) if abs(biomass[i] - biomass[i - 1]) < 1e-9),
            float('nan'),
        )
        print(f"{floor:14.2f} {final_biomass:14.4f} {final_eps:16.4f} {t_plateau:14.2f}")


def plot(results, outfile="growth_floor_sensitivity.png"):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for floor, rows in results.items():
        t = [r['global_time'] for r in rows]
        biomass = [r['fields']['biomass'] for r in rows]
        eps = [r['fields']['eps'] for r in rows]
        aa_total = [_aa_total(r['fields']['substrates']) for r in rows]
        color = COLORS[floor]
        label = f"{floor:.2f}"
        axes[0].plot(t, biomass, linewidth=2, color=color, label=label)
        axes[1].plot(t, eps, linewidth=2, color=color, label=label)
        axes[2].plot(t, aa_total, linewidth=2, color=color, label=label)

    titles = ["Biomass (gDW/L)", "Accumulated EPS[e] (mg)", "Total amino acids remaining (mM, summed)"]
    for ax, title in zip(axes, titles):
        ax.set_title(title)
        ax.set_xlabel("t (h)")
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(title="growth_floor_fraction", fontsize=8, title_fontsize=8)

    plt.tight_layout()
    plt.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")


if __name__ == "__main__":
    results = run_sweep(total_hours=5.0)
    summarize(results)
    plot(results)
