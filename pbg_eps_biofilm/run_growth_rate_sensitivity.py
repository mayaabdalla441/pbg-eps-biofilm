"""Sensitivity test A.1 (flow_phase_open_issues memory): how much does time-to-depletion and
total EPS shift if cells don't grow at their unconstrained FBA-optimal rate?

Real E. coli biology (RpoS/CsgD cascade) ties curli/cellulose matrix production to early
stationary phase, not to continuous near-maximal growth -- our default growth_floor_fraction=0.9
design assumes EPS production happens alongside near-maximal growth throughout. This script
doesn't resolve that modeling question; it quantifies how sensitive the batch-phase trajectory
is to growth rate, as a first step toward deciding whether it matters.

Sweeps EpsFBAStep's growth_rate_cap_fraction (a no-op at 1.0, reproduces the validated static-
phase result exactly) across [1.0, 0.75, 0.5, 0.25] and plots biomass/eps/aa_total together.

Guarded by `if __name__ == "__main__":` for the same reason as run_static_composite.py --
build_core()'s auto-discovery imports every top-level module in this package while registering
Process/Step classes, and must never trigger a real simulation run as a side effect of that.
"""
import warnings
warnings.filterwarnings('ignore')

from pbg_eps_biofilm.run_static_composite import run

CAP_FRACTIONS = [1.0, 0.75, 0.5, 0.25]
# Fixed categorical color order, consistent with the palette used in run_static_composite.py's
# biomass/eps/aa_total panels (blue/orange/aqua) -- extended with one more categorical hue.
COLORS = {1.0: "#2a78d6", 0.75: "#eb6834", 0.5: "#1baf7a", 0.25: "#8858d0"}


def run_sweep(total_hours=5.0):
    results = {}
    for cap in CAP_FRACTIONS:
        print(f"running growth_rate_cap_fraction={cap} ...")
        results[cap] = run(total_hours=total_hours, growth_rate_cap_fraction=cap)
    return results


def _aa_total(substrates):
    return sum(substrates.values())


def summarize(results):
    print(f"\n{'cap_fraction':>12} {'final_biomass':>14} {'final_eps (mg)':>16} {'t_plateau (h)':>14}")
    for cap, rows in results.items():
        biomass = [r['fields']['biomass'] for r in rows]
        eps = [r['fields']['eps'] for r in rows]
        t = [r['global_time'] for r in rows]
        final_biomass, final_eps = biomass[-1], eps[-1]
        # first time biomass stops changing meaningfully step-to-step (plateau onset)
        t_plateau = next(
            (t[i] for i in range(1, len(biomass)) if abs(biomass[i] - biomass[i - 1]) < 1e-9),
            float('nan'),
        )
        print(f"{cap:12.2f} {final_biomass:14.4f} {final_eps:16.4f} {t_plateau:14.2f}")


def plot(results, outfile="growth_rate_sensitivity.png"):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for cap, rows in results.items():
        t = [r['global_time'] for r in rows]
        biomass = [r['fields']['biomass'] for r in rows]
        eps = [r['fields']['eps'] for r in rows]
        aa_total = [_aa_total(r['fields']['substrates']) for r in rows]
        color = COLORS[cap]
        label = f"{cap:.2f}×"
        axes[0].plot(t, biomass, linewidth=2, color=color, label=label)
        axes[1].plot(t, eps, linewidth=2, color=color, label=label)
        axes[2].plot(t, aa_total, linewidth=2, color=color, label=label)

    titles = ["Biomass (gDW/L)", "Accumulated EPS[e] (mg)", "Total amino acids remaining (mM, summed)"]
    for ax, title in zip(axes, titles):
        ax.set_title(title)
        ax.set_xlabel("t (h)")
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(title="growth rate cap", fontsize=8, title_fontsize=8)

    plt.tight_layout()
    plt.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")


if __name__ == "__main__":
    results = run_sweep(total_hours=5.0)
    summarize(results)
    plot(results)
