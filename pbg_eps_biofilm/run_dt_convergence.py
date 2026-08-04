"""Convergence test A.3 (flow_phase_open_issues memory): is dt=0.01h numerically converged, or
is the ~4.4h plateau timing (and the depletion-cascade shape leading up to it) sensitive to
timestep coarseness? The pre-port (non-process-bigraph) version of this project found a real
~9% timing difference between dt=0.1h and dt=0.001h -- this was never re-checked after the port,
and never checked against dt=0.01h specifically (the value actually used everywhere in this
repo).

Scope note: this compares dt=0.01h against a finer dt=0.001h over the SAME window the static
model already covers (0-5h, through the full depletion cascade and plateau) -- NOT a genuine
72h check. Nothing changes in this model after the plateau (a closed batch system with no flow
yet), so a flat, unchanging tail is trivially "converged" regardless of dt and wouldn't add any
real information; the numerically interesting part is entirely within the already-covered
window.

Guarded by `if __name__ == "__main__":` for the same reason as the other run_*.py scripts.
"""
import warnings
warnings.filterwarnings('ignore')

from pbg_eps_biofilm.run_static_composite import run

DT_VALUES = [0.01, 0.001]
COLORS = {0.01: "#2a78d6", 0.001: "#eb6834"}


def run_sweep(total_hours=5.0):
    results = {}
    for dt in DT_VALUES:
        print(f"running dt={dt} ({int(total_hours / dt)} steps) ...")
        results[dt] = run(total_hours=total_hours, dt=dt)
    return results


def _nearest_row(rows, t_target):
    return min(rows, key=lambda r: abs(r['global_time'] - t_target))


def summarize(results):
    coarse, fine = results[0.01], results[0.001]
    coarse_final, fine_final = coarse[-1], fine[-1]
    print(f"\n{'':>8} {'biomass':>12} {'eps (mg)':>12}")
    print(f"{'dt=0.01':>8} {coarse_final['fields']['biomass']:12.4f} {coarse_final['fields']['eps']:12.4f}")
    print(f"{'dt=0.001':>8} {fine_final['fields']['biomass']:12.4f} {fine_final['fields']['eps']:12.4f}")
    biomass_pct = 100 * abs(coarse_final['fields']['biomass'] - fine_final['fields']['biomass']) / fine_final['fields']['biomass']
    eps_pct = 100 * abs(coarse_final['fields']['eps'] - fine_final['fields']['eps']) / fine_final['fields']['eps']
    print(f"\nfinal biomass relative diff: {biomass_pct:.2f}%")
    print(f"final eps relative diff:     {eps_pct:.2f}%")

    # plateau onset: first point where biomass stops changing step-to-step, per dt's own resolution
    print(f"\n{'':>8} {'t_plateau (h)':>14}")
    for dt, rows in results.items():
        biomass = [r['fields']['biomass'] for r in rows]
        t = [r['global_time'] for r in rows]
        tol = 1e-6 if dt == 0.01 else 1e-8  # tighter tol at finer dt, same relative step-size logic
        t_plateau = next(
            (t[i] for i in range(1, len(biomass)) if abs(biomass[i] - biomass[i - 1]) < tol),
            float('nan'),
        )
        print(f"{f'dt={dt}':>8} {t_plateau:14.3f}")

    # max relative deviation of the coarse trajectory from the fine one, at matching timepoints
    max_dev = 0.0
    max_dev_t = None
    for r in coarse:
        t = r['global_time']
        ref = _nearest_row(fine, t)
        if ref['fields']['biomass'] > 1e-9:
            dev = abs(r['fields']['biomass'] - ref['fields']['biomass']) / ref['fields']['biomass']
            if dev > max_dev:
                max_dev, max_dev_t = dev, t
    print(f"\nmax relative biomass deviation (dt=0.01 vs dt=0.001, matched by nearest time): "
          f"{100 * max_dev:.2f}% at t={max_dev_t:.2f}h")


def plot(results, outfile="dt_convergence.png"):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for dt, rows in results.items():
        t = [r['global_time'] for r in rows]
        biomass = [r['fields']['biomass'] for r in rows]
        eps = [r['fields']['eps'] for r in rows]
        color = COLORS[dt]
        axes[0].plot(t, biomass, linewidth=2, color=color, label=f"dt={dt}")
        axes[1].plot(t, eps, linewidth=2, color=color, label=f"dt={dt}")

    for ax, title in zip(axes, ["Biomass (gDW/L)", "Accumulated EPS[e] (mg)"]):
        ax.set_title(title)
        ax.set_xlabel("t (h)")
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")


if __name__ == "__main__":
    results = run_sweep(total_hours=5.0)
    summarize(results)
    plot(results)
