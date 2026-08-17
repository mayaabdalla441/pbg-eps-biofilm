"""Mass-to-coverage mapping test (flow_phase_open_issues follow-up, 2026-08-12).

Open question: the model's EPS output is an unbounded accumulated MASS (mg); Nona's
real data is a bounded/saturating projected AREA FRACTION (EPS coverage %, 2D image-
plane measurement). These behave very differently even for identical underlying
biology, so before concluding the model needs a new mechanism to explain the real
72-80h coverage decline, test whether a simple saturating transform of the EXISTING
mass trajectory can reproduce the real curve's shape -- including the decline --
using ONLY the already-validated, real-anchored calibration (initial_biomass=2.4e-7
gDW, biomass_carrying_capacity=4.17e-4 gDW, growth_rate_cap_fraction=0.095,
flow_start_time_h=3.0h -- see flow_phase_open_issues memory item H). No new
simulation mechanism, no new fitted growth/timing parameters -- only ONE new free
parameter (the mapping's `scale`) is fit here.

Mapping tested: coverage_pred(t) = 100 * (1 - exp(-eps_mass(t) / scale))
  - Monotonically increasing in mass (matches "more matrix -> more coverage").
  - Saturates toward 100% as mass grows (matches the real data's ceiling near 94-95%).
  - Critically: since eps_mass(t) is itself NOT monotonic under shedding, but even
    under the current no-shedding (eps_shed_rate=0.0) baseline eps_mass is strictly
    non-decreasing (production only, no removal term) -- so IF eps_mass is strictly
    increasing throughout, this monotonic transform CANNOT produce a coverage
    decline. That's a real, checkable prediction of the hypothesis, not assumed --
    see the printed diagnostic below.

Guarded by `if __name__ == "__main__":` for the same reason as the other run_*.py
scripts.
"""
import warnings
warnings.filterwarnings('ignore')

import csv
import numpy as np

from pbg_eps_biofilm.run_flow_scaffold import run

REAL_DATA_CSV = (
    "workspace/investigations/eps-biofilm-dfba-static/inputs/datasets/"
    "combined_eps_dataset_8_to_80h.csv"
)

# Real-anchored calibration, see flow_phase_open_issues memory item H.
CALIBRATED_PARAMS = dict(
    initial_biomass=2.4e-7,
    biomass_carrying_capacity=4.17e-4,
    growth_rate_cap_fraction=0.095,
    flow_start_time_h=3.0,
    dt=0.01,
)


def load_real_data(path=REAL_DATA_CSV):
    t_real, coverage_real = [], []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            t_real.append(float(row["Experimental time (h)"]))
            coverage_real.append(float(row["EPS coverage (%)"]))
    order = np.argsort(t_real)
    return np.array(t_real)[order], np.array(coverage_real)[order]


def run_calibrated(total_hours=80.0):
    results = run(total_hours=total_hours, **CALIBRATED_PARAMS)
    t = np.array([r["global_time"] for r in results])
    eps_mass = np.array([r["fields"]["eps"] for r in results])
    biomass = np.array([r["fields"]["biomass"] for r in results])
    return t, eps_mass, biomass


def fit_scale(t_model, eps_mass, t_real, coverage_real):
    """Grid + local refine search over `scale` minimizing SSE (no scipy dependency)."""
    eps_at_real_t = np.interp(t_real, t_model, eps_mass)

    def sse(scale):
        pred = 100.0 * (1.0 - np.exp(-eps_at_real_t / scale))
        return np.sum((pred - coverage_real) ** 2)

    # Coarse log-spaced grid across the mass trajectory's own range, then refine.
    lo, hi = eps_at_real_t[eps_at_real_t > 0].min() * 1e-3, eps_at_real_t.max() * 10
    grid = np.geomspace(lo, hi, 400)
    losses = np.array([sse(s) for s in grid])
    best = grid[np.argmin(losses)]

    # Local refine around the coarse best.
    fine = np.linspace(best * 0.5, best * 1.5, 400)
    fine = fine[fine > 0]
    losses_fine = np.array([sse(s) for s in fine])
    best_fine = fine[np.argmin(losses_fine)]

    pred = 100.0 * (1.0 - np.exp(-eps_at_real_t / best_fine))
    resid = pred - coverage_real
    rmse = np.sqrt(np.mean(resid ** 2))
    ss_res = np.sum(resid ** 2)
    ss_tot = np.sum((coverage_real - coverage_real.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot
    return best_fine, pred, rmse, r2


def plot(t_model, eps_mass, t_real, coverage_real, pred, scale, outfile):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    axes[0].plot(t_model, eps_mass, linewidth=2, color="#2a78d6")
    axes[0].set_title("Model: accumulated EPS mass (mg)")
    axes[0].set_xlabel("t (h)")
    axes[0].spines[["top", "right"]].set_visible(False)

    axes[1].scatter(t_real, coverage_real, s=14, color="#94a3b8", label="Real EPS coverage (%)", zorder=2)
    order = np.argsort(t_real)
    axes[1].plot(t_real[order], pred[order], linewidth=2, color="#eb6834",
                 label=f"Predicted: 100*(1-exp(-mass/{scale:.3g}))", zorder=3)
    axes[1].set_title("Mass-to-coverage mapping vs. real data")
    axes[1].set_xlabel("t (h)")
    axes[1].set_ylabel("EPS coverage (%)")
    axes[1].spines[["top", "right"]].set_visible(False)
    axes[1].legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(outfile, dpi=200, bbox_inches="tight")
    print(f"Saved: {outfile}")


if __name__ == "__main__":
    print("Running calibrated flow-scaffold composite for 80h (real B0/B_max)...")
    t_model, eps_mass, biomass = run_calibrated(total_hours=80.0)

    # Diagnostic: is eps_mass itself monotonic? A monotonic transform of a
    # monotonic series cannot produce the real late decline -- check this directly
    # before over-interpreting the fit.
    is_monotonic = np.all(np.diff(eps_mass) >= -1e-12)
    print(f"\neps_mass(t) strictly non-decreasing throughout 0-80h: {is_monotonic}")
    print(f"eps_mass range: {eps_mass.min():.6g} mg -> {eps_mass.max():.6g} mg")
    idx_72 = np.argmin(np.abs(t_model - 72))
    idx_80 = np.argmin(np.abs(t_model - 80))
    print(f"eps_mass @ t=72h: {eps_mass[idx_72]:.6g} mg, @ t=80h: {eps_mass[idx_80]:.6g} mg "
          f"(delta={eps_mass[idx_80]-eps_mass[idx_72]:+.6g})")

    t_real, coverage_real = load_real_data()
    print(f"\nLoaded {len(t_real)} real datapoints, t=[{t_real.min():.1f}, {t_real.max():.1f}]h, "
          f"coverage=[{coverage_real.min():.1f}, {coverage_real.max():.1f}]%")

    scale, pred, rmse, r2 = fit_scale(t_model, eps_mass, t_real, coverage_real)
    print(f"\nFitted scale = {scale:.6g} mg   RMSE = {rmse:.3f} pct-points   R^2 = {r2:.4f}")

    # Late-decline check: real data declines 72->80h (94.4% -> 89.0%, see
    # flow_phase_open_issues memory item C.7's table). Does the fitted mapping?
    mask_late = t_real >= 70
    print("\nLate-window (t>=70h) comparison:")
    print(f"{'t (h)':>8} {'real (%)':>10} {'predicted (%)':>14}")
    for tt, real_c, pred_c in zip(t_real[mask_late], coverage_real[mask_late], pred[mask_late]):
        print(f"{tt:8.2f} {real_c:10.2f} {pred_c:14.2f}")

    plot(t_model, eps_mass, t_real, coverage_real, pred, scale,
         outfile="workspace/investigations/eps-biofilm-dfba-static/inputs/mass_to_coverage_mapping.png")
