"""Re-run of run_depth_scaffold.py's calibrated_sweep() with TWO changes vs the
original (see decline_investigation_summary memory, 2026-08-17/18):

1. Uses the volume-anchoring fix (EpsFBAStep.substrate_reference_volume_L,
   wired to each layer's real physical volume in eps_biofilm_depth_scaffold.py)
   instead of the 1.0L no-op default.
2. REFITS the mass-to-coverage transform's `scale` against the volume-fixed
   no-shedding baseline, instead of reusing the old scale=10.93mg. That old
   scale was calibrated against the PRE-fix model's much larger EPS mass
   range (max ~39mg); the fixed model produces genuinely smaller totals
   (max ~3.5mg in a quick check), so reusing the stale scale gave a
   nonsensical negative R^2 -- a comparison artifact, not a real result.
   Refitting once against the no-shedding baseline (not per-shed-rate, so
   scale can't mask shape differences between runs) restores an honest
   comparison.

Prints incrementally (flushed after each run) so a partial/killed run still
leaves useful output on disk -- this machine has real, unpredictable CPU
contention (other apps/processes), and the previous attempt at this exceeded
the tool's 10-minute wall-clock cap under contention despite the underlying
compute being well within budget when uncontended.

Guarded by `if __name__ == "__main__":` for the same reason as the other
run_*.py scripts.
"""
import warnings
warnings.filterwarnings('ignore')

import sys
import numpy as np

import pbg_eps_biofilm.run_depth_scaffold as rds
from pbg_eps_biofilm.run_mass_to_coverage_mapping import load_real_data, fit_scale


def main(dt=0.1, total_hours=80.0):
    t_real, coverage_real = load_real_data()

    print(f"Running no-shedding baseline (dt={dt}h) to refit the mass-to-coverage scale...", flush=True)
    baseline = rds.run(total_hours=total_hours, dt=dt, num_layers=3,
                        surface_eps_shed_rate=0.0, **rds.REAL_CALIBRATION)
    t_base = np.array([r['global_time'] for r in baseline])
    eps_base = np.array([sum(r['layers'][str(i)]['eps'] for i in range(3)) for r in baseline])
    scale, _, rmse0, r2_0 = fit_scale(t_base, eps_base, t_real, coverage_real)
    print(f"Refit baseline scale={scale:.6g} mg  RMSE={rmse0:.3f} pct-pts  R^2={r2_0:.4f}", flush=True)
    print(f"eps_agg range: {eps_base.min():.6g} -> {eps_base.max():.6g} mg "
          f"(old pre-fix model's range was ~0 -> 39.1 mg, for comparison)", flush=True)

    print(f"\n{'half_life':>10} {'rate (1/h)':>11} {'RMSE':>10} {'R^2':>8} {'decline 70->80h?':>18}", flush=True)
    results_table = []
    for hl in rds.HALF_LIVES_H:
        rate = rds._shed_rate_for(hl)
        label = "no shedding" if hl is None else f"T1/2={hl}h"
        print(f"running surface_eps_shed_rate={rate:.5f} 1/h ({label}) ...", flush=True)
        results = rds.run(total_hours=total_hours, dt=dt, num_layers=3,
                           surface_eps_shed_rate=rate, **rds.REAL_CALIBRATION)
        t = np.array([r['global_time'] for r in results])
        eps_agg = np.array([sum(r['layers'][str(i)]['eps'] for i in range(3)) for r in results])
        eps_at_real_t = np.interp(t_real, t, eps_agg)
        pred = 100.0 * (1.0 - np.exp(-eps_at_real_t / scale))
        resid = pred - coverage_real
        rmse = np.sqrt(np.mean(resid ** 2))
        ss_res = np.sum(resid ** 2)
        ss_tot = np.sum((coverage_real - coverage_real.mean()) ** 2)
        r2 = 1.0 - ss_res / ss_tot
        idx70 = np.argmin(np.abs(t_real - 70))
        idx80 = np.argmin(np.abs(t_real - 80))
        declined = pred[idx80] < pred[idx70]
        row_label = "none" if hl is None else f"{hl}h"
        row = f"{row_label:>10} {rate:11.5f} {rmse:10.3f} {r2:8.4f} {('YES' if declined else 'no'):>18}"
        print(row, flush=True)
        results_table.append((hl, rate, rmse, r2, declined))

    best = min(results_table, key=lambda r: r[2])
    print(f"\nBest fit: half_life={best[0]}, RMSE={best[2]:.3f} pct-points", flush=True)
    any_decline = any(r[4] for r in results_table)
    print(f"Any half-life reproduces the 70-80h decline: {any_decline}", flush=True)
    return scale, results_table


if __name__ == "__main__":
    main()
