"""Driver for eps_biofilm_depth_scaffold, mirroring run_flow_scaffold.py's pattern.
See /Users/mayaabdalla/.claude/plans/pure-gliding-ocean.md for the full plan.

Three phases, run in order:
  1. smoke_test() -- N=2, diffusion off, shedding off: both layers must behave
     identically. Catches wiring/orientation bugs before trusting a long run.
  2. mechanics_run() -- short window, real-calibrated numbers, diffusion ON,
     shedding OFF: confirm the expected qualitative behavior (a real depth
     gradient develops as biomass builds up, per run_penetration_depth_check.py).
  3. calibrated_sweep() -- full 80h, real-calibrated numbers, sweep the
     surface-only eps_shed_rate across the same style of explicit half-life
     assumption used in run_eps_shedding_sensitivity.py (not one guessed
     value), aggregate EPS across layers, compare against real coverage data
     using the ALREADY-FITTED mass-to-coverage transform from
     run_mass_to_coverage_mapping.py (scale=10.93mg, not refit).

Guarded by `if __name__ == "__main__":` for the same reason as the other
run_*.py scripts -- but see run_penetration_depth_check.py's note: invoke via
`python -c "import ...; m.main()"`, not `python -m`, to avoid the __main__/
dotted-import double-module issue when this file also gets walked by
core.py's auto-discovery.
"""
import warnings
warnings.filterwarnings('ignore')

import math

import numpy as np
from process_bigraph import Composite
from process_bigraph.emitter import gather_emitter_results, emitter_from_wires

from pbg_eps_biofilm.core import build_core
from pbg_eps_biofilm.composites.eps_biofilm_depth_scaffold import eps_biofilm_depth_scaffold
from pbg_eps_biofilm.run_mass_to_coverage_mapping import load_real_data

REAL_CALIBRATION = dict(
    initial_biomass=2.4e-7,
    biomass_carrying_capacity=4.17e-4,
    growth_rate_cap_fraction=0.095,
    # flow_start_time_h intentionally omitted -- see eps_biofilm_depth_scaffold.py's
    # _eps_fba_layer() comment: every layer always uses the static/depletion-cap
    # accounting branch, "flow" is represented by the diffusion boundary instead.
)
FITTED_MASS_TO_COVERAGE_SCALE = 10.93  # mg, from run_mass_to_coverage_mapping.py -- reused, not refit

HALF_LIVES_H = [None, 6, 12, 24, 48]  # None = no shedding (control); explicit swept assumption,
                                       # same style as run_eps_shedding_sensitivity.py


def _shed_rate_for(half_life_h):
    return 0.0 if half_life_h is None else math.log(2) / half_life_h


def run(total_hours, dt, **generator_kwargs):
    core = build_core()
    state = eps_biofilm_depth_scaffold(dt=dt, **generator_kwargs)
    state['emitter'] = emitter_from_wires({
        "global_time": ["global_time"],
        "fields": ["fields"],
        "layers": ["layers"],
    })
    sim = Composite({"state": state}, core=core)
    sim.run(total_hours)
    return gather_emitter_results(sim)[('emitter',)]


def smoke_test():
    """N=2, diffusion negligible (tiny D, effectively decoupled), shedding off:
    both layers should be numerically identical throughout."""
    print("=== Smoke test: N=2, symmetric (diffusion~0, shedding off) ===")
    results = run(
        total_hours=2.0, dt=0.01, num_layers=2,
        diffusion_coeff_cm2_per_h=1e-12,  # ~decoupled: no meaningful inter-layer flux
        surface_eps_shed_rate=0.0,
        **REAL_CALIBRATION,
    )
    last = results[-1]['layers']
    b0, b1 = last['0']['biomass'], last['1']['biomass']
    e0, e1 = last['0']['eps'], last['1']['eps']
    rel_diff_b = abs(b0 - b1) / max(b0, 1e-30)
    rel_diff_e = abs(e0 - e1) / max(e0, 1e-30)
    print(f"final biomass: layer0={b0:.6g} layer1={b1:.6g} (rel diff {rel_diff_b:.2e})")
    print(f"final eps:     layer0={e0:.6g} layer1={e1:.6g} (rel diff {rel_diff_e:.2e})")
    ok = rel_diff_b < 1e-6 and rel_diff_e < 1e-6
    print("PASS" if ok else "FAIL -- layers diverged with no real asymmetry active")
    return ok


def mechanics_run(total_hours=80.0, dt=0.05):
    """Real-calibrated numbers, diffusion ON, shedding OFF: confirm a real
    depth gradient develops as biomass builds up (per run_penetration_depth_check.py)."""
    print("\n=== Mechanics run: N=3, real-calibrated, diffusion ON, shedding OFF ===")
    results = run(total_hours=total_hours, dt=dt, num_layers=3,
                   surface_eps_shed_rate=0.0, **REAL_CALIBRATION)
    print(f"{'t (h)':>8} {'B[surf]':>12} {'B[mid]':>12} {'B[basal]':>12} "
          f"{'AA[surf]':>10} {'AA[basal]':>10}")
    for r in results[::max(1, len(results) // 12)]:
        layers = r['layers']
        aa_surf = r['fields']['substrates']['EX_ser__L_e']['0']
        aa_basal = r['fields']['substrates']['EX_ser__L_e']['2']
        print(f"{r['global_time']:8.2f} {layers['0']['biomass']:12.6g} "
              f"{layers['1']['biomass']:12.6g} {layers['2']['biomass']:12.6g} "
              f"{aa_surf:10.4f} {aa_basal:10.4f}")
    return results


def calibrated_sweep(total_hours=80.0, dt=0.01):
    print("\n=== Full real-calibrated 80h run, surface-shedding half-life sweep ===")
    t_real, coverage_real = load_real_data()

    sweep_results = {}
    for hl in HALF_LIVES_H:
        rate = _shed_rate_for(hl)
        label = "no shedding" if hl is None else f"T1/2={hl}h"
        print(f"running surface_eps_shed_rate={rate:.5f} 1/h ({label}) ...")
        results = run(total_hours=total_hours, dt=dt, num_layers=3,
                       surface_eps_shed_rate=rate, **REAL_CALIBRATION)
        t_model = np.array([r['global_time'] for r in results])
        eps_agg = np.array([sum(r['layers'][str(i)]['eps'] for i in range(3)) for r in results])
        sweep_results[hl] = (t_model, eps_agg)

    print(f"\n{'half_life':>10} {'rate (1/h)':>11} {'RMSE (pct-pts)':>15} {'R^2':>8} "
          f"{'decline 70->80h?':>18}")
    best = None
    for hl, (t_model, eps_agg) in sweep_results.items():
        eps_at_real_t = np.interp(t_real, t_model, eps_agg)
        pred = 100.0 * (1.0 - np.exp(-eps_at_real_t / FITTED_MASS_TO_COVERAGE_SCALE))
        resid = pred - coverage_real
        rmse = np.sqrt(np.mean(resid ** 2))
        ss_res = np.sum(resid ** 2)
        ss_tot = np.sum((coverage_real - coverage_real.mean()) ** 2)
        r2 = 1.0 - ss_res / ss_tot
        idx_70 = np.argmin(np.abs(t_real - 70))
        idx_80 = np.argmin(np.abs(t_real - 80))
        declined = pred[idx_80] < pred[idx_70]
        label = "none" if hl is None else f"{hl}h"
        print(f"{label:>10} {_shed_rate_for(hl):11.5f} {rmse:15.3f} {r2:8.4f} "
              f"{'YES' if declined else 'no':>18}")
        if best is None or rmse < best[1]:
            best = (hl, rmse)

    print(f"\nBest fit: half_life={best[0]}, RMSE={best[1]:.3f} pct-points")
    return sweep_results


def main():
    ok = smoke_test()
    if not ok:
        print("\nSmoke test FAILED -- stopping before the longer runs. Check layer "
              "wiring/orientation before trusting mechanics_run/calibrated_sweep.")
        return
    mechanics_run()
    calibrated_sweep()


if __name__ == "__main__":
    main()
