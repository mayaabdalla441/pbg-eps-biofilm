# EPS modeling: what we've found (biology summary for Nona)

Short summary of what the genome-scale metabolic model tells us biologically, and where we're stuck — written for discussion, not as paper prose.

## What the model shows

**EPS production is a real growth-tradeoff, not a free byproduct.** The model can redirect substantial metabolic flux into matrix production (cellulose, colanic acid, PNAG) specifically when growth is constrained. Nitrogen limitation is the strongest single trigger — carbon stays abundant, nitrogen scarce — which is the classic biological cue for producing storage/matrix polymers instead of biomass. Cellulose dominates total matrix mass under maximal-EPS conditions.

**Gene-level structure checks out.** A knockout screen confirms wcaJ/cpsG, bcsA, and pgaB each independently and cleanly kill only their own matrix type (colanic acid, cellulose, PNAG respectively) — the three pathways don't cross-talk, matching known biology. The c-di-GMP degradation redundancy (5 isozymes) also came through correctly as a clean negative control.

**Growth and matrix production run on separate clocks.** Fitting the model against your real imaging data, cell growth follows a clean logistic curve — but EPS coverage doesn't track it at all; it has its own multi-phase shape (early dip, flat period, late surge, peak, decline). That's consistent with real biofilm-maturation biology being staged (attachment → microcolony → maturation → dispersion), not just "more cells, more matrix."

**Real growth is ~12x slower than the model's unconstrained optimum.** Calibrating against your Cell occupancy (%) data gives a real growth rate far below what the metabolism alone would allow — consistent with cells inside a growing biofilm being physiologically throttled relative to free-living optimal growth, not a modeling artifact.

**Quantitatively, the model explains most of the real trajectory.** Once accumulated EPS mass (the model's native output) is mapped through a simple saturating transform into a projected-area coverage fraction (your measurement's actual units), it fits your real coverage data at R²=0.90 across the full 8–80h window. That's a real, validated link between predicted metabolic capacity and what you're actually imaging.

## Where we're stuck

The one thing the model consistently cannot reproduce, however we've tried to build it, is the **real decline in EPS coverage at t≈72–80h** (94%→89% in your data). We've tested this three separate ways — a global production throttle, a first-order removal term on the whole matrix pool, and a spatial model with the removal term isolated to just a flow-exposed surface layer — and all three fail for a related reason: production never actually plateaus in the model, so no constant-rate removal term can ever outpace it.

We think this isn't a tuning problem — it's a scope boundary. A stoichiometric metabolic model can tell you how much matrix *can* be produced; it structurally cannot represent active dispersal, enzymatic matrix degradation (e.g. BcsZ), or shear-driven erosion, because those aren't metabolic reactions. So the late decline is very likely real dispersal-stage biology, not something more flux-balance tuning will capture.

## The actual ask

**We need a real flush-out/removal mechanism, and we don't have a citable rate for one.** Concretely, this could be:
- A real matrix turnover/degradation rate (enzymatic, e.g. BcsZ-mediated) if one exists for MG1655 or a close relative.
- A shear-erosion rate or qualitative onset (does erosion visibly track the continuous flow rate, or does it look more like an active, triggered dispersal event around t≈70h?).
- Anything from the z-stack data (even unanalyzed/preliminary) that shows the matrix thinning or detaching at the surface specifically, rather than uniformly — that would tell us whether "surface-only erosion" is the right structural picture or not.

Any of these would let us build a mechanistically real removal term instead of a swept assumption, which is currently the missing piece.
