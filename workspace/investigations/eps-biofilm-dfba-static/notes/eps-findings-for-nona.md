# EPS modeling: what we've found (biology summary for Nona)

Short summary of what the genome-scale metabolic model tells us biologically, and where we're stuck — written for discussion, not as paper prose.

## What the model shows

**EPS production is a real growth-tradeoff, not a free byproduct.** The model can redirect substantial metabolic flux into matrix production (cellulose, colanic acid, PNAG) specifically when growth is constrained. Nitrogen limitation is the strongest single trigger — carbon stays abundant, nitrogen scarce — which is the classic biological cue for producing storage/matrix polymers instead of biomass. Cellulose dominates total matrix mass under maximal-EPS conditions.

**Gene-level structure checks out.** A knockout screen confirms wcaJ/cpsG, bcsA, and pgaB each independently and cleanly kill only their own matrix type (colanic acid, cellulose, PNAG respectively) — the three pathways don't cross-talk, matching known biology. The c-di-GMP degradation redundancy (5 isozymes) also came through correctly as a clean negative control.

**Growth and matrix production run on separate clocks.** Fitting the model against your real imaging data, cell growth follows a clean logistic curve — but EPS coverage doesn't track it at all; it has its own multi-phase shape (early dip, flat period, late surge, peak, decline). That's consistent with real biofilm-maturation biology being staged (attachment → microcolony → maturation → dispersion), not just "more cells, more matrix."

**Real growth is ~12x slower than the model's unconstrained optimum.** Calibrating against your Cell occupancy (%) data gives a real growth rate far below what the metabolism alone would allow — consistent with cells inside a growing biofilm being physiologically throttled relative to free-living optimal growth, not a modeling artifact.

**Quantitatively, the model explains most of the real trajectory.** Once accumulated EPS mass (the model's native output) is mapped through a simple saturating transform into a projected-area coverage fraction (your measurement's actual units), it fits your real coverage data at R²=0.90 across the full 8–80h window. That's a real, validated link between predicted metabolic capacity and what you're actually imaging.

## Where we're stuck

The one thing the model consistently cannot reproduce, however we've tried to build it, is the **real decline in EPS coverage at t≈72–80h** (94%→89% in your data). We've now tested this four separate ways — a global production throttle, a first-order removal term on the whole matrix pool, a spatial model with removal isolated to a flow-exposed surface layer, and (after finding and fixing a real unit-conversion bug) that same spatial model re-run correctly — and all four fail for a related reason: production never actually plateaus toward zero in the model (it settles at a reduced but still-positive rate), so no *constant-rate* removal term can ever outpace it.

We also did the physics math on your real chamber directly: at 5µL/min through your geometry, the wall shear stress is only about **0.011 Pa**. A real paper measuring E. coli biofilm mechanical strength found matrix doesn't even start to yield until **~43-92 Pa** (though that was mature, agar-grown biofilm, not your young, liquid-grown one — so this isn't a final answer, but it's a real caution against assuming flow alone is shearing cells off).

We think this isn't a tuning problem — it's a scope boundary. A stoichiometric metabolic model can tell you how much matrix *can* be produced; it structurally cannot represent active dispersal, enzymatic matrix degradation (e.g. BcsZ), or a *triggered* mechanical/biological event, because those aren't metabolic reactions and aren't constant rates. So the late decline is very likely a real, triggered event, not something more flux-balance tuning (or a bigger constant rate) will capture.

## The actual ask

**We need a real flush-out/removal mechanism, and we don't have a citable rate for one.** Four things would each genuinely help, roughly in order of how cheap they are for you to answer:

1. **Shape of the decline, from frames you already have:** looking at t≈65-80h, does the coverage decline look like a smooth, gradual decrease the whole way through, or does it look like it starts more abruptly at some specific point — like a real event happened, not a steady ongoing process? This is the single most useful thing you could tell us — it tells us which whole family of mechanism to pursue.
2. **A visual crowding read, also from frames you already have:** right before the decline starts (~t=65-70h), does the biofilm look completely packed — no visible gaps — compared to earlier timepoints? We have a real hypothesis that physical crowding itself could be triggering an active response (there's a real paper showing self-imposed mechanical stress alone can flip E. coli into a biofilm-associated state), and this would be a cheap first check of whether that's even plausible.
3. **The volumetric/thickness estimate you mentioned being able to produce** (2D coverage × chamber height, or an assumed effective thickness): could you actually generate that trajectory? It would tell us how thick the biofilm gets relative to your 150µm channel by t≈70-80h, which directly determines whether the flow's shear even reaches the biofilm surface at all (real fluid physics: shear is near-zero right at the wall and increases with height, so a thin biofilm may never leave the low-shear zone regardless of flow rate or maturity). This also gives us a real number to check our spatial model against, instead of just coverage %.
4. **Anything you have** — a real matrix turnover/degradation rate (enzymatic, e.g. BcsZ-mediated), the z-stack data (even unanalyzed), or your own hypothesis (a gene, regulatory pathway, or something from other biofilm work) for why this might be an active dispersal response rather than passive erosion.

Any of these would let us build a mechanistically real removal term instead of a swept assumption, which is currently the missing piece.
