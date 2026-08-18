# Legacy GEM analysis

This is the original, notebook-based genome-scale metabolic modeling work for
the E. coli EPS/biofilm project, copied over from the `string-analysis-ecoliGEM`
repo (kept there too, unchanged — this is a one-way copy, not a migration).

It predates and forms the basis for the process-bigraph port in
`pbg_eps_biofilm/` (the active, current model — see that package for the
real `EpsFBAStep` two-stage dFBA implementation, the real-data calibration,
and the ongoing flow-phase/depth-resolved investigation).

## What's here

- **`notebooks/`** — the original COBRApy/Jupyter analysis:
  - `eps_analysis.ipynb` — master six-condition FBA analysis, all 3 EPS pathways (cellulose, colanic acid, PNAG)
  - `eps_media_comparison.ipynb` — M9 vs LB medium comparison
  - `eps_model_KO.ipynb` — 8-gene knockout screen
  - `eps_ecoli_model.ipynb` — the original dFBA prototype (Static Optimization Approach), historical basis for `EpsFBAStep`
  - `eps_competition_analysis.ipynb` — investigation into whether the 3 EPS pathways compete for shared precursors when optimized jointly vs. independently
- **`models/`** — iML1515 with the 3 EPS pathways manually added (cellulose/CELSYNTH, colanic acid/COLASYNTH, PNAG/PUACGAMex), in a few variants used across the notebooks above. Does NOT include the raw unmodified iML1515 base model (publicly downloadable from BiGG, not unique work product).
- **`results/`** — CSV outputs from the notebooks above (six-condition tables, knockout screen results, M9/LB comparisons, etc.)
- **`report/`** — `project_summary.txt` is the authoritative written summary of this original computational work; has real content relevant to the paper's §2.2 "Modeling Strategy" section. `eps_all_pathways_report.html` is a rendered report.
- **`images/`** — plots referenced by `project_summary.txt` and the notebooks (pathway diagrams, condition comparisons, production envelopes).

## Why keep this instead of just the port

The port (`pbg_eps_biofilm/`) reimplements the two-stage FBA logic but
doesn't carry over the original six-condition matrix, M9/LB comparison, gene
knockout screen, or pathway-competition investigation as runnable
notebooks — those results only exist here. Kept for provenance and because
they're directly needed for parts of the paper write-up.
