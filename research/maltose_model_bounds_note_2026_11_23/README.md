# Reproducible supplement: parameter-bound-dependent model selection for maltose formation

Stefan Traublinger | 23 November 2026 | Computational follow-up to Getreideamylasen und Staerkehydrolyse (2026, DOI https://doi.org/10.5281/zenodo.24378474).

## Purpose and provenance

This folder supports the draft article Parameter-Bound-Dependent Model Selection for Maltose Formation in Wheat Flour. It reanalyses two previously published, technically digitised time courses from Hackenberg et al. (2018), Carbohydrate Polymers 185, 153–158, https://doi.org/10.1016/j.carbpol.2017.12.064. Both curves derive from one study family, not independent studies. Original initial means and digitised figure values were already included in the public monograph and baseline repository.

Input file (already present in this repository): ../zertifizierte_fehlergrenzen/data/original/Hackenberg.csv
Input SHA-256: d40867e07ec1333931a46df1cbe88a0eef42833782d1714f9dab0f9e2d0f7e5e
The original article and technical readings retain their respective source rights; they are not relicensed here.

## Models and evaluation

E(t) = P0 + delta [1 − exp(−k t)]
W(t) = P0 + delta [1 − exp(−(k t)^n)]

The published t=0 mean P0 is fixed. Four early values (5, 10, 20, 60 min) determine the unweighted least-squares fits. The late values (120 and 180 min) are excluded from fitting and used only for retrospective continuation-error assessment. Delta is profiled analytically. Multistart L-BFGS-B bounded optimisation is used. A zero of late RMSE(W) − RMSE(E) defines a conditional numeric crossover, NOT a physical material threshold. Global optimality is not mathematically certified.

Four conditions are evaluated: nominal original values, all readings rounded to 0.01 mg/g, and two deterministic scenarios shifting both late readings together by +0.236 and −0.236 mg/g. These are not confidence intervals, independent laboratory measurements or probability estimates.

## Re-run

From the repository root, install requirements with:
python -m pip install -r research/maltose_model_bounds_note_2026_11_23/requirements.txt

Run:
python research/maltose_model_bounds_note_2026_11_23/reproduce.py

The script generates a CSV, JSON and three SVG figures inside a generated subfolder. A full run can take several minutes. A plot-only run is possible with: python research/maltose_model_bounds_note_2026_11_23/reproduce.py --figures-only

The checked results snapshot is in crossover_scenarios.csv. The documented environment used Python 3.13.5, NumPy 2.3.5 and SciPy 1.17.0. Last decimals can vary slightly across numerical environments.

## Limitations and disclosure

No new experiments or baking trials. The two curves originate from one source family, and the later assessment comprises only two data points per curve. The earlier monograph already demonstrated rank reversals at selected discrete parameter limits. This supplement locates the nominal two-crossover window and examines illustrative sensitivity of its position. Neither model has been independently validated for new batches; no broad biological or technological model-superiority claim is made. The editor must judge the novelty and acceptability of overlap with the pre-existing monograph and separately submitted review.
