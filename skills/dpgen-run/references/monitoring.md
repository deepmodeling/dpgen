# DP-GEN Monitoring and Troubleshooting

Load this reference after each iteration or whenever training accuracy stalls.
Here "accuracy" means energy/force train and validation RMSE; model deviation
measures sampling uncertainty, not accuracy.

## Record every iteration

- train/validation energy and force RMSE from `lcurve.out` or equivalent logs;
- model-deviation distributions and the low/high trust-window counts;
- candidate, FP-success, FP-failure, and newly labeled structure counts;
- current `iter.*`, `record.dpgen`, pending jobs, and whether the next iteration
  actually contains the new data.

## Diagnose in this order

| Symptom                                              | First checks                                                                   | Focused adjustment                                                                                    |
| ---------------------------------------------------- | ------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------- |
| Train and validation RMSE both stay high             | Units, labels, `type_map`, data coverage, descriptor capacity                  | Fix data/schema first; then review descriptor, fitting net, learning rate, and steps                  |
| Train RMSE falls but validation RMSE stalls or rises | Split leakage, duplicates, outliers, validation coverage                       | Deduplicate/re-split; add representative labels; reduce overfitting only after data is sound          |
| RMSE is good but model deviation is high             | Unseen compositions/geometries, temperature/pressure range, initial structures | Broaden `model_devi_jobs`, temperatures, ensembles, or initial structures; label the candidate window |
| Candidate count is zero                              | Trust thresholds, MD stability, trajectory frequency, model diversity          | Check thresholds and sampling settings before changing the model                                      |
| Candidate count is excessive or mostly failed FP     | Thresholds too loose, unstable MD, bad structures, FP convergence              | Tighten/repair sampling and FP inputs; do not flood training with failed labels                       |
| New labels do not improve the next round             | FP failures, wrong paths, stale `record.dpgen`, data not appended              | Verify stage completion, merged systems, and restart identity                                         |

## Data checks

Check units and signs, atom ordering, `type_map.raw`, duplicate structures,
energy/force/stress consistency, abnormal magnitudes, and FP convergence.
Ensure train/validation data cover every composition, geometry, and target
temperature/pressure rather than only adding more near-duplicates.

## Training checks

Compare train versus validation curves, loss components, learning-rate decay,
step budget, batch size, model capacity, ensemble seeds, and backend/version
compatibility. Change one major factor at a time and retain the previous curve
for comparison.

## Sampling and labeling checks

Inspect MD stability, temperatures, pressures, ensembles, `nsteps`, `trj_freq`,
and trust-window placement. Confirm selected structures are labeled successfully,
new data are included in the next training systems, and old/new data remain
balanced. Record conclusions in the project journal; never retune silently.
