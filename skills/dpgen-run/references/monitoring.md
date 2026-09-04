# DP-GEN Monitoring and Troubleshooting

Load this reference after each iteration or whenever accuracy or stage progress
stalls. Here "accuracy" means energy/force train and validation RMSE; model
deviation and FP failure ratio are separate signals.

## Evidence loop

For each iteration, record the parameter/configuration identity, systems and
temperature/pressure conditions, train/validation RMSE, model-deviation
accurate/candidate/failed counts, FP success/failure counts, current stage, and
whether new labels reached the next training input. Report per-system and
per-condition results before any global aggregate.

## Diagnose in this order

| Symptom                                              | First checks                                                                 | Focused adjustment                                                          |
| ---------------------------------------------------- | ---------------------------------------------------------------------------- | --------------------------------------------------------------------------- |
| Train and validation RMSE both stay high             | Units, labels, type map, coverage, descriptor capacity                       | Fix data/schema first; then review fitting net, learning rate, and steps    |
| Train RMSE falls but validation RMSE stalls or rises | Split leakage, duplicates, outliers, validation coverage                     | Re-split or deduplicate; add representative labels before changing capacity |
| RMSE is good but model deviation is high             | Unseen systems/geometries or conditions, initial structures                  | Broaden targeted sampling and labeling the candidate window                 |
| Candidate count is zero                              | Trust thresholds, MD stability, trajectory frequency, model diversity        | Check thresholds and sampling settings before changing training             |
| Candidate count is excessive or FP failures rise     | Thresholds, unstable MD, malformed structures, FP convergence                | Repair inputs/sampling; do not add failed labels                            |
| Stage does not advance or task count is wrong        | record state, schedule length, sys_idx, generated input, submission metadata | Reconcile state and actual tasks; do not restart blindly                    |
| Effective steps differ from requested steps          | Reuse/iteration override in generated input                                  | Compare root and task-level values; fix the override deliberately           |

## Data and label checks

Check units and signs, atom ordering, type maps, duplicates, abnormal
magnitudes, and train/validation coverage. Every training frame must have
finite, parseable energy and force; include virial/stress only when configured
and verified. Reject incomplete, duplicated, or unconverged labels before
merging them into training data.

## Training checks

Compare train and validation curves, loss components, learning-rate decay, step
budget, batch size, model capacity, ensemble seeds, and backend/version
compatibility. A good global metric must not hide a failing system or condition.
Change one major factor at a time and retain the previous evidence.

## Sampling and decision checks

Inspect MD stability, temperatures, pressures, ensembles, `nsteps`, `trj_freq`,
trust-window placement, and per-condition distributions. Distinguish a model
problem from missing coverage: good RMSE with high deviation usually calls for
targeted sampling, while high RMSE on covered data calls for data or training
review.

Monitoring is an evidence-based decision stage, not only status reporting.
Before extending sampling, changing thresholds, increasing steps, or restarting,
record the evidence, exact parameter delta, expected effect, and controller
state. Afterward record the new task identity and outcome; never retune or
restart silently.
