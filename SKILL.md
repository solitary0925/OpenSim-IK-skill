---
name: opensim-vicon-ik
description: Run and review marker-based OpenSim IK from dynamic Vicon C3D and a scaled original-knee CGM24 model. Includes portable preflight, blocked-input reporting with partial-batch continuation, full-trial error locations, timing maps and marker overlays. Does not perform scaling, ID, RRA, SO, CMC or JRA.
---

# Vicon dynamic trials to OpenSim IK

Use a Python environment importing `opensim`, `numpy` and `matplotlib` (when diagnostic plots are enabled). Inputs are an already scaled subject model and matching dynamic C3D files; force data and VSK are not required for marker-based IK. Use the original single-rotation-knee model unless the user explicitly changes this policy. Never substitute the `_knee3rot` version or modify knee axes automatically.

## Inspect before execution

Confirm subject/session and dynamic trial list, sampling rate, units, time span, model marker definitions and C3D point completeness. Do not treat every C3D point channel as a physical marker. Use the included name-alignment helper; preserve raw C3D labels and rename only a working model copy. Block the affected trial on missing/ambiguous model-to-data correspondences and report extra experimental points requiring manual model additions if used. The alias table changes names only; do not move model points or change body assignments during IK.

## Unobservable coordinates

The user authorized locking coordinates without marker observations. For the validated CGM24 profile, all foot markers belong to calcaneus, and torso/toes have no marker descendants. The runner tests marker-position sensitivities at three poses; it locks zero-sensitivity MTP/lumbar coordinates at the input model defaults, with zero default speed, on `ik_model.osim`. It records values as assumptions, not measured angles. It refuses unexpected unobservable coordinates or remaining rank deficiency for review. A body without its own marker may still be constrained through downstream markers; **do not lock ankle/subtalar or a joint merely because its body has no direct marker**. Source models and previously generated scaling artifacts remain untouched.

## Execute

Prefer the portable JSON entry point. Copy `assets/config.example.json` into the study directory; put paths, time window, marker weights/map, rotation, exclusions and visualization limits in that file. All relative paths resolve from the configuration location. Read [configuration and timing contract](references/configuration.md) when configuring a study or diagnosing preflight failures.

```text
python scripts/run_config.py --config study.json --check-only
python scripts/run_config.py --config study.json
python scripts/run_config.py --config study.json --review-only
```

Preflight checks runtime, model/knee/observability, marker names/completeness, units, windows, original C3D frame metadata, rotation, weights, exclusions and output/reuse policy before any IK. Check-only writes nothing unless `--report-dir` is supplied. Normal execution retains resolved configuration, source configuration SHA256 and preflight JSON/Markdown. Review-only upgrades diagnostic sidecars for existing results without solving or retroactively applying new solver settings. Do not ship personal study configurations/data in the shared skill.

## Blocked inputs and partial batches

The unified configuration entry and batch CLI preflight each trial before starting its solver. A trial with unreadable data, unresolved marker names, missing selected samples, inconsistent timing, invalid window or conflicting existing-output policy is `blocked_preflight`; record it and continue other runnable trials. Environment/model/observability or global configuration/exclusion errors stop the whole invocation (`blocked_global` where trials have been enumerated). Never bypass these checks to maximize completion counts. A blocked input is not a user-discarded trial, a solver failure or a completed trial with high errors.

Produce `blocked_inputs.json`, `输入受阻明细.md` and `输入受阻帧范围.csv`. Report all unresolved model names/candidates, all affected marker channels, full-trial and selected-window missing counts and contiguous ranges. Use full-source sample ordinal/total, original acquisition frame labels and adapter seconds as distinct fields; only fill acquisition labels after timing verification. List fully observed ranges for diagnosis, without automatically choosing one. Missing samples wholly outside an explicitly configured solver window do not block that window. Never interpolate, relabel anonymous points, move model markers or introduce additional joint locks to work around a block.

Batch accounting/selection must separately retain `blocked_preflight`, `blocked_global`, `failed_or_unavailable`, user exclusions and completed-trial review classifications. Skip blocked solvers and blocked-trial plots. Preserve historical files, but do not audit/reuse an old MOT as a current completed result when its trial is currently blocked. A repaired input must pass a new preflight; normal reuse/overwrite rules still apply. Report affected subject/trial IDs and concise reasons in the conversation, alongside counts; use the report files for long missing-frame details rather than loading them all into context.

`ready` means every non-excluded trial passed; `can_run` means at least one trial is runnable with no global error. A mixed batch continues even when `ready=false`. Check-only exits 2 for any block; normal batch execution exits 2 after producing partial results if any trial was blocked or failed. Exit 0 means no block/failure (numerical review flags can still exist). Review-only does not solve blocked trials. See [blocked-input outputs and execution contract](references/configuration.md#blocked-inputs-and-execution-status).

Single trial:

```text
python scripts/run_ik.py --subject SUBJECT --model scaled.osim --c3d dynamic.c3d --geometry-dir Geometry --output ik_results/trial
```

Model marker weights default to 1. An explicit `--weights` JSON can assign positive weights by original model name or matched C3D label. Do not automatically delete markers/downweight anatomical points to lower errors. No automatic filtering, interpolation or static-model recalibration. Default rotation is model XYZ = lab (-X,Z,Y); verify for a changed laboratory/session. `--start` / `--end` must lie in the actual data interval, not snap silently outside. `--marker-map` provides confirmed original-name -> C3D-name mappings. A nonempty output requires authorized `--overwrite`.

After a successful representative trial and output verification, batch the same policy:

```text
python scripts/batch_ik.py --subject SUBJECT --model scaled.osim --data session_folder --geometry-dir Geometry --output ik_results --workers 2
```

Batch includes top-level C3D files except names containing `static`. Inspect this selection first; it is not proof of experimental condition or complete gait cycles. `--reuse-existing` requires matching model/C3D hashes, implementation hashes, OpenSim version, weights/map/rotation and time/lock/accuracy policy plus required artifacts. Legacy or changed-policy results require review or explicitly authorized replacement, not silent reuse. `--overwrite` is only for authorized replacement of generated results. Failures retain logs and do not stop unaffected trials. Batch forwards configured weights/map/rotation/time windows to every included trial; shared windows must fit each trial. Split different window policies into separate configs.

## Mandatory post-IK review and user exclusions

`run_ik.py` automatically audits saved results with `review_ik.py`; `batch_ik.py` additionally creates a batch selection/review summary after all trials finish. Existing results can be audited without solving again:

```text
python scripts/review_ik.py --results ik_results --exclude-list excluded_trials.json
python scripts/review_ik.py --trial ik_results/trial
```

The excluded-trial file is a JSON array of **user-confirmed discarded C3D basenames**. Batch execution accepts `--exclude-list`, or loads `excluded_trials.json` from its output folder when present. Validate that exclusions belong to this batch and display them; skip their IK execution and subsequent scientific analysis. Never infer exclusions from a failed threshold or a swap hypothesis, never silently include previously discarded trials, and do not delete source data or old outputs. No S02-specific filenames are embedded in this skill.

The review verifies CSV point/RMS/Max consistency, native error time/frame agreement, finite coordinates, persisted coordinate locks, and source provenance. It localizes contiguous per-point Max>40 mm and frame-RMS>20 mm episodes, writes their start/end adapter time, original source frames, frame count and peak value to `异常时间段.csv`, and records angular-change screens. Report whether anomalies occupy a few frames or persist through a trial; do not automatically reject the whole trial for a short event.

Always locate abnormalities within the **full source trial**: e.g. sample 10–20 / 200, inclusive, plus sample-index percentages and original acquisition frame labels separately. The sample ordinal starts at 1 at the source C3D's first stored sample; it does not reset to 1 for a cropped IK window and is not gait-cycle percent. Show computed IK coverage and merge overlapping point/angle flags when counting unique abnormal samples. If source timing is unverified, label full-trial position unavailable rather than substituting an IK-window denominator. Add a full-trial error timeline for included flagged trials and show positions in batch/per-trial Markdown, CSV/JSON and overlay titles/index.

Optional per-trial `analysis_windows` define an explicitly confirmed later analysis/retention interval using inclusive full-source sample ordinals and a recorded reason. This is a review scope, separate from the solver's `window` seconds. Mark each episode inside, outside or crossing that interval; report unique abnormal frames and RMS/Max for the retained interval separately while preserving whole-trial metrics/status. Outside-only errors are background for the configured downstream analysis, not grounds alone to discard that trial. Persist full evidence, never silently crop/relabel/reweight data, and never re-include user-discarded trials. Without a confirmed interval/event definition, do not infer force-plate contact, complete gait cycles or dispensable lead-in/out samples from ordinal position or large errors alone. Plate-contact phase and a complete gait cycle must not be treated as interchangeable. See the [configuration contract](references/configuration.md#analysis-interval-and-trial-position).

For included `needs_review` trials, automatically generate deterministic model/experimental marker overlays under `异常帧可视化/`, with a Markdown index. Select peak frames from abnormal episodes and angular jumps, prioritize Max error and cap at three frames/trial by default (configurable 1–20). Show three projections of all tracked points and a local view around the worst residual, with point names, time, original frame, RMS and Max. Verify reconstructed residuals against saved CSV before plotting. Preserve all original outputs and report plotting failures explicitly. Existing user-excluded trials retain error evidence but get no new plots by default. The fixed local script needs no image-generation calls and does not require loading every PNG into the conversation; visually inspect a representative output when changing the renderer.

Always preserve C3D acquisition frames and window offsets in `timing.json` and `time_mapping.csv`; newly generated TRC frame labels retain original source frame numbers. Prefer verified TRIAL frame parameters over the separately retained header range. Keep adapter time unchanged through window selection. Store sample offset, point rate, analog samples per point frame and the declared `(source_frame-1)/rate` time convention. Review of older results adds mapping sidecars without rewriting their TRC/MOT. Unsupported processor formats or inconsistent/unverified timing must be reported, never guessed. Frame-number time is not a wall clock or proof of external force synchronization; a future ID workflow must verify force timing/units/axes/plate conventions separately.

For eligible high-error feet, diagnose the FMH/VMH swap hypothesis with a proper-rotation rigid fit of five same-body foot points using the selected experimental TRC. Only report a hypothesis when original rigid RMS>15 mm, swapped RMS<50% of original, and improvement>10 mm. These are engineering screens, not anatomy/label validation. Skip the diagnostic if required points are missing or span different bodies (e.g. a multi-segment foot). Do not make any actual swap, move points, change weights or crop samples without a separately confirmed correction.

Per-trial review outputs: `核查报告.md`, `review.json`, `异常时间段.csv`, timing sidecars, and diagnostic overlays when flagged. Batch outputs: `核查汇总.md`, `review_status.json`, `trial_selection.json` and blocked-input reports. Classification is `blocked_preflight`, `blocked_global`, `excluded_by_user`, `needs_review`, `no_numerical_flags`, or `failed_or_unavailable`. Preserve diagnostic evidence for existing excluded results, but exclude them from the counts of remaining usable/review trials. **No numerical flags means only that these tests did not trigger; never call it anatomically validated or guaranteed problem-free.** The user-confirmed selection manifest governs downstream analysis, while raw batch metrics remain a historical solver record.

## Validate and report

Verify output time/frame count, model reload and unchanged source SHA256; all locked output coordinates must equal their defaults, and all coordinates must be finite. Independently reconstruct marker errors from saved `ik.mot` and compare with OpenSim native errors. Report overall RMS (sqrt over all frames/markers), mean and peak frame RMS, global max/time/marker, per-marker errors and coordinate ranges. Keep metrics definitions distinct from static scaling RMS.

The native `*_ik_marker_errors.sto` contains OpenSim error summaries and `*_ik_model_marker_locations.sto` contains predicted marker locations. Independent `marker_errors_mm.csv` and `metrics.json` facilitate comparison. Report any frames RMS>20 mm or max>40 mm and the engineering screen for >20 degree frame-to-frame angular changes; these are review prompts, not acceptance criteria. Inspect movement/angles before calling the study validated. Locked MTP angles cannot answer a toe-motion research question. Do not infer RRA/SO/CMC/JRA readiness from good IK error alone.

If a few foot points have persistent very large errors while other points fit, check the compatibility of experimental within-foot geometry with the model. A proper-rotation rigid registration can compare the observed labeling with a specific FMH/VMH swap hypothesis. Label this as a diagnostic hypothesis, preserve the original data and IK results, and require confirmation of the physical marker definitions before applying a swap. Never change labels solely to improve fit. S02 testing found this pattern in some Run10 trials; it is not a universal swap rule.

Outputs include the locked working model, aligned model, TRC, IK configuration/MOT, native error/marker-location files, independent errors, lock/observability metrics, name-alignment reports and a Chinese Markdown report. Do not publish subject data or results when sharing the skill.

Official references: [IK workflow](https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53084226/Overview%2Bof%2BOpenSim%2BWorkflows), [IK best practices](https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim33/pages/53674534), [IK API](https://opensim-org.github.io/opensim-moco-site/docs/1.1.0/html_user/classOpenSim_1_1InverseKinematicsTool.html).
