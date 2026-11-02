# Configuration and timing contract

Copy `assets/config.example.json` to the study directory and edit it. All relative paths resolve from that JSON file, independent of the shell working directory. Do not distribute study data, absolute personal paths or exclusion lists in the reusable package.

```text
python scripts/run_config.py --config study.json --check-only
python scripts/run_config.py --config study.json --check-only --report-dir preflight_report
python scripts/run_config.py --config study.json
python scripts/run_config.py --config study.json --review-only
```

Check-only starts no IK and writes nothing unless `--report-dir` is supplied. Runtime imports and temporary timing checks occur in a temporary directory. Normal execution saves resolved configuration, source config SHA256, preflight JSON/Markdown and an action record under `_configuration` (batch) or a sibling `TRIAL_configuration` (single trial). Review-only audits saved results and generates timing/visualization sidecars; it does not retroactively apply new solver settings or rerun IK. A batch review requires `batch_metrics.json`.

A mixed batch now continues runnable trials while separately reporting trial-local blocks. Global environment/model/configuration errors still stop execution. A nonzero exit for a mixed batch is an incomplete-result signal, not proof that all trials failed; inspect the status counts below.

| Field | Meaning |
|---|---|
| `schema_version` | Required integer 1; unknown fields block to catch misspellings. |
| `subject`, `model`, `output` | Required subject label, scaled original-knee model, separate generated-output directory. |
| `data` or `c3d` | Exactly one: directory of top-level non-static C3Ds for batch, or one file. |
| `geometry_dir` | Optional mesh directory; absent meshes warn but do not prevent marker-only IK/overlays. |
| `window.start`, `window.end` | Optional C3D adapter seconds. Null uses each trial's actual boundary; shared explicit bounds must fit every included trial. Split different window policies into separate configurations. |
| `analysis_windows` | Optional per-trial confirmed downstream review/retention intervals; inclusive sample ordinals in the full source C3D, with reason. Does not crop IK or establish contact events. |
| `weights` | Object mapping original model names or matched C3D labels to finite positive weights; omitted markers have weight 1. |
| `marker_map` | Confirmed original model marker name -> exact experimental label. Name changes only. |
| `rotation` | Nine row-major numbers, proper orthogonal matrix, lab-to-model; default model XYZ=lab(-X,Z,Y). Confirm laboratory convention before using elsewhere. |
| `exclusions` | Unique user-confirmed C3D basenames. Existing output/excluded_trials.json is inherited if omitted. Explicit config cannot silently re-include persisted exclusions; reconcile the study manifest deliberately. |
| `workers` | 1 or 2 subprocesses. |
| `reuse_existing` | Batch only. Reuse requires exact implementation, source, OpenSim version and solver policy fingerprints plus required files. Legacy results can be reviewed but are not silently reused for a solve. |
| `overwrite` | False by default; enable only when replacing generated results is authorized. |
| `visualization.enabled` | True by default; plot only included `needs_review` trials. |
| `visualization.max_frames_per_trial` | 1–20, default 3. Select peak-error frames from contiguous abnormal episodes and angular jump frames, then prioritize by Max error. Capping is documented in the image index. |

Preflight checks runtime (`opensim`, `numpy`, and `matplotlib` when plots enabled), source model/knees/observability, marker alignment, finite nonzero selected samples, units, proper rotation, weights, exclusions, requested time window, C3D original frame metadata and output/reuse policy. It does not establish anatomical accuracy or force synchronization.

## Blocked inputs and execution status

`run_config.py` and `batch_ik.py` run per-trial preflight before submitting IK jobs. `run_ik.py` remains the low-level solver with its own input checks; use the unified entry for reusable blocked-input reporting.

| Status | Behavior |
|---|---|
| `runnable` (preflight) / `completed` (execution) | Submit the trial, or reuse only an exactly compatible result. Completed results are separately classified by numerical review. |
| `blocked_preflight` | Do not start this trial. Causes include missing/ambiguous names, selected-window gaps, invalid/unreadable C3D/timing, window errors, or incompatible/nonempty output. Continue other runnable trials. |
| `blocked_global` | Stop the invocation for environment, model/knee/observability or global config/exclusion problems. Trial rows can be absent if failure occurs before enumeration. |
| `excluded_by_user` | Only from the explicit user-confirmed exclusion list; skip data parsing/IK and preserve that decision. |
| `failed` (execution) / `failed_or_unavailable` (review) | Solver failed or required output is unavailable; retain logs and continue independent trials. |
| `needs_review` / `no_numerical_flags` | Completed-trial numerical screening results, distinct from input blocks and user exclusions. |

Preflight JSON keeps `errors` (all errors), `global_errors`, `counts`, `ready` (all passed), `can_run` (at least one runnable trial, no global error), and per-trial `status`, `errors`, `input_details`. A partial batch may have `ready=false, can_run=true`; normal execution still runs the valid subset. All-blocked batches produce accounting/reports without invoking IK. Input blocks are never appended to exclusions automatically.

Outputs at the batch root (and in the saved preflight snapshot):

- `blocked_inputs.json`: structured blocked rows, errors, unresolved names/candidates, all missing marker channels, full-source and selected-window counts/ranges, verified acquisition frame bounds and complete-marker ranges.
- `输入受阻明细.md`: readable per-trial cause/count summary.
- `输入受阻帧范围.csv`: per-marker contiguous ranges, with `scope=full_source` or `selected_window`. Do not sum the two scopes: selected gaps are also part of full-source gaps. `unresolved_name` rows identify missing/ambiguous correspondences without inventing sampled positions.
- `batch_metrics.json`: execution status for every selected trial; the current `review_status.json`/`trial_selection.json` separately retain blocked, excluded, failed and completed states.

Single-trial configuration writes preflight/block reports under the sibling `TRIAL_configuration` directory and does not start a blocked solver. Check-only writes reports only when explicitly given `--report-dir`. An unreadable C3D has an error row even when frame/marker detail cannot be recovered. A timing failure leaves acquisition labels empty, while sample ordinals can still describe readable rows.

Missing sample detection uses the validated finite/nonzero marker check. It records every affected marker, rather than stopping after the first. Full-source sample ordinals begin at 1 and do not reset for selected windows; original acquisition frame labels and adapter times are separate columns. Missing samples outside the explicitly selected solver window are documented but do not block that otherwise valid window. Complete-marker windows are diagnostic candidates, not automatic gait-cycle selection or authorization to crop.

Currently blocked or failed records never become completed merely because historical metrics/MOT files still exist. Preserve those historical files, omit them from current completed analysis, and show their existence in `previous_results_preserved`. Review-only does not repair inputs or rerun IK. After an input is repaired, run a fresh preflight; incompatible old results still require the usual authorized replacement or a separate output directory.

Exit codes: **0** when there are no input/global blocks or solver failures; **2** for check-only blocks, a blocked single trial/global failure, or a normal batch that finishes with blocked/failed trials. Numerical review flags alone do not make solver execution fail. A batch with only user exclusions starts no solver and can exit 0.

## Timing outputs

`timing.json` preserves header and TRIAL frame ranges, processor, point rate, analog samples per point frame, selected offset and explicit time basis. `time_mapping.csv` has one row per IK sample:

- `ik_row`: output row, starting at 1.
- `ik_time_s` and `c3d_adapter_time_s`: unshifted adapter time. A selected window starting at 0.2 s still starts at 0.2 s in MOT/TRC.
- `source_sample_index_0`: sample offset from the C3D file's first stored point frame.
- `source_frame`: original frame number, preferring TRIAL:ACTUAL_START_FIELD/END_FIELD when present and verified against row count.
- `source_elapsed_s`: source sample index / point rate.
- `frame1_based_time_s`: (source frame - 1) / point rate, an explicitly declared frame-number convention, not an absolute clock.

TRC frame labels retain the selected original frame numbers. OrigDataStartFrame/OrigNumFrames preserve the full source file start/count. Existing older TRCs are not rewritten during review; time_mapping.csv supplies their original frame correspondence. Abnormal episode CSV adds original start/end frame columns; overlay titles and index include source frame/time.

The dependency-free metadata reader supports Intel and MIPS C3D. Unsupported DEC, inconsistent frame ranges/rates and unverified source provenance block timing verification rather than inventing offsets. Two 16-bit TRIAL words are combined as an unsigned 32-bit value (low + high shifted 16 bits); header and TRIAL values are retained separately. Header-only fallback is explicitly labeled.

Do not infer external force latency or synchronization from analog samples per point frame. A future ID workflow must check analog/force sample timing, units, coordinate rotation and plate/COP conventions against this mapping.

## Diagnostic overlays

`异常帧可视化/index.md` links deterministic PNGs. Blue circles represent experimental markers, orange crosses represent reconstructed model markers, and connecting lines show residuals. Three model-coordinate projections show all markers and a local view around the largest residual. Coordinates are millimeters; the model's positive anatomical directions must be interpreted from its axes. These are marker overlays, without inferred bones, surfaces or automatic corrections. The renderer verifies residuals against saved independent error CSV before exporting. Figures are diagnostic, not publication figures.

## Analysis interval and trial position

Every abnormal episode reports full-source sample range/total (e.g. **10–20 / 200**, inclusive), percent of source sample count, original C3D acquisition frame labels and IK row range/row count separately. Source ordinal 1 is the C3D file's first stored point sample, even when the acquisition frame label is 689. If only a later IK window was solved, the denominator remains the complete source trial and coverage is explicitly reported. Percentages are `sample ordinal / source sample count * 100`, not normalized gait-cycle time. Unverified original mapping leaves these fields unavailable rather than inventing a denominator.

`trial_position.abnormal_ranges` merges overlapping marker/RMS flags and both endpoints of angular jumps, preventing double counting. `异常时间段.csv` retains per-marker episodes and adds ordinal/total/percent, IK row and analysis-relation columns. `异常帧可视化/error_timeline.png` locates RMS/Max on the full trial, with flagged samples and an analysis interval when configured. Uncomputed samples have no error curve.

Example study config addition (replace this illustrative filename/interval with confirmed study values):

```json
{
  "analysis_windows": {
    "dynamic.c3d": {
      "start_frame": 21,
      "end_frame": 180,
      "reason": "User-confirmed retained gait interval from synchronized events"
    }
  }
}
```

Names are case-insensitive C3D basenames; bounds must be integers, start at 1, and lie fully within both the original trial and computed IK coverage. `reason` is required. Each trial can specify its own interval. Empty/omitted means unconfigured. Preflight blocks unknown trials or invalid bounds. The optional standalone `--analysis-windows windows.json` argument to run/review/batch scripts reads just this filename-to-interval object.

`analysis_context` records inside/outside unique abnormal counts, retained-interval RMS/Max and a separate numerical screening classification. Episodes can be `inside_analysis_window`, `outside_analysis_window` or `crosses_analysis_boundary`; angular jumps count both involved samples so a boundary jump is not hidden. Whole-trial metrics/classification and user exclusions remain intact. Outside-only errors can be treated as background for that explicitly configured downstream scope, while input/lock/provenance verification failures still require review. An interval with no flags is not anatomical validation.

No contact detection or actual cropping occurs. A future event/force workflow or explicit study decision must provide the correct retained interval. An early/late abnormality does not itself prove the foot was off the plate; force-contact phase alone does not define a complete gait cycle.

Reference: [C3D TRIAL frame parameters](https://www.c3d.org/HTML/Documents/thetrialgroup.htm), [C3D processor/endianness](https://www.c3d.org/HTML/Documents/parameterheader.htm).
