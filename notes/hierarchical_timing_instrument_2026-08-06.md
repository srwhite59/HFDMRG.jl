# Hierarchical Timing Instrument

Date: 2026-08-06

Status: private instrumentation implementation validated on a compact fixture;
package and verifier results are recorded before commit.

## Outcome

HFDMRG now has one task-local hierarchical wall-clock instrument, ported and
compacted from GaussletBases `TimeG`. Collection and live output are disabled by
default. When enabled, repeated phase labels merge with exact call counts, and
every node reports inclusive elapsed time and self time after subtracting its
recorded children.

The sweep tree covers initialization, backend initialization, window
construction, local SCF, Fock builds, eigensolves, environment SVD, absorption,
and frozen occupation audits. No public solver keyword, backend operation, or
numerical policy changed.

## Diagnostics reconciliation

The following `_RunDiagnostics` fields were deleted because `TimeG` supersedes
their wall-clock or call-count role:

- `mode` and its timing-only mode;
- `eigsym_seconds` and `eigsym_calls`;
- `aufbau_seconds` and `aufbau_calls`;
- the per-window derived `eigsym_seconds` entry.

`aufbau_bytes` was also deleted. Allocation instrumentation is intentionally
outside this wall-clock-only facility, and retaining one isolated `@timed` byte
counter would preserve the same ad hoc accounting this change removes.

The retained correctness/observability fields are `spectra`, `windows`,
`frozen_rebuilds`, `frozen_max_generations`, `frozen_live_generations`,
`frozen_recurrence_error`, `frozen_zero_up_windows`, `frozen_zero_dn_windows`,
`frozen_zero_both_windows`, and `frozen_rebuild_projector_error`.

The historical dead-field report does not match current source:

- `frozen_recurrence_error` now has a real detailed-mode writer in
  `_frozen_store_growth!`, comparing the recursive direct field, both exchange
  fields, and scalar energy with direct reconstruction. Its numerical test is
  therefore retained.
- frozen absorption/window byte counters are absent from the live structure;
  they had already been removed.
- every other retained field has a writer and a live numerical or lifecycle
  consumer. No vacuous assertion was found in the current tree.

## Additive one-sweep tree

The deterministic one-thread common-H UHF fixture had 14 windows and four
forced local iterations per window. A representative enabled run produced:

```text
initialization: elapsed=0.000124791 s, self=0.000010166 s, calls=1
  backend initialization: elapsed=0.000002500 s, calls=2
  environment SVD: elapsed=0.000018623 s, calls=7
  absorption: elapsed=0.000093502 s, calls=7
sweep: elapsed=0.001069083 s, self=0.000045961 s, calls=1
  window construction: elapsed=0.000012248 s, calls=14
  local SCF: elapsed=0.000903584 s, self=0.000154297 s, calls=14
    Fock build: elapsed=0.000118414 s, calls=70
    eigensolve: elapsed=0.000630873 s, calls=112
  environment SVD: elapsed=0.000026166 s, calls=14
  absorption: elapsed=0.000081124 s, calls=14
```

The sweep children plus sweep self time equal `0.001069083 s`; the local-SCF
children plus local-SCF self time equal `0.000903584 s`. Recursive tests enforce
this identity rather than relying on rounded rendered values.

Each window performed `n=4` local updates and exactly `70/14=5=n+1` backend
Fock builds. UHF performed `112/14=8=2n` spin eigensolves.

## Inertness and overhead

The same authenticated runner was executed from detached production
`main@ac5acee` and the instrumentation branch. With timing disabled:

- energy was exactly `2.012179501149096` on both trees;
- alpha and beta orbital SHA-256 values were both
  `5d5001d8a008f89ba07ea519769c33141242d3e8e6b323e031dd2fe7f18b3413`;
- the sweep-energy trajectory was exactly `[2.012179501149096]`;
- warmed allocation was exactly `2,498,640` bytes on both trees;
- median wall time was `1.209292 ms` on main and `1.125167 ms` with disabled
  instrumentation. This small reverse difference is ordinary timing noise; it
  establishes no performance claim.

With timing enabled, the warmed median was `1.187334 ms`, 5.5% above the
disabled branch median. Allocation was `2,584,016--2,584,096` bytes, about 3.4%
above the disabled solve. Numerical output remained exactly equal. The enabled
overhead is diagnostic and is not paid by default solves.

Runner:
`/Users/srw/dmrgtmp/hfdmrg_hierarchical_timing_20260806/compact_timing_gate.jl`

Runner SHA-256:
`b3a32898295b30098961de4dcbeeb8445f072fe130b0777ae19e2f7f02aeb182`

## Scope

This pass does not add allocation profiling, a public timing keyword, solver or
backend API, or any numerical behavior. It runs no scientific H-chain, WL,
history, or frozen-occupied calculation.

The next useful measurement is one compact frozen-occupied workflow with exact
Fock-call counts and the additive tree enabled, before any large-system timing.
That measurement is not part of this implementation.

## Validation

- The complete package suite passes under `Pkg.test()`.
- The fixed/ragged cached-versus-projection verifier passes. Its maximum window
  Fock difference is `5.10702591327572e-15`, and its maximum sweep-energy
  difference is `1.0658141036401503e-14 Ha`.
- `git diff --check` passes.
