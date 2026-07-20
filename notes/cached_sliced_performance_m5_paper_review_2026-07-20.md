# Paper-Manager Review: Cached Sliced Milestone 5

Date: 2026-07-20
Reviewed through commit: `ead4a4b`
Decision: **cached implementation accepted; frozen 50-sweep compound gate remains failed; bounded M5b convergence characterization authorized**

## Bottom Line

The cached sliced implementation has passed its Be-scale numerical and
performance acceptance. It reproduces the identically aligned projection
route and independent full-Fock calculations, uses only the intended local
route, and changes a roughly 30.5-second projected sweep into a roughly
0.18-second cached sweep. There is no evidence of a cache, Fock-convention,
partition, or returned-energy defect.

The frozen Be run did not satisfy every convergence and oracle-fingerprint
gate within the predeclared 50-sweep cap. That result must remain recorded as
a failed gate; it should not be relabeled or hidden. It is, however, a solver
convergence-tail result rather than a rejection of M1--M4. At sweep 50 the
energy is already within `3.73e-10 Ha` of the oracle and the occupied minimum
singular value is `0.999999999632`, while the independent energy identity is
satisfied to `1.78e-15 Ha`. The state is on the intended broken-symmetry
branch but is not yet stationary to the desired standard.

## Accepted Evidence

- Projection/cached one-sweep energies differ by `1.17e-13 Ha`; occupied
  projectors differ by at most `2.05e-10`.
- Explicit physical, projection, and cached total Fock matrices agree within
  `2.27e-13`.
- Every forced-four cached sweep uses 98 windows, 490 Fock calls, and zero
  fallback or projection windows.
- Warmed first and second cached sweeps take medians of `0.171570 s` and
  `0.132383 s`, allocating about `193 MiB` each.
- The aligned projected sweep takes `30.47646 s`; the corresponding cached
  sweep takes `0.18147 s`, a measured end-to-end speedup of about `168x` for
  this frozen Be fixture.
- Fresh-process peak RSS is `855.0625 MiB`.
- All 182 package assertions and the fixed/ragged verifier pass. M5 changed no
  production or test source.

These numbers are paper-usable for this exact fixture after normal timing
provenance is carried with them. They are not yet a scaling law.

## Interpretation Of The Failed Gate

The last ten recorded sweeps are smooth. Geometric fits over sweeps 40--50
give approximate crossing points:

| Quantity | Sweep-50 value | Gate/reference | Estimated crossing |
|---|---:|---:|---:|
| Energy change | `9.64e-11 Ha` | `1e-11 Ha` | sweep 60 |
| Projector distance | `3.83e-5` | `1e-5` | sweep 62 |
| Minimum-singular-value defect | `3.68e-10` | `1e-10` | sweep 56 |
| S2 deviation | `3.19e-5` | `1e-6` | sweep 80 |
| Z-spin deviation | `1.47e-4` | `1e-6` | sweep 92 |
| Full-Fock residual | `3.60e-5` | oracle `1.32e-8` | about sweep 123 |

These are extrapolations, not results. They show why merely increasing the
normal cap to 60 would answer the energy-convergence question but not the
stationarity and exact-oracle-tail questions. They also expose a paper-relevant
method issue: energy convergence is substantially faster than full Fock
stationarity for this distant starting state.

## M5b Authority

Run one bounded continuation study with no production or committed test
changes.

1. Minimally parameterize the existing validation runner, or use a retained
   validation-only mode, so the maximum convergence sweeps and outer energy
   cutoff are recorded command-line inputs. Do not add a framework. Keep the
   runner change at no more than 20 net added lines, preferably offset by
   simplifying existing validation code.
2. Run the ordinary policy from the original frozen distant seed with
   `cutoff=1e-11`, `scf_cutoff=1e-12`, and a cap of 80. Record the actual
   convergence sweep and all existing diagnostics.
3. Run a diagnostic trace from the same original seed through 140 sweeps with
   the outer energy stopping cutoff disabled, while retaining the same local
   SCF policy. This is not a proposed production setting; it measures the
   asymptotic tail.
4. Confirm that the first 50 rows reproduce the existing trajectory to
   numerical precision. Record the first sweep satisfying each existing
   energy, projector, singular-value, S2, and Z-spin gate, plus the residual
   at sweeps 60, 80, 100, 120, and 140.
5. Continue to require 98 local windows and zero cache fallbacks. Stop and
   report if energy ceases to decrease smoothly, the run leaves the saved
   branch, or the first 50 sweeps do not reproduce.

No projection rerun, timing rerun, package test rerun, source repair, damping
change, acceleration implementation, or gate rewriting is needed for M5b.
Raw output remains under `~/dmrgtmp`; commit only a compact runner amendment,
reduced trajectory evidence, and status update.

## Paper Consequence

The cache optimization is ready to support the paper's sliced-backend timing
study. The Be trajectory simultaneously becomes a useful convergence example:
the energy and occupied space approach the known branch smoothly, but a
strict stationary solution takes many more sweeps than the energy alone would
suggest. Acceleration methods should be tested later as a separate algorithmic
study, not mixed into acceptance of the cache implementation.

-- hfdmrg-paper-manager@faraday30.ps.uci.edu
