# Paper-Manager Cross-Host Acceptance: Cached Sliced Milestones 1 And 2

Date: 2026-07-19
Accepted through commit: `6cb865d`
Decision: **Milestones 1 and 2 accepted; Milestone 3 authorized**

## Cross-Host Result

Faraday independently ran the ordinary package entry point at `6cb865d`:

```text
julia version 1.12.6
BLAS: OpenBLAS 0.3.29, ILP64, 10 threads
julia --project=. test/runtests.jl
result: all 95 assertions passed
```

The branch contained only the pre-existing untracked `JuliaStyle.md`.
`git show --check 6cb865d` passed, and the commit changed only the split-H
projector tolerance in `test/runtests.jl` from `1e-10` to the documented
cross-platform `1e-9` gate.

The earlier faraday mismatch was a one-sweep trajectory amplification, not a
cache defect. The split-H energy disagreement was `7.98e-12 Ha`; projector
disagreement fell from `6.85e-10` after one fixed sweep to `1.12e-11`,
`5.65e-14`, and `1.65e-15` after subsequent sweeps. Direct cache and Fock
oracles retain their machine-precision gates.

## Accepted Milestones

Milestone 1's slice-aligned partition and Milestone 2's incremental cached
absorption are accepted. This acceptance includes the corrective
`nblockcenter >= 3` offset, anchored-side eligibility, left/right off-span and
partial-slice fallbacks, density-density partition semantics, deterministic
trajectory fixture, and the previously reviewed numerical and cache-chain
performance evidence.

This is engineering acceptance. It is not yet an end-to-end Be timing or paper
scaling claim.

## Milestone 3 Authority

Execute Milestone 3 only: replace the aligned cached backend's global sliced
Fock scans with the approved window-local Coulomb and exchange kernel.

Required gates remain:

- isolate all nine ordered `LL`, `LC`, `LR`, `CL`, `CC`, `CR`, `RL`, `RC`,
  and `RR` sectors with unsymmetrized interactions and symmetric,
  non-idempotent RHF/UHF densities;
- preserve explicit/projection parity, with the frozen workflow window error
  at most `1e-12` where summation order permits;
- delegate the entire split-window interaction to projection rather than keep
  duplicate global fallback algebra;
- make warmed aligned Fock at least `20x` faster than projection;
- require `time_Fock(S=52) / time_Fock(S=8) <= 1.5` at fixed ranks and center;
- require zero steady allocation in the aligned backend Fock method and zero
  fallback windows in the accepted structured route;
- extend the retained `scripts/bench_sliced.jl` harness to report cache-chain,
  window-assembly, and Fock timing/allocation with batch count and statistic;
- materially delete the superseded global RHF/UHF sliced Fock paths and the
  stale scratch/metadata listed in the status memo; and
- stay within the existing M3 and combined M3/M4 line budgets without adding
  a source file, public type, generic cache framework, or backend API change.

Stop after the M3 commit and complete validation. Do not begin M4 duplicate
Fock reuse or M5 Be acceptance until paper-manager review.

-- hfdmrg-paper-manager@faraday30.ps.uci.edu
