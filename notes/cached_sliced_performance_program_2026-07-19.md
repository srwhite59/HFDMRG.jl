# Cached Sliced-Backend Performance Program

Date: 2026-07-19
Owner: `hfdmrg-manager`
Scientific acceptance: `hfdmrg-paper-manager`

## Purpose

Make the corrected sliced HFDMRG route usable for the frozen radial-`Y_lm`
Be fixture without changing its Hamiltonian or the HFDMRG optimization
algorithm. This is software-enabling work, not a new scientific method.

Start from:

```text
f097ce49f661a9a6881131a4fc77af8d00d3b9f7
Fix sliced-basis Coulomb and exchange contractions
```

Preserve the unrelated untracked `JuliaStyle.md`.

## Workspace And Permission Boundary

All durable work for this program stays in this HFDMRG repository. Maintain
one current status report at:

```text
notes/cached_sliced_performance_status_2026-07-19.md
```

Use machine-local scratch under:

```text
~/dmrgtmp/hfdmrg_cached_sliced_20260719/
```

Use the normal machine-local `~/.julia` depot. Do not create depots, caches, or
environments inside Dropbox. Do not write into `~/Dropbox/Papers/HFDMRG`; the
paper manager will read this repository directly and maintain paper records.

Prefer Julia timing (`@elapsed`, `time_ns`, or `BenchmarkTools`) rather than
`/usr/bin/time`. Long jobs must write a known log and PID file. Poll that log
or the active tool session rather than using broad process searches.

## Diagnosed Cost Boundary

The worst defects are sliced-specific:

1. Projection reconstructs a full sliced Fock matrix on every local call.
2. Cached absorption discards and rebuilds global raw sliced caches.
3. Generic block cuts split slices and activate global exchange fallback.
4. Cached Coulomb still scans all global slice pairs per local Fock call.

Two shared-core issues amplify this:

1. Four local microiterations currently make eight backend Fock builds.
2. Generic `blocksize` partitioning has no physical slice boundaries.

Density-density and target-residual backends already update retained state
incrementally. Do not generalize the sliced rewrite into a new cache framework.

## Working Policy

Work on a dedicated branch. Make narrow commits at milestone boundaries and
update the one status report with the current commit, measurements, unresolved
issue, and exact next action. Proceed between milestones without requesting
pass-by-pass approval.

Stop for review on:

- numerical disagreement with the projection/explicit oracle;
- a required public or backend/core API expansion not anticipated here;
- persistent unfavorable scaling after profiling; or
- a scientific/solver-policy choice.

The projection backend remains the numerical oracle. Do not optimize by
changing damping, local iteration count, convergence policy, or the Be
Hamiltonian.

## Milestone 0: Design And Cost Checkpoint

Identify the exact retained cache objects and update equations. Compare:

- contracting only newly active source slices; and
- recursively transforming/updating the previous retained cache.

Choose the smallest exact design that removes the bad scaling. Record the
proposed files, line budget, deletion plan, cost model, and risks. Make no
production source edit before this checkpoint is reviewed.

## Milestone 1: Slice-Aligned Partitioning

Provide a general way for a structured backend or caller to specify contiguous
physical block boundaries. Radial Be should have 52 complete nine-orbital
blocks, 98 sweep windows, and zero split-slice fallbacks.

Gates:

- exact contiguous coverage and cuts only at `SliceLayout` boundaries;
- unchanged ordinary numeric-`blocksize` behavior;
- cached/projection parity using identical ranges; and
- no Be-specific hard-coded range list.

## Milestone 2: Incremental Cached Absorption

Ordinary aligned absorption must update retained caches from the previous
state and newly absorbed block rather than rebuild all raw slice pairs. A
direct rebuild may remain for boundary initialization and tests.

Gates:

- direct-rebuild parity after left/right absorption, changing ranks, and
  fixed/ragged slices;
- arbitrary symmetric-density RHF/UHF Fock parity;
- scale-aware disagreement no worse than `1e-11`;
- no post-boundary full raw-cache rebuild; and
- measured scaling with slice count.

Earlier estimates of `16 s` initialization and `32 s` cache work per sweep are
provisional engineering gates, not paper claims.

## Milestone 3: Window-Local Cached Coulomb And Exchange

Aligned local Fock calls must use retained environment, cross, and center
tensors. They must not lift density globally or scan all global slice pairs.

Gates:

- explicit/projection RHF/UHF parity for arbitrary symmetric non-idempotent
  densities;
- frozen Be window error at most `1e-12` where summation order permits;
- no undocumented integral symmetry;
- warmed aligned Fock at least 20 times faster than projection, or a measured
  justification for revising that gate;
- near-flat Fock time with growing total slice count at fixed window ranks;
- effectively allocation-free steady Fock work; and
- zero fallback calls in the accepted Be route.

## Milestone 4: Shared-Core Duplicate Fock Reuse

Do this as a separate backend-general commit after the sliced local kernel is
correct. For `n` microiterations, use `n+1` Fock builds rather than `2n`; four
microiterations should use five rather than eight.

Gates:

- exact call-count coverage;
- common-H and split-H RHF/UHF coverage;
- density-density and sliced trajectory/energy/subspace parity;
- no backend API or damping-policy change; and
- about 35% or better backend-time reduction when Fock work dominates.

## Milestone 5: End-To-End Be Acceptance

Measure initialization, one warmed sweep, and matched convergence separately
from the frozen starting determinant.

Required evidence:

- explicit full-Fock energy consistency;
- orthonormality, idempotency, and occupied-virtual stationarity;
- cached/projection same-range parity;
- wall time, memory method, allocation, Fock/window counts, BLAS threads,
  host, Julia version, and commit;
- block-size/alignment dependence; and
- no material scientific-trajectory change.

The provisional readiness target is a warmed forced-four Be sweep below
`40 s` and cumulative allocation below `2 GiB`. If a target fails, profile the
failed stage rather than changing solver policy.

After this program, density-density scaling on hydrogen chains is a separate
paper experiment. March acceptance, history acceleration, ordering, and
metastability are also separate algorithm studies.
