# Paper-Manager Review: Cached Sliced Milestone 4

Date: 2026-07-19
Reviewed through commit: `65f3c0d`
Decision: **Milestone 4 accepted; Milestone 5 authorized below**

## Bottom Line

The shared solver core now performs one initial Fock build and one fresh Fock
build after each density update. Four local microiterations therefore require
five backend calls rather than eight. The retained post-density Fock matrix is
exactly the matrix needed for the next diagonalization; it is not stale and is
not modified between uses.

No change was made to damping, local convergence, energy evaluation, orbital
expansion, observer notification, or the backend additive-Fock contract. No
correctness, architecture, or minimality blocker remains.

## Independent Faraday Verification

At `65f3c0d`, Julia 1.12.6 with ten-thread OpenBLAS:

```text
julia --project=. test/runtests.jl
    all 182 assertions passed, including all 13 M4 assertions

julia --project=. scripts/verify_cached_vs_projection.jl
    maximum reported Fock disagreement: 3.55e-15
    maximum reported sweep-energy disagreement: 1.78e-15 Ha
```

The committed counting wrapper proves five calls per forced-four window and
two calls for one immediately accepted update. It covers density-density RHF,
common-H UHF, split-one-body UHF, and aligned cached-sliced UHF. Its 69 test
lines are live-contract instrumentation rather than production scaffolding.

To independently test the old-versus-new trajectory, faraday ran one fixed
two-sweep fixture at pre-M4 `b3362f3` and current `65f3c0d` for six routes:

```text
density-density: RHF, common-H UHF, split-H UHF
cached sliced:   RHF, common-H UHF, split-H UHF
```

For all six routes, returned energies, full alpha/beta orbital matrices, and
every two-sweep observer snapshot were bit-for-bit identical. Review outputs
are retained at:

```text
~/dmrgtmp/hfdmrg_m4_faraday_review_old.jls
~/dmrgtmp/hfdmrg_m4_faraday_review_new.jls
```

This closes the only evidence gap in the committed tests, which compare
counted and uncounted current execution rather than carrying an old solver
inside the test suite.

## Scope And Claim Boundary

M4's `+86/-20` accounting is correct. `src/core.jl` shrank by three lines, no
API or source file was added, and the whole performance program remains 213
lines smaller.

Paper-safe now: four forced local updates structurally use five rather than
eight delegated interaction builds, a 37.5% work-count reduction, without
changing the tested trajectory.

Engineering-only: the macmini `0.736171 -> 0.460419 s` timing and
`80.212 -> 50.133 MiB` cumulative-allocation comparison. Those numbers do not
establish complete-sweep, Be, peak-memory, or scientific speedup.

## Milestone 5 Authority

Execute a validation-only frozen radial-`Y_lm` Be campaign. Do not change
production source, tests, solver policy, the Hamiltonian, damping, iteration
counts, or convergence criteria to meet a gate.

### Frozen Fixture

Use the existing packet, copying it machine-to-machine only if necessary:

```text
radial_be_lmax2_v1.ser
SHA256 fba29651c02cde3cdde90850d90ba30a1abd120ebfa8ac8b64cf2e5a92fea629
dimension 468 = 52 radial slices * 9 angular functions
Nalpha = Nbeta = 2
```

The faraday copy is at:

```text
/Users/srw/dmrgtmp/hfdmrg_radial_be_20260718/radial_be_lmax2_v1.ser
```

Verify the checksum before running. Stop if the packet or its frozen distant
seed and current-Hamiltonian UHF oracle cannot be recovered exactly.

### Durable Runner And Evidence

Retain one concise parameterized runner, preferably
`scripts/bench_radial_be.jl`, with a preferred cap of 220 lines. It may accept
the packet path and output directory but must not embed a machine-specific
path. Do not add a framework, public API, or committed test.

Keep raw logs, serialized checkpoints, and PID files under
`~/dmrgtmp/hfdmrg_cached_sliced_20260719/m5_be/`. Commit only the runner,
compact reduced TSV/Markdown evidence, and replacement status text. Record
the exact command, commit, host, CPU, Julia and BLAS threads, packet checksum,
warmup/trial count, timing statistic, allocation method, and memory method.

Use Julia timing and allocation measurements. Record absolute process peak
RSS using `Sys.maxrss()` in a fresh-process run, clearly distinguishing it
from cumulative Julia allocation.

### Required Runs

1. **Fixture/oracle check.** Recompute the starting and oracle determinant
   energies from the physical full sliced interaction. Confirm dimensions,
   occupations, ordering, and the frozen oracle branch.
2. **Same-range parity.** From the same distant seed, run one forced-four
   sweep with projection and cached backends using the identical
   `block_partition=layout`. Compare energy and occupied projectors, and check
   representative explicit full-Fock windows.
3. **End-to-end timing.** Warm the current cached method, then measure at least
   three same-process trials. Separate initialization, first sweep, and a
   steady second sweep. Report window, Fock, initialization, incremental,
   local/projection, and fallback call counts.
4. **Matched convergence.** Run the aligned cached route from the frozen
   distant seed with the established physical settings
   (`cutoff=1e-11`, `scf_cutoff=1e-12`) until convergence or 50 sweeps. Record
   every sweep's energy, independently recomputed energy, full-Fock
   occupied-virtual residual, occupied angle to the oracle, orthonormality,
   idempotency, and branch fingerprints.
5. **Alignment control.** On the current commit, compare one sweep of legacy
   numeric `blocksize=9` projection, aligned projection, and aligned cached
   execution from the same seed. Keep this a geometry/runtime control rather
   than an algorithm ranking.

### Acceptance Gates

- representative explicit/projection/cached Fock disagreement at most
  `1e-12` where summation order permits;
- independently recomputed returned energy disagreement at most `1e-10 Ha`;
- aligned projection/cached one-sweep energy disagreement at most `1e-9 Ha`
  and occupied-projector norm disagreement at most `1e-8`;
- orbital orthonormality and density idempotency errors at most `1e-11`;
- exactly 98 windows per aligned sweep, five Fock calls per forced-four
  window, and zero cached fallback windows;
- warmed forced-four cached first and steady sweeps each below `40 s`;
- cumulative allocation for one measured cached sweep below `2 GiB`; and
- recovery of the established broken-symmetry Be branch without a material
  scientific-trajectory change.

If a timing, allocation, parity, or scientific gate fails, profile and report
the failed stage. Do not repair source or change solver/scientific policy in
M5. Stop after the evidence commit for paper-manager review.

-- hfdmrg-paper-manager@faraday30.ps.uci.edu
