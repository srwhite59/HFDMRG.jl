# Paper-Manager Review: Cached Sliced Milestone 3

Date: 2026-07-19
Reviewed through commit: `22f7c07`
Decision: **Milestone 3 accepted; Milestone 4 authorized under the amended budget below**

## Bottom Line

The aligned cached sliced backend now evaluates Coulomb and exchange from the
nine local left/center/right sectors without a global density lift or global
slice-pair scan. Unaligned or split windows delegate their complete interaction
calculation to the projection backend. No correctness or architecture blocker
was found.

Milestone 3 is engineering acceptance, not yet an end-to-end Be or paper
scaling result. The several-thousand-fold ratios compare warmed interaction
kernels on already constructed windows. Cache construction and window
assembly are measured separately; one-body construction, eigensolves, SVDs,
and the complete sweep remain for Milestone 5.

## Independent Faraday Verification

Environment:

```text
host: faraday30.ps.uci.edu
commit: 22f7c07
Julia: 1.12.6
Julia threads: 1
BLAS: OpenBLAS 0.3.29, ILP64, 10 threads
```

Results:

```text
julia --project=. test/runtests.jl
    all 169 assertions passed

julia --project=. scripts/verify_cached_vs_projection.jl
    maximum reported Fock disagreement: 3.55e-15
    maximum reported sweep-energy disagreement: 1.78e-15 Ha

julia --project=. scripts/bench_sliced.jl
    S52 RHF kernel speedup: 3044x
    S52 UHF kernel speedup: 3796x
    cached S52/S8: 0.979 RHF, 0.991 UHF
    cached steady Fock allocation: zero
    S52 cache initialization: 0.09685 s
    S52 two-way cache sweep: 0.20008 s
    S52 cached window assembly: 0.273 ms
```

These independent values differ in ordinary host-dependent timing from the
macmini record but pass every predeclared numerical, speed, scaling,
allocation, and routing gate by a wide margin. `rh310l.ps.uci.edu` in the
macmini evidence is the machine's canonical hostname, not a second host.

## Numerical And Scope Review

The shared sector kernel implements `2J(D)-K(D)` for RHF and
`J(D_alpha+D_beta)-K(D_sigma)` for UHF. The nine sector calls use distinct
`W/Wswap`, `VLR/VRL`, raw center tensors, and conventional retained
`Vijkl` orientations. The committed oracle activates one unsymmetrized raw
sector at a time and separately recovers direct and exchange contractions.
Compressed fixed/ragged windows and the retained verifier provide independent
transformed-basis parity.

The old global cached RHF/UHF paths and their `S`, `Srho`, `rho_slice`, `J`,
and temporary matrices are gone. Cached source shrank from 886 to 422 lines.
M3 added no source file, public type, export, backend API operation, or generic
cache framework.

The exact macmini speedups and cache timings remain engineering provenance.
Do not describe them as complete-solver acceleration, Be timing, peak memory,
or scientific scaling. Milestone 5 must supply those claims.

## Milestone 4 Budget Amendment

M3 used all but three additions in the previous combined M3/M4 gross-addition
cap while deleting 829 lines and leaving the program 279 lines smaller. Three
gross additions cannot support a properly tested backend-general M4. Forcing
the implementation into that number would reward artificial compression or
missing tests.

This review therefore supersedes only the combined three-line remainder:

- restore the original M4 individual cap of at most 90 counted additions;
- retain the M4 deletion target of 15 lines;
- require the final production `src/core.jl` line count not to increase;
- require the overall performance program to remain a substantial net
  shrinkage after M4; and
- continue excluding replacement-only status and this governance memo from
  implementation accounting.

No other M3 accounting or anti-bloat rule changes. Do not fold unrelated dead
code or style cleanup into M4 merely to manufacture deletions.

## Milestone 4 Authority

Execute M4 only: reuse the post-density Fock matrix as the next local
microiteration's pre-diagonalization Fock matrix.

Required behavior and evidence:

- `n+1` backend Fock calls for `n` local microiterations, including exactly
  five calls for four iterations rather than eight;
- common-H RHF/UHF and split-one-body UHF coverage;
- density-density and cached-sliced energy, trajectory, and occupied-projector
  parity against the pre-M4 behavior;
- unchanged density damping, local convergence, energy evaluation, and
  observer behavior;
- exact-zero or side-effect-free backend semantics must not be assumed beyond
  the documented additive Fock contract;
- no new public keyword, type, source file, backend API method, cache layer, or
  compatibility wrapper; and
- a compact timing check showing the expected 37.5% reduction in backend Fock
  calls/work for four forced iterations. Do not require or claim a comparable
  complete-sweep speedup now that the aligned sliced kernel is cheap.

Commit this review before source edits. Stop after the M4 implementation and
full validation for paper-manager review. Do not begin M5 in the same pass.

-- hfdmrg-paper-manager@faraday30.ps.uci.edu
