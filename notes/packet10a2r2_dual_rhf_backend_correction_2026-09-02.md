# Packet 10A2R2: dual RHF backend correction

Date: 2026-09-02

## Outcome

Packet 10A2R2 passes its complete gate and is accepted for a new linear
corrective commit based on accepted Packet 10A2 commit
`18afb6f476c1817727e59f636e4f33fb912d10f1`. The correction was built and
validated in an isolated machine-local worktree; the primary transfer
worktree and its provisional Packet 10A3 files were not used.

HFDMRG now preserves two scientifically distinct implementations behind the
one exported `solve_hfdmrg` name:

- typed `BandedOneBody, UnitCellInteraction` calls enter the compact
  block-native unit-cell RHF solver and return `HFDMRGResult`;
- matrix, sliced-layout, sliced-backend, cached-backend, cached-Coulomb, and
  target-residual calls enter the historical dense/sliced solver and return
  `(psiup, psidn, energy)` with global occupied coefficients.

There is no fallback between the implementations. The public method table has
19 methods and no ambiguity: one compact RHF method, six historical RHF
methods, and twelve unchanged historical UHF methods.

## Restored historical capability

The correction restores the six RHF routes retired by Packet 10A2:

| Route | Dispatch and result |
|---|---|
| dense interaction | `solve_hfdmrg(H, V, psiup0; ...)`; historical tuple |
| sliced backend | `solve_hfdmrg(H, SlicedBasisBackend, psiup0; ...)`; historical tuple |
| sliced layout | `solve_hfdmrg(H, SliceLayout, V, psiup0; ...)`; historical tuple |
| cached sliced backend | `solve_hfdmrg(H, SlicedBasisBackendCached, psiup0; ...)`; historical tuple |
| cached Coulomb backend | `solve_hfdmrg(H, _SlicedBasisBackendCachedCoulomb, psiup0; ...)`; historical tuple |
| target-residual backend | `solve_hfdmrg(H, DensityDensityTargetResidualBackend, psiup0; ...)`; historical tuple |

The historical private RHF-history driver is restored. Dense initial occupied
coefficients, global coefficient output, legacy controls, and corresponding
history behavior retain their historical semantics. The twelve UHF public
methods and their numerical behavior were not changed.

The only method-table correction beyond restoration is a truthful type bound
on the sliced-layout interaction argument. It accepts exactly the three
representations already accepted by the constructor: fixed
`Array{Float64,6}`, `SlicedVeeRagged`, or the ragged vector-of-vector
four-tensor form. This removes intersections with split-one-body backend
signatures for nonsensical mixed argument lists; it does not remove a
constructible sliced input.

## Compact route preservation

All eight compact production files are byte-identical to commit `18afb6f`:

- `src/unit_cell_interaction.jl`
- `src/pairs.jl`
- `src/state.jl`
- `src/one_body.jl`
- `src/interaction_blocks.jl`
- `src/interaction_centers.jl`
- `src/fixed_state_lifecycle.jl`
- `src/nonlinear_rhf.jl`

The compact route retains one exact-sized collection and live root,
block-native dimer initialization, terminal-complete centers and fragments,
monotone one-Fock updates, finite selection/replay semantics, zero later
coefficient-row reads, and no global coefficient panel, `N x N` owner,
second collection, reflected operator copy, or PPP dependency.

## Production ownership and size

Physical source-line accounting is:

| Category | Lines |
|---|---:|
| compact unit-cell RHF engine | 2,793 |
| historical dense/sliced engine and history | 5,509 |
| shared module/timing interface | 279 |
| installed optional target-residual extension | 183 |
| total loaded production | 8,764 |

The former 5,160-line single-engine target is superseded. The physical total is
slightly above the earlier 7,000--8,000 estimate because capability
preservation retains all three historical history owners and the cached,
cached-Coulomb, target-residual, frozen-occupied, and sliced-audit support.
There is still one package and one exported solver name, but two deliberately
different numerical engines and two honest result contracts.

Future bond-parameterized or cut-dependent compact interactions, compact
pair-response compression, compact history, and disk-backed blocks remain
deferred. They were not implemented or counted as present compact features.

## Test restoration and validation

The original historical suite has all 816 assertions restored. This includes
the 112 RHF assertions retired by Packet 10A2 and all 39 private
history-accelerated RHF assertions. No test is skipped. The complete suite
also runs 648 compact RHF assertions and 47 dual-dispatch assertions, for
1,511 assertions total.

Julia 1.12.6 with the caller's four BLAS threads passed:

- the complete ordinary 1,511-assertion suite;
- the complete `Pkg.test()` checked-bounds 1,511-assertion suite;
- a fresh Manifest-free isolated package copy with all 1,511 assertions;
- all historical UHF, backend, utility, timing, allocation, and history tests;
- all compact RHF foundation, lifecycle, terminal-fragment, nonlinear,
  allocation, and atomicity tests;
- six executable historical RHF dispatch boundaries and their corresponding
  historical UHF availability checks;
- exports, 19-method inventory, zero ambiguities, sole runtime dependency,
  caller-BLAS preservation, and zero-byte warmed compact preflight;
- `git diff --check`.

A Julia 1.10 exploratory run was not counted: it used the wrong system Julia
and differed on one historical damping-path assertion. No code was changed in
response; the supported Julia 1.12 runs all passed.

## GaussletBases handshake

The unchanged GaussletBases production operation
`run_atomic_injected_angular_hfdmrg_hf` was exercised read-only with its
restricted three-matrix payload. GaussletBases was the explicit active project,
and HFDMRG was loaded from the isolated Packet 10A2R2 source path rather than
ambient discovery or the dirty primary checkout.

The constructed angular adapter had basis dimension 840. Its unchanged call
returned global alpha/beta coefficient matrices and finite energy
`13.02625601311897` Ha. The restricted alpha and beta outputs were equal.
GaussletBases source, tests, dependency declarations, and working tree were not
modified.

## Repository state and non-actions

The correction is based directly on parent `18afb6f`. Historical main,
remotes, dependency declarations, PPP, GaussletBases, and the eight
pre-existing untracked historical files are unchanged. Packet 10A3 remains
paused in the primary worktree with its provisional files preserved by
recorded SHA-256 fingerprints; none participate in this correction or its
validation.

Packet 10A3 was not resumed. Historical UHF was not altered, no chemistry
extension was added, and no dependency declaration was changed. The accepted
landing is a new corrective commit after `18afb6f`; existing transfer history
is not rewritten.
