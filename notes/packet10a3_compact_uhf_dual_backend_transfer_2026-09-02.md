# Packet 10A3: compact UHF with preserved historical UHF

Date: 2026-09-02

## Outcome

Packet 10A3 passes its complete implementation and validation gate as an
unstaged, uncommitted review candidate based directly on Packet 10A2R2 commit
`19e11fc2b12a142138dc3324f585d8aa92d46098`.

HFDMRG now retains two scientifically distinct UHF implementations behind the
one exported `solve_hfdmrg` name:

- typed `BandedOneBody, UnitCellInteraction` calls with `spin=:uhf` enter the
  compact independent-spin unit-cell engine and return `UHFDMRGResult`;
- all twelve matrix, sliced, cached, cached-Coulomb, and target-residual UHF
  signatures continue to enter the historical dense/sliced engine and return
  `(psiup, psidn, energy)` with global occupied coefficients.

There is no fallback between the implementations. The complete public method
table remains 19 methods with no ambiguity: one typed compact RHF/UHF method,
six historical RHF methods, and twelve historical UHF methods.

## Compact UHF transfer

The compact addition installs:

- independently sized alpha and beta roots and reversible state links;
- dual projected-H1 records and independent packed exchange-pair payloads;
- one shared Hartree connector and a distinct factored internal cross-spin
  Hartree record;
- a spin-generic two-sided center with a thin joint UHF transaction;
- block-native atomic-Neel initialization, exact terminal fragments, one
  collection and one root object, and live-root turnaround;
- simultaneous one-Fock proposals, safe full-step acceptance, one-common-step
  geodesic fallback, atomic rejection, invariant represented-energy
  convergence, and exact/finite replay semantics.

The four compact UHF production files and the two substantive focused test
files are byte-identical to the frozen trimmed/HFDMRG donor. Only narrow module
includes/exports, the existing typed facade, dispatch-isolation tests, and the
test include list were adapted in the recipient.

The compact path owns no common-union pair field, persistent common-to-spin
maps, dense cross-spin pair map, global alpha/beta coefficient panel, `N x N`
owner, second collection, reflected operator copy, PPP callback, or historical
backend callback.

## Historical capability preservation

Packet 10A3 does not remove or redirect historical UHF. All twelve historical
methods retain arbitrary dense and sliced inputs, dense initial occupied
coefficients, global tuple output, cached backends, split-one-body support, and
historical convergence/history machinery. The six restored historical RHF
methods from Packet 10A2R2 also remain unchanged.

The unchanged GaussletBases production adapter was exercised with explicit
candidate source loading and GaussletBases as the active project. Two bounded
840-dimensional matrix-dispatch smoke calls passed:

- restricted three-matrix RHF: `13.02625601311897` Ha;
- one-electron open-shell four-matrix UHF: `4.916191887700468` Ha.

Both smoke calls returned the historical global coefficient tuple contract.
The old Pass 583 value is a separate three-electron handshake with incomplete
retained controls and is not the governing preservation certificate.

The governing Packet 10A3R0 preservation certificate is the matched
clean-versus-candidate two-alpha, one-beta comparison at
`39.38381696539275` Ha. It has zero energy and projector differences and
identical alpha/beta coefficient and projector hashes. Neither GaussletBases
nor its dependency declarations were modified.

## Numerical and ownership evidence

The compact tests cover H2/H10/H20 exact and finite UHF, both starting
orientations, spin swap, nonorthogonal and unequal spin spaces, zero-rank spin
sectors, center energy/Fock/vector/panel actions, terminal fragments,
four-half-sweep fixed-state traversal, monotone full steps, common-geodesic
fallback, forced rejection, update-disabled replay, failure atomicity,
hostile capacity, allocation, timing, and BLAS preservation.

The complete Julia 1.12.6 results are:

- ordinary direct suite: 1,938/1,938 assertions;
- checked-bounds `Pkg.test()`: 1,938/1,938 assertions;
- fresh Manifest-free isolated `Pkg.test()`: 1,938/1,938 assertions;
- `solve_hfdmrg` methods: 19; ambiguities: 0;
- warmed compact-UHF preflight allocation: 0 bytes;
- caller BLAS threads: 4 before and after compact tests;
- compact-to-historical numerical callbacks: 0.

An initial exploratory command used the system Julia 1.10 executable and
reproduced Packet 10A2R2's already documented historical damping assertion
difference. It is excluded from the gate. No source was changed in response;
all governing Julia 1.12.6 suites passed.

No H1000 calculation was run. Producer evaluation, optional compact operator
representations, pair-response compression, compact history, and disk storage
remain outside this packet.

## Production size

Physical loaded source-line accounting is:

| Category | Lines |
|---|---:|
| compact common and RHF | 2,799 |
| compact independent-spin UHF addition | 1,785 |
| historical dense/sliced engines and history | 5,509 |
| shared module and timing utilities | 287 |
| optional installed historical target-residual feature | 183 |
| total loaded production | 10,563 |

The increase from Packet 10A2R2 is 1,799 production lines: 1,785 lines of
compact UHF numerical ownership plus fourteen lines of module/facade
integration. Capability preservation takes precedence over speculative
cross-model consolidation.

## Files and repository state

Production additions are `src/uhf_blocks.jl`, `src/uhf_centers.jl`,
`src/uhf_lifecycle.jl`, and `src/nonlinear_uhf.jl`. Narrow adaptations are
limited to `src/HFDMRG.jl` and `src/nonlinear_rhf.jl`.

Focused tests are `test/compact_uhf_foundation.jl`,
`test/compact_nonlinear_uhf.jl`, and `test/compact_timing.jl`, with include and
dispatch-boundary adaptations in `test/runtests.jl` and
`test/compact_rhf_compatibility.jl`.

The candidate is on branch `packet10a3-compact-uhf-dual-backend` in an
isolated machine-local worktree. It is unstaged and uncommitted. Historical
main, remotes, PPP, GaussletBases, the accepted GaussletBases relocation, and
the eight pre-existing historical untracked files are unchanged. The earlier
provisional Packet 10A3 work remains separately preserved in the primary
worktree with its recorded fingerprints unchanged.
