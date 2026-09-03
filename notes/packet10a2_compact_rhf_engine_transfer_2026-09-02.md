# Packet 10A2: compact RHF engine transfer candidate

> Superseded in architecture by Packet 10A2R2. Packet 10A2R2 retains this
> compact engine unchanged but restores the scientifically distinct historical
> dense/sliced RHF engine behind disjoint `solve_hfdmrg` signatures. The
> migration boundaries and single-engine line target below describe the
> accepted intermediate commit, not the corrected public contract.

## Outcome

Packet 10A2 installs the standalone terminal-complete compact RHF lifecycle
and monotone solver behind the sole public solver name, `solve_hfdmrg`.  The
candidate is deliberately unstaged and uncommitted for review on
`packet10a1-compact-rhf-foundation`, whose parent is Packet 10A1 commit
`afbe940d0f785835be93c8c086c6dad491fb256e`.

The compact route owns one exact-sized collection and one live RHF root.  It
uses block-native dimer initialization, exact-sized terminal fragments for
remainders 0--17, physical/reflected traversal, live-root turnaround, exact or
finite state selection, one-Fock proposals, bounded largest-safe geodesic
fallback, invariant represented-energy convergence, and update-disabled replay
diagnostics.  It does not own a global coefficient panel or an `N x N` array.

The historical UHF engine and its still-required sliced/cached backend
utilities remain unchanged in this transition packet.  Consequently the
temporarily loaded source tree is 8,751 lines, while the standalone compact RHF
production slice is 2,793 lines.  Packet 10A3 transfers independent-spin UHF;
the accepted projection remains approximately 5,000--5,400 lines after the
old transitional engine is removed in the bounded later steps.

## Public dispatch and compatibility boundaries

There are 19 `solve_hfdmrg` methods: one compact RHF method, six deterministic
RHF migration boundaries, and twelve unchanged historical UHF methods.  Every
RHF signature that formerly entered the global-coefficient engine now either
runs the compact unit-cell method or throws before allocation.  The private
history-accelerated RHF driver is also fail-closed, while its common policy,
Pulay, model-audit, and UHF history utilities remain loaded for the unchanged
UHF path.

| Historical route/signature | Former behavior | Packet 9A4B disposition and Packet 10A2 test treatment |
|---|---|---|
| `solve_hfdmrg(H, V, psiup0; ...)` | constructed `DensityDensityBackend` and entered restricted global-coefficient sweeps | superseded dense RHF route; replaced by compact facade tests plus deterministic migration-boundary tests |
| `solve_hfdmrg(H, backend::SlicedBasisBackend, psiup0; ...)` | entered restricted global-coefficient sweeps using projected sliced contractions | superseded sliced RHF route; solver boundary retained, numerical backend/Fock tests retained independently |
| `solve_hfdmrg(H, layout::SliceLayout, V, psiup0; ...)` | constructed a sliced backend and forwarded to the sliced RHF solver | superseded sliced-layout route; compatibility boundary retained, layout/backend utilities retained |
| `solve_hfdmrg(H, backend::SlicedBasisBackendCached, psiup0; ...)` | entered restricted global-coefficient sweeps with cached environments | superseded cached-sliced RHF route; solver boundary retained, cache seed/growth/window tests retained |
| `solve_hfdmrg(H, backend::_SlicedBasisBackendCachedCoulomb, psiup0; ...)` | entered restricted global-coefficient sweeps with packed Coulomb caches | superseded cached-Coulomb RHF route; solver boundary retained, packed backend arithmetic tests retained |
| `solve_hfdmrg(H, backend::DensityDensityTargetResidualBackend, psiup0; ...)` | forwarded zero residual to dense RHF or entered restricted target-residual sweeps | superseded target-residual RHF route; solver boundary retained, target-space Fock/action tests retained |

No boundary silently redirects legacy `maxiter`, `cutoff`, `scf_cutoff`, or
`environment_cutoff` to the compact state and convergence controls.  The error
names `BandedOneBody`, `UnitCellInteraction`, `state_maxdim`, `state_cutoff`,
`energy_tolerance`, and `maximum_half_sweeps` explicitly.

The compact RHF result is `HFDMRGResult`, containing convergence status and
reason, represented energy, half-sweep and update counts, rank/discarded-weight
diagnostics, replay diagnostics, and a compact state handle.  The former
three-matrix tuple is not emulated because doing so would materialize the
global occupied coefficients.  Historical UHF calls retain their existing
tuple result in Packet 10A2.  Thus RHF compatibility is an explicit migration
boundary, not a numerically different silent adapter.

## Exact file disposition in this step

New compact production owners are `src/fixed_state_lifecycle.jl` and
`src/nonlinear_rhf.jl`.  `src/HFDMRG.jl` is adapted to load/export the compact
route and describe the transitional dispatch accurately.  `src/core.jl`,
`src/backends/sliced_basis.jl`, `src/backends/sliced_basis_cached.jl`,
`src/backends/sliced_basis_cached_coulomb.jl`, and
`src/backends/density_density_target_residual.jl` retain their UHF and backend
arithmetic but disconnect their restricted solver overloads.
`src/history_accelerated_rhf.jl` retains common utilities required by UHF while
its private RHF driver fails closed.

No tracked historical source file is removed in this step because the
unchanged historical UHF engine still consumes the old core/backend owners.
The tests added are `test/compact_rhf_lifecycle.jl`,
`test/compact_terminal_fragments.jl`, `test/compact_nonlinear_rhf.jl`, and
`test/compact_rhf_compatibility.jl`; `test/runtests.jl` is adapted to replace
only superseded RHF solver coverage.  All other historical source and test
files are retained unchanged.  The new validation owner and evidence record
are under `validation/packet10a2/`.

## Test disposition

The historical baseline executed 816 assertions.  Packet 10A2 retains 702
unchanged historical assertions, repurposes two existing private-RHF-history
assertions as deterministic boundary checks, and retires 112 assertions whose
sole contract was a superseded RHF solver route.  It adds 607 compact RHF
lifecycle/fragment/nonlinear assertions and 18 public compatibility-boundary
assertions.  Together with the 41 Packet 10A1 compact-foundation assertions,
each complete suite executes 1,370 assertions.

No test is conditionally skipped.  Retired RHF assertions were removed from
mixed test sets without changing their UHF or backend-utility assertions.  In
particular, sliced Fock/window arithmetic, cached absorption, cached Coulomb,
target-residual arithmetic, timing, and UHF history tests remain active.

## Validation

Julia 1.12.6 with caller-selected four-thread BLAS passed:

- all 1,370 ordinary package assertions;
- all 1,370 assertions with `--check-bounds=yes`;
- all 1,370 assertions in a fresh `/tmp` checkout instantiated without the
  repository Manifest;
- the six-route public compatibility boundary (18 assertions);
- export, method-table, sole-dependency, caller-BLAS, and zero-allocation
  warmed-preflight checks in `validation/packet10a2/check_transfer.jl`;
- `git diff --check`.

After the complete checked-bounds and isolated runs, the final source delta was
limited to truthful module wording plus the private history-RHF fail-closed
boundary.  The final ordinary 1,370-assertion suite was rerun, and the public
six-route plus private-history boundary was rerun under checked bounds in the
updated isolated checkout.

The compact H2 audit converged in five half-sweeps to represented energy
`-0.44914425911533473` Ha.  Warmed nonlinear preflight allocated zero bytes and
left the four-thread BLAS setting unchanged.  Existing permanent compact tests
also cover zero-allocation lifecycle/center preparation, publication, root/link
operations, failure atomicity, terminal fragments, exact/finite controls, and
physical/reflected replay.  Historical `TimeG` identity/count coverage remains
unchanged and passed as part of the 124 timing assertions.

No H1000 run was performed; the accepted five-step manifest reserves bounded
H1000 RHF/UHF extraction equivalence for Packet 10A5.

## Repository state

Historical `main`, remotes, PPP, and the eight pre-existing untracked files
remain unchanged.  Their SHA-256 values are recorded in
`validation/packet10a2/packet10a2_evidence.toml`.  Packet 10A2 is unstaged and
uncommitted.  The tracked historical Manifest is unchanged; its removal remains
assigned to the later manifest/dead-file audit rather than this transition.

The checked-in README and manual still describe the historical public surface;
the accepted transfer manifest assigns their complete facade update to Packet
10A4.  The module docstring now states the transitional RHF/UHF dispatch
truthfully so method help does not advertise the removed RHF routes.
