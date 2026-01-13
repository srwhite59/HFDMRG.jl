# HFDMRG

HFDMRG is a Hartree-Fock solver that uses DMRG-style sweeps. The algorithm
builds a local L-C-R superblock basis from left/right environment blocks and
bare center sites, runs a few SCF micro-iterations (build rho, build F, solve,
occupy, mix), then truncates the block bases via SVD of the occupied orbitals.
Sweeps proceed left-to-right and right-to-left until the energy converges.

The core is interaction-agnostic: it owns the sweep schedule, basis transforms,
and SCF loop, while backends manage interaction caches and add mean-field
contributions to the Fock matrix.

Entry points:
- `solve_hfdmrg(H, V, psiup0; ...)` for density-density RHF.
- `solve_hfdmrg(H, V, psiup0, psidn0; ...)` for density-density UHF.
- `solve_hfdmrg(H, backend::SlicedBasisBackend, psiup0; ...)` for sliced-basis RHF
  (projection-based).
You can also pass `layout::SliceLayout` and `V6` directly; the backend is constructed
for you.

## Repository layout
- `src/core.jl`: generic HF-DMRG sweep engine.
- `src/backend_api.jl`: backend API contract for interactions.
- `src/backends/density_density.jl`: density-density backend for `V::Matrix`.
- `src/backends/sliced_basis.jl`: sliced-basis backend (projection-based).
- `src/slice_layout.jl`: slice layout helper for fixed-size slices.
- `reference/HF_dmrg_legacy.jl`: legacy solver used for regression tests.

## How to run tests
```
~/codexhome/cjulia -e 'using Pkg; Pkg.test()'
```
For quick performance sanity, see `scripts/bench_sliced.jl`.

## How to add a backend
Checklist:
1. Define a backend type and any block/window state you need.
2. Implement `vee_init_block`, `vee_absorb_block`, `vee_window`,
   `vee_add_fock_r!`, and `vee_add_fock!`.
3. Add a regression test that compares against a known solution.

The density-density backend is implemented today; a sliced-basis backend is
implemented via a projection-based window path (correctness-first, slow).
