# HFDMRG

HFDMRG is a Hartree-Fock solver that uses DMRG-style sweeps. The algorithm
builds a local L-C-R superblock basis from left/right environment blocks and
bare center sites, runs a few SCF micro-iterations (build rho, build F, solve,
occupy, mix), then truncates the block bases via SVD of the occupied orbitals.
Sweeps proceed left-to-right and right-to-left until the energy converges.

The core is interaction-agnostic: it owns the sweep schedule, basis transforms,
and SCF loop, while backends manage interaction caches and add mean-field
contributions to the Fock matrix.

## Repository layout
- `src/core.jl`: generic HF-DMRG sweep engine.
- `src/backend_api.jl`: backend API contract for interactions.
- `src/backends/density_density.jl`: density-density backend for `V::Matrix`.
- `reference/HF_dmrg_legacy.jl`: legacy solver used for regression tests.

## How to run tests
```
~/codexhome/cjulia -e 'using Pkg; Pkg.test()'
```

## How to add a backend
Checklist:
1. Define a backend type and any block/window state you need.
2. Implement `vee_init_block`, `vee_absorb_block`, `vee_window`,
   `vee_add_fock_r!`, and `vee_add_fock!`.
3. Add a regression test that compares against a known solution.

The density-density backend is implemented today; a sliced-basis backend is
future work.
