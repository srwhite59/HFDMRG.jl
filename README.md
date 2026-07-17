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
- `solve_hfdmrg(Hup, Hdn, V, psiup0, psidn0; ...)` for UHF with
  spin-dependent one-body terms and density-density interactions.
- `solve_hfdmrg(H, layout::SliceLayout, V6, psiup0; ...)` for sliced-basis RHF
  with uniform slice sizes.
- `solve_hfdmrg(H, layout::SliceLayout, Vblocks, psiup0; ...)` for sliced-basis
  RHF with ragged slice sizes.
- `solve_hfdmrg(H, backend::SlicedBasisBackend, psiup0; ...)` for sliced-basis RHF
  (projection-based).
The corresponding sliced UHF forms take `psiup0, psidn0`, and split one-body
UHF forms take `Hup, Hdn` before the interaction/backend arguments.

## Usage

All solver entry points return `(psiup, psidn, energy)`, where `psiup`/`psidn` are
the optimized orbitals (orthonormal columns) and `energy` is the final HF energy.

### Density-density interaction (V::Matrix)

This backend expects an N×N interaction matrix `V`, interpreted as a
density-density interaction (see `src/backends/density_density.jl`).

```julia
using LinearAlgebra
using Random
using HFDMRG

rng = MersenneTwister(0)
N = 16
H = randn(rng, N, N); H = (H + H') / 2
V = randn(rng, N, N); V = (V + V') / 2

Nup = 4
psiup0 = Matrix(qr(randn(rng, N, Nup)).Q)

psiup, psidn, energy = solve_hfdmrg(H, V, psiup0; maxiter = 10, blocksize = 4)
```

UHF uses separate initial orbitals:

```julia
Ndn = 3
psidn0 = Matrix(qr(randn(rng, N, Ndn)).Q)

psiup, psidn, energy = solve_hfdmrg(H, V, psiup0, psidn0; maxiter = 10, blocksize = 4)
```

Spin-dependent one-body UHF passes separate up/down one-body Hamiltonians while
leaving the two-body interaction unchanged:

```julia
Hup = H + Diagonal(field)
Hdn = H - Diagonal(field)

psiup, psidn, energy = solve_hfdmrg(Hup, Hdn, V, psiup0, psidn0;
    maxiter = 10, blocksize = 4)
```

The existing `solve_hfdmrg(H, V, psiup0, psidn0; ...)` path remains the
common-H fast path with one projected one-body cache. Passing the same matrix
object as both `Hup` and `Hdn` also routes to that path; only genuinely split
one-body operators allocate split one-body caches.

### Density-density plus a target-space residual

A small signed four-index correction can be added in a fixed orthonormal
target space without constructing a global four-index interaction:

```julia
backend = HFDMRG.DensityDensityTargetResidualBackend(V, Q, residual_pair)
psiup, psidn, energy = solve_hfdmrg(H, backend, psiup0, psidn0;
    maxiter = 10, blocksize = 16)
```

`Q` is `N x m`, with orthonormal target orbitals as columns.
`residual_pair` is a symmetric `P x P` matrix, `P=m(m+1)/2`, storing the
signed residual `(pq|rs)` in canonical order `p=1:m, q=1:p`. Pair coordinates
are not `sqrt(2)` normalized. The corresponding RHF and split-one-body UHF
forms use the same backend argument position as the sliced backends. An
exactly zero residual routes through the original density-density solver.

### External sliced Hamiltonians

For block- or slice-structured Hamiltonians from another package, such as an
atomic sliced export, the normal HFDMRG route is the sliced solver:

```julia
using HFDMRG

# Produced by the consumer package in HFDMRG's slice-ordered orbital basis:
# H::AbstractMatrix          dense one-body term
# dims::Vector{Int}          orbitals per slice
# Vblocks                    ragged sliced two-body tensor
# psiup0                     initial occupied orbitals in the full basis
layout = HFDMRG.SliceLayout(dims)

psiup, psidn, energy = HFDMRG.solve_hfdmrg(H, layout, Vblocks, psiup0;
    maxiter = 10, blocksize = 16)
```

This path does not densify the two-body interaction. The public solver still
accepts the one-body term as a dense global matrix `H`, while the two-body term
stays in the native sliced representation (`V6` or `Vblocks`). There is not a
separate high-level conserved-`m` Hamiltonian entry point; map that structure
onto `H + SliceLayout + V6/Vblocks` before calling `solve_hfdmrg`.

When all slice sizes are equal, prefer the fixed-size `V6` representation over
ragged `Vblocks`:

```julia
layout = HFDMRG.SliceLayout(fill(nj, ns))
psiup, psidn, energy = HFDMRG.solve_hfdmrg(H, layout, V6, psiup0;
    maxiter = 10, blocksize = 16)
```

`V6` is the faster sliced format for uniform slices and should usually be the
first choice in that case.

For spin-dependent one-body terms with sliced interactions, put `Hup, Hdn`
before the sliced layout or backend:

```julia
psiup, psidn, energy = HFDMRG.solve_hfdmrg(Hup, Hdn, layout, V6,
    psiup0, psidn0; maxiter = 10, blocksize = 16)
```

### Sliced-basis interaction (fixed per-slice size)

For fixed slice size `nj`, the sliced interaction is stored as
`V6[a,b,c,d,n,m]` representing `(p q| r s)` with
`p=(n,a), q=(m,c), r=(n,b), s=(m,d)`.

```julia
using LinearAlgebra
using Random
using HFDMRG

rng = MersenneTwister(1)
ns = 4
nj = 2
layout = HFDMRG.SliceLayout(fill(nj, ns))
N = layout.offs[end]

H = randn(rng, N, N); H = (H + H') / 2
V6 = randn(rng, nj, nj, nj, nj, ns, ns)

Nup = 2
psiup0 = Matrix(qr(randn(rng, N, Nup)).Q)

psiup, psidn, energy = HFDMRG.solve_hfdmrg(H, layout, V6, psiup0; maxiter = 5, blocksize = 2)
```

### Sliced-basis interaction (ragged per-slice sizes)

For variable slice sizes, store a ragged 4-tensor for each slice pair:
`Vblocks[n][m]` has size `dims[n] × dims[n] × dims[m] × dims[m]` and uses the
same `(p q| r s)` mapping as the fixed case.

```julia
using LinearAlgebra
using Random
using HFDMRG

rng = MersenneTwister(2)
dims = [1, 2, 3, 2, 2]
layout = HFDMRG.SliceLayout(dims)
ns = length(dims)
N = layout.offs[end]

H = randn(rng, N, N); H = (H + H') / 2
Vblocks = [[randn(rng, dims[n], dims[n], dims[m], dims[m]) for m in 1:ns]
           for n in 1:ns]

Nup = 2
psiup0 = Matrix(qr(randn(rng, N, Nup)).Q)

psiup, psidn, energy = HFDMRG.solve_hfdmrg(H, layout, Vblocks, psiup0; maxiter = 5, blocksize = 2)
```

### Cached sliced backend

The cached sliced backend uses precomputed block tensors for speed. It shares
the same interaction conventions as the projection backend and can be constructed
directly from `layout` and `V6`/`Vblocks`:

```julia
backend = HFDMRG.SlicedBasisBackendCached(layout, V6)
psiup, psidn, energy = HFDMRG.solve_hfdmrg(H, backend, psiup0; maxiter = 5, blocksize = 2)
```

For ragged interactions, pass `Vblocks` instead of `V6`.

## Interaction representations

Supported representations:
- Density-density: `V::Matrix` (N×N).
- Fixed sliced: `V6::Array{Float64,6}` with shape `(nj, nj, nj, nj, ns, ns)`.
- Ragged sliced: `Vblocks::Vector{Vector{Array{Float64,4}}}` with per-slice sizes.

Backend constructors:
- `HFDMRG.SlicedBasisBackend(layout, V6)` or
  `HFDMRG.SlicedBasisBackend(layout, Vblocks)`
- `HFDMRG.SlicedBasisBackendCached(layout, V6)` or
  `HFDMRG.SlicedBasisBackendCached(layout, Vblocks)`

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

## Developer commands
- `make test`: run the package tests via `cjulia`.
- `make bench`: run `scripts/bench_sliced.jl`.
- `make verify`: run `scripts/verify_cached_vs_projection.jl`.

## How to add a backend
Checklist:
1. Define a backend type and any block/window state you need.
2. Implement `vee_init_block`, `vee_absorb_block`, `vee_window`,
   `vee_add_fock_r!`, and `vee_add_fock!`.
3. Add a regression test that compares against a known solution.

The density-density backend is implemented today; a sliced-basis backend is
implemented via a projection-based window path (correctness-first, slow).
