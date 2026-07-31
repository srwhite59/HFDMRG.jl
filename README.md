# HFDMRG

HFDMRG solves restricted or unrestricted Hartree-Fock problems using
DMRG-style sweeps over an ordered one-particle basis. The name describes how
the calculation is organized: HFDMRG still returns a single Slater
determinant, not a correlated matrix-product state.

The core is interaction-agnostic: it owns the sweep schedule, basis transforms,
and SCF loop, while backends manage interaction caches and add mean-field
contributions to the Fock matrix.

## Installation

HFDMRG is currently tested with Julia 1.12. From a source checkout, instantiate
the project and verify the installation with:

```sh
cd /path/to/hfdmrg
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using HFDMRG; println("HFDMRG loaded")'
```

To use the checkout from another Julia project, replace the paths below and
develop HFDMRG into that environment:

```sh
julia --project=/path/to/my-project -e \
  'using Pkg; Pkg.develop(path="/path/to/hfdmrg")'
```

## Quick start

This deterministic restricted-HF example builds a small one-dimensional model,
runs up to four complete HFDMRG sweeps, and checks the returned orbitals. It is
also available as [`examples/quickstart.jl`](examples/quickstart.jl).

```julia
using HFDMRG
using LinearAlgebra
using Random

rng = MersenneTwister(7)
N = 12
Nocc = 2

onsite = collect(range(-1.0, 1.0; length = N))
hopping = diagm(1 => fill(-0.1, N - 1), -1 => fill(-0.1, N - 1))
H = Matrix(Diagonal(onsite) + hopping)
V = [0.4 / (1 + abs(i - j)) for i = 1:N, j = 1:N]
psi0 = Matrix(qr(randn(rng, N, Nocc)).Q)[:, 1:Nocc]

psiup, psidn, energy = solve_hfdmrg(H, V, psi0;
    maxiter = 4, blocksize = 2, cutoff = 1e-9,
    scf_cutoff = 1e-9, verbose = false)

@assert psiup == psidn
@assert norm(psiup' * psiup - I) < 1e-10
println("energy = ", round(energy; digits = 12))
```

Run the checked-in example from the repository root with:

```sh
julia --project=. examples/quickstart.jl
```

Representative output on Julia 1.12 is:

```text
energy = -2.41361165233
```

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

## How HFDMRG works

An ordinary Hartree-Fock iteration updates all occupied orbitals in one global
step. HFDMRG instead moves a smaller working problem through the ordered
physical basis. At each position it divides that basis into

```text
retained left environment | explicit center sites | retained right environment
```

The center contains bare physical sites. The permanently fixed first and last
physical blocks use complete coordinate-identity bases, so startup cannot
discard unoccupied boundary directions. Interior environments grown during
initialization or a sweep are represented by smaller bases obtained from the
occupied orbitals on that side. This is the DMRG-like part of the algorithm;
there is no user-selected many-body bond dimension.

One calculation proceeds as follows:

1. **Build the initial environments.** The physical basis is divided into
   contiguous chunks, and the two fixed edge chunks are seeded with literal
   identity bases over their complete physical ranges. `blocksize` controls
   the requested chunk size, while `nblockcenter` controls how many center
   chunks remain explicit.
2. **Assemble one moving window.** The one-body Hamiltonian is transformed into
   the current left-center-right basis. The active backend supplies the needed
   block and window interaction caches from its density-density, sliced, or
   target-residual representation.
3. **Solve the local mean-field problem.** HFDMRG builds the spin densities
   and Fock matrices, diagonalizes the Fock matrices, occupies the lowest
   `Nup` and `Ndn` orbitals, and mixes the updated densities when needed. A
   window receives at most four of these SCF micro-iterations.
4. **Move and rebuild an environment.** An SVD of the occupied-orbital
   coefficients supplies the new environment basis as the center advances.
   `environment_cutoff` controls which singular vectors are retained. The
   solver transforms the one-body and backend caches into that basis; it does
   not carry a tunable DMRG bond dimension.
5. **Complete both directions.** A left-to-right pass followed by a
   right-to-left pass is one sweep. After the sweep, HFDMRG expands the
   orbitals back into the full physical basis, tests the energy change, and
   calls the optional observer.

If a fixed edge chunk has physical dimension `d`, its complete retained rank
is also `d`. The density-density and cached sliced routes therefore retain an
edge interaction term with `O(d^4)` storage. Grown interior environments remain
occupied-SVD compressed, so this cost is localized to the two fixed edges;
unusually large edge chunks can nevertheless be expensive.

This organization is intended for ordered bases whose interactions fit one of
HFDMRG's compact backend representations. It is still Hartree-Fock: the
orbital occupations remain fixed during a solve, different initial guesses
can reach different self-consistent determinants, and energy convergence does
not prove that the global Hartree-Fock minimum was found.

## Usage

All solver entry points return `(psiup, psidn, energy)`, where `psiup`/`psidn` are
the optimized orbitals (orthonormal columns) and `energy` is the final HF energy.

### Solver keywords

All `solve_hfdmrg` entry points accept the same solver keywords:

| Keyword | Default | Meaning |
|---|---:|---|
| `blocksize` | `200` | Requested number of physical sites in each interior block chunk. HFDMRG may reduce it for small systems to obtain more than four blocks. This is not a DMRG bond dimension. |
| `block_partition` | `nothing` | Optional `SliceLayout` defining one complete physical slice per block for any backend. It overrides `blocksize` and requires `nblockcenter >= 1` and enough slices for the sweep; when the backend owns a `SliceLayout`, the two layouts must match exactly. |
| `nblockcenter` | `1` | Number of consecutive physical chunks or slices kept explicit between the left and right environments in each moving window. Larger values enlarge the local variational space while increasing local work. |
| `maxiter` | `1000` | Maximum number of complete left-to-right plus right-to-left sweeps. |
| `cutoff` | `1e-11` | Sweep-energy convergence tolerance. Common-H routes use `abs(Eold - Enew) < cutoff`; genuinely split-H routes use `abs(Eold - Enew) < cutoff * max(1, abs(Enew))`. |
| `scf_cutoff` | `nothing` | Local SCF energy tolerance, using the same absolute common-H or scale-aware split-H comparison as `cutoff`. `nothing` means `cutoff` on common-H routes and `cutoff / 10` on genuinely split-H routes. Each window performs at most four SCF micro-iterations. |
| `environment_cutoff` | `1e-10` | Absolute cutoff for singular values of the restricted occupied-orbital coefficient matrix used to build grown environments. Vectors with `sigma > environment_cutoff` are retained; at least one is always kept. Must be finite and nonnegative. |
| `frozen_occupied` | `:off` | `:roundoff_exact` removes only roundoff-certified completed occupied modes from density-density environments while retaining recoverable fields and energy. |
| `observer` | `nothing` | Callable notified after every complete sweep. Return `true` to stop, or `false`/`nothing` to continue. See below. |
| `verbose` | `false` | Print block decomposition and sweep-energy progress. |

“Common-H” includes RHF, UHF with one shared `H`, and split-form calls where
`Hup === Hdn`; those calls use the common-H fast path. Passing distinct `Hup`
and `Hdn` objects selects the split-H convergence rules, even if their entries
are numerically equal.

F1 `:roundoff_exact` supports density-density RHF/UHF and zero target residuals;
sliced routes reject it, and its dense Aufbau audit is not a long-chain method.

### Convergence and troubleshooting

HFDMRG declares convergence from the end-of-sweep energy change described in
the keyword table. Reaching `maxiter` is not an error: if the threshold has
not been met, the solver returns the orbitals and energy from the last
completed sweep. Use an [observer](#per-sweep-observer) when the caller needs
to distinguish those outcomes programmatically:

```julia
last_info = Ref{Any}(nothing)
observer = info -> (last_info[] = info; false)

psiup, psidn, energy = solve_hfdmrg(H, V, psiup0;
    maxiter = 20, blocksize = 16, observer)

last_info[].converged || @warn "HFDMRG reached maxiter without convergence"
```

`verbose = true` is useful for an interactive run; it prints the block
decomposition and the old and new end-of-sweep energies. An observer is the
better choice for saving a durable energy history or checkpoint.

Tune one control at a time:

- `nblockcenter` is a scientific window-size control. Larger values expose a
  larger local variational space and may reduce sweep count or improve basin
  and stationarity behavior. They slightly reduce the number of windows but
  increase local Fock, eigensolve, and cache-window work. Cached sliced runs
  still require `block_partition` to align with the interaction slices.
- `environment_cutoff` controls environment compression independently of
  `cutoff` and `scf_cutoff`, which govern sweep-energy and local-SCF stopping.
  Lower values retain more restricted occupied-orbital directions and increase retained
  ranks; compare energy and stationarity as well as runtime and memory.
- Interpret `cutoff` in the units of the supplied Hamiltonian. A looser value
  is useful for exploratory runs; tighten it and compare the final orbitals
  and energy for production work.
- Increasing `blocksize` makes each explicit center chunk larger and usually
  makes each window more expensive, but can reduce the number of windows in a
  sweep. Decreasing it has the opposite tradeoff. HFDMRG may reduce a requested
  size on small systems so that it can form enough blocks.
- Increase `maxiter` only when the saved sweep history is still approaching a
  stable value. More sweeps do not repair inconsistent inputs or necessarily
  escape a poor basin of attraction.
- For a consequential result, repeat with a different initial guess and at
  least one reasonable `blocksize`. Hartree-Fock can have multiple stationary
  solutions, so compare their converged energies and scientific diagnostics.

Common symptoms and checks:

| Symptom | What to check |
|---|---|
| Energy oscillates or rises between sweeps | End-of-sweep energies are not guaranteed to be monotone. HFDMRG reduces subsequent density-update weights after a local rise; if oscillation persists, try a better initial determinant and compare a larger `blocksize`. |
| The run reaches `maxiter` | Inspect the observer history before simply extending the run. A decreasing energy change may justify more sweeps; a flat oscillation usually calls for a different initial guess or window size. |
| Different initial guesses give different energies | This is possible for RHF and especially UHF. Check that the requested spin occupations are the same, converge each run to the same tolerance, and compare the resulting self-consistent energies and scientific diagnostics. |
| Returned orbitals are not orthonormal | Check `norm(psiup' * psiup - I)` and the corresponding beta quantity. A large error indicates invalid input, numerical failure, or mutation by an observer; ordinary roundoff should be small. |
| An independently recomputed energy differs from `energy` | First make sure the final sweep converged. During damped, unconverged iterations the internal density can be a mixture of determinants, so exact equality with the returned determinant is not promised. Also include the same interaction convention and exclude any nuclear or consumer offset. |
| A sliced run is slow or uses too much memory | Confirm that uniform slices use `V6`, production runs explicitly construct `SlicedBasisBackendCached`, and the solve passes `block_partition=layout`. See [Choosing an interaction backend](#choosing-an-interaction-backend). Observer work is also part of solver wall time. |

Before tuning, verify the input contract: one-body matrices and the
density-density `V` should be finite and symmetric; every orbital matrix must
have `N` rows and orthonormal columns; UHF column counts set the fixed alpha and
beta occupations. For sliced inputs, `sum(layout.dims)` must equal `N` and the
interaction block sizes must match the layout; the current sliced constructors
require `Float64` interaction arrays. For a target residual, `V`, `Q`, and
`residual_pair` must share one `AbstractFloat` element type, the columns of `Q`
must be orthonormal, and `residual_pair` must be symmetric.

Energy convergence is a stopping test, not a complete proof of
self-consistency or uniqueness. High-consequence calculations should also
check orbital orthogonality, electron counts, stability with respect to the
initial guess and window size, and any problem-specific Fock or residual
diagnostic available to the caller.

### Numerical conventions

HFDMRG uses real occupied-orbital matrices with orthonormal columns and expects
real symmetric one-body matrices; the density-density matrix backend likewise
expects a real symmetric `V`. It performs no unit conversion: the returned
energy has the units of the supplied one- and two-body operators. In RHF,
`psiup0` has one column per occupied spatial orbital and represents the density
of one spin, so `Nocc` columns mean `2Nocc` electrons. UHF supplies the alpha
and beta occupied orbitals separately.

Write `D_alpha` and `D_beta` for the one-spin density matrices used by the
local SCF. Without damping, `D_sigma = C_sigma * C_sigma'`; with damping it can
be a mixture of consecutive determinant densities. The Fock matrices include
their one-body terms. Every interaction backend must follow the core energy
convention

```text
E_UHF = 1/2 Tr[D_alpha (F_alpha + H_alpha)]
      + 1/2 Tr[D_beta  (F_beta  + H_beta)]

E_RHF = Tr[D (F + H)]
```

Thus `energy` is the electronic energy for the Hamiltonian passed to HFDMRG.
It does not include nuclear repulsion or a consumer-defined scalar offset.

For the density-density matrix backend, the implied four-index interaction is
`(ij|kl) = delta(i,j) delta(k,l) V[i,k]`. Define
`q = diag(D_alpha + D_beta)` and let `.*` denote elementwise multiplication.
The unrestricted Fock matrices and energy are

```text
F_sigma = H_sigma + Diagonal(V * q) - V .* D_sigma

E_UHF = Tr[D_alpha H_alpha] + Tr[D_beta H_beta]
      + 1/2 q' V q
      - 1/2 sum_ij V[i,j] (D_alpha[i,j]^2 + D_beta[i,j]^2).
```

For RHF, put `D_alpha = D_beta = D` and `d = diag(D)`:

```text
F_RHF = H + 2 Diagonal(V * d) - V .* D

E_RHF = 2 Tr[D H] + 2 d' V d - sum_ij V[i,j] D[i,j]^2.
```

The local SCF may damp densities when an update raises the energy. Away from
convergence, an energy independently recomputed from the returned determinant
therefore need not agree bit-for-bit with the reported local-SCF energy. Sliced
and target-residual backends use the same core trace convention; the latter's
signed pair-storage formulas are recorded in the
[target-residual design](notes/target_space_interaction_residual_design_2026-07-17.md#fock-and-energy-formulas).

### Per-sweep observer

Every solver entry point accepts an optional callable `observer`. It is called
once after each complete left-to-right plus right-to-left sweep, including the
final converged sweep:

```julia
function observer(info)
    println("sweep=$(info.sweep) energy=$(info.energy)")
    write_checkpoint(info.psiup, info.psidn)
    return user_stop_condition(info) # true stops after this sweep
end

psiup, psidn, energy = solve_hfdmrg(H, V, psiup0, psidn0;
    maxiter = 10, blocksize = 4, observer)
```

The qualified `HFDMRG.SweepInfo` passed to the observer contains `sweep`,
`energy`, full physical-basis `psiup` and `psidn`, and the built-in
`converged` flag. The orbitals are expanded after the local SCF update that
produces the reported end-of-sweep energy. Returning `true` requests a clean
stop; `false` or `nothing` continues unless the built-in convergence test has
succeeded. Observer exceptions propagate. The orbital arrays must be treated
as read-only and copied if mutable retained state is needed.

HFDMRG does not prescribe checkpoint formats or scientific measurements.
Observers own their I/O, diagnostics, and domain-specific stop rules, and their
work is included in total solver wall time. With `observer = nothing` (the
default), no observer is called.

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
psiup, psidn, energy = HFDMRG.solve_hfdmrg(H, backend, psiup0; maxiter = 5, block_partition = layout)
```

For ragged interactions, pass `Vblocks` instead of `V6`. With
`block_partition=layout`, whole-slice windows contract all nine ordered sectors
with zero steady allocation and fixed-rank work independent of total slices.
Any split window delegates its complete interaction to projection.
Each block folds absorbed slices into one exact pair-grouped block interaction
and keeps the two ordered channel orientations only for still-exterior slices.
Those channels use separate packed matrices with range-local ragged offsets;
no interaction symmetry is assumed.

## Choosing an interaction backend

Let `N` be the full basis size, `ns` the number of slices, `nj` a uniform
slice size, and `dims[s]` a ragged slice size. For a target residual, let `m`
be the target dimension, `P = m(m+1)/2`, and `w` the moving-window dimension.
The storage column describes the supplied global interaction representation;
all routes also build temporary and retained block/window caches.

| Route | Supplied-interaction storage | Window behavior | Recommended use |
|---|---:|---|---|
| Density-density `solve_hfdmrg(H, V, ...)` | `O(N^2)` | Projects the matrix interaction into dense retained-block tensors; never stores a global `N^4` tensor. | Default when the interaction is exactly density-density. This is the simplest and most mature route. |
| Density plus target residual `HFDMRG.DensityDensityTargetResidualBackend(V, Q, residual_pair)` | `O(N^2 + N*m + P^2)` | Adds target projection/lift work `O(w^2*m + w*m^2)` and dense-pair work `O(m^4)` to the base backend. | A fixed, small orthonormal target space needs a signed four-index correction. Do not use it as an arbitrary global interaction container. |
| Fixed sliced projection `HFDMRG.SlicedBasisBackend(layout, V6)` | `O(ns^2*nj^4)` | Reconstructs full-`N` density/Fock intermediates for each Fock build; correctness-first and relatively slow. | Uniform slice sizes, validation, small calculations, or comparison against the cached route. |
| Ragged sliced projection `HFDMRG.SlicedBasisBackend(layout, Vblocks)` | `O((sum_s dims[s]^2)^2)` | Uses the same full-basis projection strategy as fixed sliced data. | Slice sizes genuinely vary and a correctness/reference path is more important than speed. |
| Cached sliced `HFDMRG.SlicedBasisBackendCached(layout, V)` | Same sliced input storage, plus retained block contractions and local cross tensors. | Whole-slice windows use local ordered sectors with zero steady Fock allocation; split windows delegate entirely to projection. | Production sliced calculations with `block_partition=layout`. Pass `V6` when slices are uniform; use `Vblocks` only for ragged layouts. |

The convenience forms `solve_hfdmrg(H, layout, V, ...)` construct
the projection backend. Construct `HFDMRG.SlicedBasisBackendCached` explicitly
when the cached production path is intended. An exactly zero target residual
routes directly through the density-density backend.

## Advanced examples

The checked-in advanced examples use small deterministic models and validate
their own outputs:

- [`examples/observer_checkpoint.jl`](examples/observer_checkpoint.jl) records
  sweep history, writes each checkpoint through a temporary file, requests a
  clean stop, validates the saved orbitals and energy on readback, and uses the
  saved orbitals for a one-sweep restart. It uses Julia `Serialization` only as
  a local illustration; consumer projects should choose their own durable
  checkpoint format.
- [`examples/target_residual.jl`](examples/target_residual.jl) constructs a
  signed pair residual in a two-orbital target space and runs UHF through
  `DensityDensityTargetResidualBackend`.
- [`examples/cached_sliced.jl`](examples/cached_sliced.jl) builds a symmetric
  fixed-slice interaction and explicitly selects `SlicedBasisBackendCached`.

Run them from the repository root:

```sh
julia --project=. examples/observer_checkpoint.jl
julia --project=. examples/target_residual.jl
julia --project=. examples/cached_sliced.jl
```

## Repository layout
- `src/HFDMRG.jl`: module definition and include wiring.
- `src/core.jl`: generic HF-DMRG sweep engine.
- `src/backend_api.jl`: backend API contract for interactions.
- `src/backends/density_density.jl`: density-density backend for `V::Matrix`.
- `src/backends/density_density_target_residual.jl`: bounded target-residual
  correction on top of density-density interactions.
- `src/backends/sliced_basis.jl`: sliced-basis backend (projection-based).
- `src/backends/sliced_basis_cached.jl`: cached production sliced backend.
- `src/slice_layout.jl`: contiguous fixed or ragged slice layout helper.
- `docs/backend_architecture.md`: backend lifecycle, numerical contract,
  scaling policy, and extension checklist.
- `examples/`: executable introductory and advanced workflows.
- `test/`: compact numerical and end-to-end regression tests.
- `reference/HF_dmrg_legacy.jl`: historical implementation retained for
  reference; it is not an active regression oracle.

## How to run tests

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

For quick performance sanity, see `scripts/bench_sliced.jl`.

## Developer commands
- `make test`: run the package tests via `cjulia`.
- `make bench`: run `scripts/bench_sliced.jl`.
- `make verify`: run `scripts/verify_cached_vs_projection.jl`.

## Developing a backend

Read the [backend architecture and extension guide](docs/backend_architecture.md)
before changing an interaction route. It defines the block/window lifecycle,
RHF and UHF energy contract, state ownership, expected scaling analysis,
solver routing, and validation ladder.

New interactions should use that five-operation seam without copying the sweep
engine or constructing a global four-index tensor. General frameworks require
a demonstrated second consumer and a separate design review.
