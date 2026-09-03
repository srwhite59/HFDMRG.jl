# HFDMRG

HFDMRG solves restricted and unrestricted Hartree--Fock problems with
DMRG-style sweeps over an ordered one-particle basis. It returns one Slater
determinant, not a correlated matrix-product state.

At each position, the solver forms a direct-sum working space

```text
retained left environment | explicit physical center | retained right environment
```

runs a short local SCF calculation, and absorbs one center chunk into the next
environment. Interaction backends supply compact density-density,
target-residual, or sliced contractions without requiring a global four-index
tensor.

## Installation

HFDMRG is currently tested with Julia 1.12. From a source checkout:

```sh
cd /path/to/hfdmrg
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using HFDMRG; println("HFDMRG loaded")'
```

To use the checkout from another Julia project:

```sh
julia --project=/path/to/my-project -e \
  'using Pkg; Pkg.develop(path="/path/to/hfdmrg")'
```

## Quick start

```julia
using HFDMRG
using LinearAlgebra
using Random

rng = MersenneTwister(7)
N, Nocc = 12, 2
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

The checked-in `examples/quickstart.jl` prints:

```text
energy = -2.41361165233
```

Historical dense and sliced routes return `(psiup, psidn, energy)`, with
full-physical-basis occupied orbitals and the electronic energy of the supplied
Hamiltonian. Historical RHF returns identical alpha and beta orbitals.

The compact unit-cell route is selected only by
`solve_hfdmrg(one_body::BandedOneBody, operator::UnitCellInteraction; ...)`.
It returns an `HFDMRGResult` containing represented working energy,
convergence diagnostics, and a compact block/root state handle. It does not
materialize global occupied coefficients. The compact represented energy and
the historical dense/sliced electronic energy are distinct result contracts.

## Supported routes

- `solve_hfdmrg(one_body::BandedOneBody,
  operator::UnitCellInteraction; ...)`: compact block-native unit-cell RHF.
- `solve_hfdmrg(H, V, psiup0; ...)`: density-density RHF.
- `solve_hfdmrg(H, V, psiup0, psidn0; ...)`: common-one-body UHF.
- `solve_hfdmrg(Hup, Hdn, V, psiup0, psidn0; ...)`: split-one-body UHF.
- `solve_hfdmrg(H, layout, V6_or_blocks, ...; ...)`: fixed or ragged sliced
  projection.
- `HFDMRG.SlicedBasisBackendCached(layout, V)`: cached production sliced
  contractions when used with `block_partition=layout`.
- `HFDMRG.DensityDensityTargetResidualBackend(V, Q, residual_pair)`: a small
  signed four-index correction in a fixed orthonormal target space.

`solve_hfdmrg` is the only exported name. Supported layout and backend types
are accessed with the `HFDMRG.` qualifier.

Argument types form a strict dispatch boundary: typed compact unit-cell calls
never enter the dense/sliced solver, and matrix/sliced calls never enter the
compact solver.

## Algorithm in one page

1. Partition the ordered physical basis into contiguous chunks. The two fixed
   edge chunks retain complete coordinate-identity spans.
2. Build compressed interior environments from the restricted occupied-orbital
   SVD. `environment_cutoff` controls their retained singular directions.
3. Transform the one-body operator and construct the backend state for one
   left--center--right window.
4. Perform at most four local RHF or UHF microiterations, with density damping
   after energy rises.
5. Expand the updated determinant, absorb the departing center chunk, and
   construct the one successor environment as the sweep advances.
6. Complete both directions, test the end-of-sweep energy change, and call the
   optional complete-sweep observer.

The one-particle spaces are combined by direct sum, not by a many-body tensor
product. Ordinary operation updates environments on the fly and does not run a
separate global block/cache rebuilding sweep.

`nblockcenter=1` is the default. Overlapping multi-block centers can improve
propagation per nominal sweep, but larger local work means they are not
universally faster.

## Convergence

The built-in `converged` flag tests end-of-sweep energy change. Reaching
`maxiter` without satisfying that test returns the last completed-sweep
determinant rather than throwing an exception.

Energy convergence is not, by itself, proof of full-space Fock stationarity,
Aufbau ordering, uniqueness, or the global HF minimum. Consequential runs
should also inspect occupied--virtual Fock residuals, Gram errors, gaps or
stability information, and sensitivity to the initial determinant and window
geometry.

## Manual

HFDMRG now has a local Documenter manual with six primary sections:

- [Home](docs/src/index.md)
- [Manual](docs/src/manual/index.md)
- [Algorithms](docs/src/algorithms/index.md)
- [Examples](docs/src/examples/index.md)
- [Reference](docs/src/reference/index.md)
- [Developer Notes](docs/src/developer/index.md)

Especially useful pages are:

- [Getting started](docs/src/manual/getting_started.md)
- [Solver controls](docs/src/manual/solver_controls.md)
- [Numerical conventions](docs/src/manual/numerical_conventions.md)
- [Moving-window sweep](docs/src/algorithms/sweep.md)
- [Local RHF and UHF](docs/src/algorithms/local_scf.md)
- [Interaction backends](docs/src/algorithms/interaction_backends.md)
- [Private/provisional history acceleration](docs/src/algorithms/history_acceleration.md)
- [Backend architecture](docs/src/developer/backend_architecture.md)

Build the site locally; no deployment is configured:

```sh
julia --project=docs docs/make.jl
```

## Examples

The recommended sequence is:

```sh
julia --project=. examples/quickstart.jl
julia --project=. examples/basic_uhf.jl
julia --project=. examples/observer_checkpoint.jl
julia --project=. examples/target_residual.jl
julia --project=. examples/cached_sliced.jl
```

The [Examples](docs/src/examples/index.md) page explains what each workflow
demonstrates. The basic UHF example independently reconstructs the physical
Focks, energy, residuals, gaps, and Aufbau projectors.

## Tests and developer checks

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=. scripts/verify_cached_vs_projection.jl
```

The Makefile also provides `make test`, `make verify`, and `make bench`.
Backend development must follow the five-operation seam and numerical contract
in the [backend architecture guide](docs/src/developer/backend_architecture.md),
rather than copying the sweep engine or adding a global four-index interaction.

Private history acceleration is documented for maintainers but has no public
calling syntax. Frozen-occupied operation is not presented as a recommended
general optimization. License, citation metadata, CI, hosted documentation,
and a public acceleration interface remain separate release decisions.
