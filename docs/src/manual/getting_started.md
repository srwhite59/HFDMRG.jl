# Getting started

## Installation

HFDMRG is currently tested with Julia 1.12. From a source checkout:

```sh
cd /path/to/hfdmrg
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using HFDMRG; println("HFDMRG loaded")'
```

To use the checkout from another Julia environment:

```sh
julia --project=/path/to/my-project -e \
  'using Pkg; Pkg.develop(path="/path/to/hfdmrg")'
```

## Quick start

This deterministic RHF problem is also checked in as
`examples/quickstart.jl`; run the local file from the repository root.

```jldoctest quickstart
julia> using HFDMRG, LinearAlgebra, Random

julia> rng = MersenneTwister(7); N = 12; Nocc = 2;

julia> onsite = collect(range(-1.0, 1.0; length = N));

julia> hopping = diagm(1 => fill(-0.1, N - 1), -1 => fill(-0.1, N - 1));

julia> H = Matrix(Diagonal(onsite) + hopping);

julia> V = [0.4 / (1 + abs(i - j)) for i = 1:N, j = 1:N];

julia> psi0 = Matrix(qr(randn(rng, N, Nocc)).Q)[:, 1:Nocc];

julia> psiup, psidn, energy = solve_hfdmrg(H, V, psi0;
           maxiter = 4, blocksize = 2, cutoff = 1e-9,
           scf_cutoff = 1e-9, verbose = false);

julia> (round(energy; digits = 12), psiup == psidn,
        norm(psiup' * psiup - I) < 1e-10)
(-2.41361165233, true, true)
```

Run the file directly with:

```sh
julia --project=. examples/quickstart.jl
```

All public solve routes return `(psiup, psidn, energy)`. In RHF the two orbital
matrices are identical. In UHF their column counts independently fix
``N_\alpha`` and ``N_\beta``.

## Supported interaction choices

- `solve_hfdmrg(H, V, psi0; ...)` selects density-density RHF.
- Adding `psidn0` selects density-density UHF; distinct `Hup, Hdn` select the
  split-one-body UHF route.
- `SliceLayout` plus fixed or ragged sliced tensors selects the projection
  sliced route.
- Explicit `SlicedBasisBackendCached` construction selects the production
  cached sliced route.
- `DensityDensityTargetResidualBackend` adds a small signed target-space
  correction to density-density interactions.

See [Solvers and controls](@ref) for signatures and [Interaction backends](@ref)
for conventions and costs.
