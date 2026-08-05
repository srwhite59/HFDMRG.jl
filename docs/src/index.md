# HFDMRG.jl

HFDMRG solves restricted and unrestricted Hartree--Fock problems by sweeping a
moving variational window through an ordered one-particle basis. It returns one
Slater determinant. The DMRG analogy concerns the moving window and compressed
environment bases; HFDMRG does **not** form a correlated matrix-product state.

The package separates two responsibilities:

- the solver core owns block geometry, the two sweep directions, local SCF,
  occupied-orbital compression, convergence, and observation; and
- an interaction backend owns its compact interaction representation and the
  block/window contractions needed to add Fock terms.

## Where to begin

- [Getting started](@ref) gives installation and a checked RHF calculation.
- [Solver controls](@ref) explains the public keywords and how to assess a run.
- [Moving-window sweep](@ref) describes the algorithm and its direct-sum
  geometry.
- [Interaction backends](@ref) compares the density, target-residual, sliced
  projection, and cached sliced representations.
- [Examples](@ref) orders the executable examples from basic use through
  checkpointing and specialized interactions.
- [Solvers and controls](@ref) and [Layouts and backends](@ref) are the curated
  API reference.

## Current status

The public surface is intentionally narrow: `solve_hfdmrg` is the only exported
name, while supported layout and backend constructors are accessed with the
`HFDMRG.` qualifier. History acceleration is implemented as a private,
provisional research path and is documented only to explain its numerical
contract. It has no public calling syntax.

HFDMRG presently has no hosted documentation deployment. Build this manual
locally with:

```sh
julia --project=docs docs/make.jl
```

The generated site is written to `docs/build/`.
