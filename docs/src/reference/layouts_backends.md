# Layouts and backends

## Slice layouts

```@docs
HFDMRG.SliceLayout
HFDMRG.nslices
HFDMRG.slice_dims
HFDMRG.orb_range
HFDMRG.slice_local
```

`sum(layout.dims)` must match the physical one-body dimension. Passing a
`SliceLayout` as `block_partition` creates one physical chunk per slice.

## Sliced interaction storage

```@docs
HFDMRG.SlicedVeeRagged
```

For uniform slices, pass a six-dimensional `Array{Float64,6}` directly. For
ragged slices, either pass `Vblocks` to the constructor or create
`SlicedVeeRagged`. The ordered tensor convention is documented in
[Interaction backends](@ref).

## Projection sliced backend

Construct the reference/correctness route as:

```julia
backend = HFDMRG.SlicedBasisBackend(layout, V6)
# or
backend = HFDMRG.SlicedBasisBackend(layout, Vblocks)
```

The convenience `solve_hfdmrg(H, layout, V, ...)` forms this backend. It
supports fixed and ragged slices and windows that split a slice by expanding
the density and Fock contribution through the full physical basis.

## Cached sliced backend

```@docs
HFDMRG.SlicedBasisBackendCached
```

Construct it explicitly and pass `block_partition=layout` for aligned
production windows:

```julia
backend = HFDMRG.SlicedBasisBackendCached(layout, V6)
psiup, psidn, energy = solve_hfdmrg(H, backend, psiup0;
    block_partition = layout)
```

A split-slice window uses the projection route for the complete interaction.

## Density plus target residual

```@docs
HFDMRG.DensityDensityTargetResidualBackend
```

`Q` has shape `N × t` and orthonormal columns. `residual_pair` is symmetric
with shape `P × P`, `P=t(t+1)/2`, in canonical unnormalized pair order
`(p,q)`, `q ≤ p`. The residual is signed. `V`, `Q`, and `residual_pair` must
share one `AbstractFloat` element type. An exactly zero residual routes through
the ordinary density-density solver.
