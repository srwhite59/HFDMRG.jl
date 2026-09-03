# Interaction backends

## Status

Four interaction families share the public sweep: density-density, density plus
a target-space residual, sliced projection, and cached sliced. Backend types
are qualified with `HFDMRG.`. The compact input/result records and the
`solve_hfdmrg` and `producer_energy` operations form the exported surface.

## Common inputs and outputs

A solve-wide backend holds a compact global representation. A block state holds
the representation projected into one retained environment. A window state
combines left and right blocks with the explicit center and may own Fock
scratch. The core passes a symmetric local density and the backend adds the
interaction term in place to one or two Fock matrices.

The output must obey the trace convention in [Numerical conventions](@ref).
Backends do not compute an independent total energy, choose occupations, stop
the sweep, or write checkpoints.

## Density-density

### Representation and algebra

An ``N\times N`` matrix ``V`` represents

```math
(ij\mid kl)=\delta_{ij}\delta_{kl}V_{ik}.
```

The RHF and UHF formulas are given in [Numerical conventions](@ref). Blocks
retain physical-to-pair projections, outside couplings, and the complete
retained ``k^4`` interaction; no global ``N^4`` tensor is constructed.

### Cost

- global input: ``O(N^2)``;
- typical block state: ``O(|r|k^2+|r_V|k^2+k^4)``;
- local dense interaction work depends on ``k_L,c,k_R`` rather than a global
  four-index tensor.

This is the simplest and most mature route when the interaction is exactly
density-density.

## Density plus target residual

### Spaces and representation

Let ``Q\in\mathbb{R}^{N\times t}`` have orthonormal columns. Canonical
unnormalized symmetric pair coordinates ``(p,q)``, ``q\le p``, number
``P=t(t+1)/2``. A symmetric signed ``P\times P`` matrix stores the residual
``(pq\mid rs)``. `V`, `Q`, and the residual must share one floating element
type; the large `V` is referenced while the compact target data are owned.

At a window, the target basis is projected recursively. The local density is
projected into target coordinates, contracted in pair space, and lifted back
to the moving window. Direct and exchange signs and factors follow the same
RHF/UHF trace convention as the base backend.

### Cost and non-goals

- global input: ``O(N^2+Nt+P^2)``;
- target block state: ``O(kt)`` beyond the density backend;
- window projection/lift: ``O(w^2t+wt^2)``;
- dense target-pair work: ``O(t^4)``.

This route is for a fixed small target correction. It is not an arbitrary
additive backend framework and never expands to a global four-index tensor. An
exactly zero residual routes through the original density backend.

## Sliced projection

### Spaces and tensor convention

For fixed slice size ``d``, `V6[a,b,c,d,n,m]` stores an unsymmetrized ordered
sector with

```math
p=(n,a),\quad q=(m,c),\quad r=(n,b),\quad s=(m,d),
```

so that the element is ``(pq\mid rs)``. Ragged
`Vblocks[n][m]` uses dimensions
``d_n\times d_n\times d_m\times d_m`` with the same orientation. Distinct
ordered orientations must not be identified without a validated physical
symmetry.

The projection backend stores the retained physical basis in each block. At
every Fock call it expands the window density to the full physical basis,
contracts the sliced interaction, and projects the Fock contribution back.
This is a direct correctness route and supports windows that cut through a
slice.

### Cost

Let ``S_2=\sum_s d_s^2``. The supplied interaction stores ``O(S_2^2)`` values
(``O(n_s^2d^4)`` for uniform slices). A Fock call uses full-physical
``O(N^2+Nw)`` temporaries and loops over all ordered slice pairs.

## Cached sliced

The cached route shares the identical unsymmetrized input convention. For
slice-aligned blocks it folds absorbed slices into one exact pair-grouped block
interaction and stores separate packed `W` and `Wswap` channel orientations
only for still-exterior slices. Absorption updates the old/new sectors before
discarding newly internal channels.

At a whole-slice window, all nine ordered left/center/right sectors are
assembled from the compact block states. Local Fock work is independent of the
total slice count at fixed retained ranks and center. A window that splits a
slice delegates its **entire** interaction to projection; cached and projection
sectors are never mixed in one window.

With ``S_{2,\mathrm{ext}}`` restricted to exterior slices, representative block
storage is

```math
O(|r|k+k^4+k^2S_{2,\mathrm{ext}}).
```

Whole-slice hot Fock kernels have zero steady allocation in the accepted route,
but the backend spends more persistent memory than projection. Use
`block_partition=layout` for the aligned production path. Complete edge ranks
can make their ``k^4=d^4`` state important.

## Lifecycle pseudocode

```text
construct backend from compact global interaction
initialize complete fixed edge block states
initialize compressed interior block states
for each moving window
  combine left block, center, and right block into one window state
  reuse that window's scratch for all local Fock builds
  after local SCF, absorb one center chunk with final recurrence maps
  return one replacement block state
end
```

## Invariants

- Caller arrays are read-only.
- Block caches describe the final stored reorthogonalized environment basis.
- RHF uses ``2J-K``; UHF uses total-density direct and same-spin exchange.
- Damped non-idempotent densities are valid inputs.
- Sliced order and orientation are preserved explicitly.
- No supported backend creates a global dense ``N^4`` interaction.

## Code map

- `src/backends/density_density.jl`
- `src/backends/density_density_target_residual.jl`
- `src/backends/sliced_basis.jl`
- `src/backends/sliced_basis_cached.jl`
- `src/slice_layout.jl`
- `src/backend_api.jl`

## Current deviations and non-goals

- Projection sliced is a correctness/reference path, not the preferred large
  production route.
- Cached split-slice windows intentionally fall back as a whole.
- Target residual is bounded to one small orthonormal target space.
- There is no generic additive-backend, factorization, MPO, or global
  four-index interface.
