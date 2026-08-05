# Moving-window sweep

## Status

This is the public solver algorithm shared by every supported interaction
backend. Backends do not copy or own the sweep engine.

## Spaces and dimensions

HFDMRG partitions an ordered physical one-particle basis of dimension ``N``
into contiguous left, center, and right regions. If the retained environment
bases have ranks ``k_L`` and ``k_R`` and the explicit center has ``c`` physical
functions, the local one-particle space has dimension

```math
w=k_L+c+k_R.
```

Its embedding in the physical basis is

```math
B=\begin{bmatrix}\Phi_L&0&0\\0&I_c&0\\0&0&\Phi_R\end{bmatrix},
```

with ordered local coordinates `[left; center; right]`.

This construction is a **direct sum of one-particle spaces**. It is not the
tensor product of many-body block Hilbert spaces used by conventional DMRG.
The state remains one Slater determinant and there is no user-selected
many-body bond dimension.

## Inputs and outputs

The sweep receives one or two real symmetric one-body matrices, an interaction
backend, orthonormal occupied orbitals, a physical partition, and solver
controls. It returns full-physical-basis alpha and beta occupied orbitals and
the end-of-sweep electronic energy.

The first and last physical chunks are permanent boundary blocks. Their bases
are complete coordinate identities, so unoccupied edge directions cannot be
discarded at initialization. Grown interior environments are compressed from
the occupied-orbital coefficients with the requested `environment_cutoff`.

## Sweep pseudocode

```text
partition the ordered physical basis into chunks
build left and right initial block families
  fixed edge blocks keep complete identity spans
  grown interior blocks use occupied-orbital SVD compression
project the initial determinant into the first window

for each complete sweep
  for each left-to-right window, then each right-to-left window
    assemble the projected one-body operator
    construct one backend window state
    run at most four local RHF or UHF microiterations
    expand the post-SCF determinant into the physical basis
    compute an occupied-orbital SVD for the departing center chunk
    reorthogonalize the stored new environment basis
    absorb the chunk and construct the one successor block/cache
    transform the determinant into the next window coordinates
  end
  test end-of-sweep energy change
  notify the observer with the matching full-basis determinant and energy
end
```

The environments are replaced on the fly as the window advances. Ordinary
operation does **not** perform a separate global block/cache rebuilding sweep.
Initialization constructs block families at solver entry; each later move
constructs only the successor block required by that move.

## Linear algebra

For a growing side, the restricted occupied coefficient matrix is factorized
by an SVD. Every left singular vector with
``\sigma>\mathtt{environment\_cutoff}`` is retained, with at least one vector
when none pass. A small reorthogonalization produces the final stored
``\Phi_{\mathrm{new}}``. The recurrence maps supplied to the one-body and
interaction caches describe this final basis, not the preliminary SVD basis.

The fixed edges deliberately do not use this occupied-only truncation.

## Invariants

- ``\Phi_L^T\Phi_L=I`` and ``\Phi_R^T\Phi_R=I`` to numerical precision.
- The center basis is the physical coordinate identity.
- The expanded occupied projector is preserved by a coordinate move, up to
  the declared environment truncation and roundoff.
- One-body and interaction cache recurrences use the same final stored basis.
- Observer orbitals correspond to the local SCF update that produced the
  reported end-of-sweep energy.
- A backend owns interaction state only; it does not choose occupations,
  sweep direction, damping, or stopping.

## Convergence semantics

One left-to-right pass plus one right-to-left pass is a sweep. The public
convergence flag tests the change between end-of-sweep energies. It does not
prove full-space Fock stationarity, Aufbau ordering, uniqueness, or global
minimality. See [Local RHF and UHF](@ref).

`nblockcenter=1` is the default. Increasing it retains consecutive chunks in
the explicit center, so adjacent windows overlap. This can improve propagation
per nominal sweep, but larger local systems and additional local updates can
make it slower overall. It is a scientific geometry control, not a universal
performance knob.

## Cost and scaling

At one window, the solver-independent dense eigensolve cost is at most
``O(w^3)``; occupied-subset diagonalization can reduce that constant when its
guard applies. Density construction and projected one-body work are at least
``O(w^2 n_{\mathrm{occ}})``. Backend-specific Fock and cache costs are listed
in [Interaction backends](@ref).

The retained ranks are data dependent. For a fixed edge of physical dimension
``d``, completeness forces rank ``d``; backends with a retained ``k^4`` block
term can therefore require ``O(d^4)`` edge storage. Interior ranks follow the
occupied-orbital singular spectrum rather than an independent bond dimension.

## Code map

- `src/core.jl`: partitioning, initialization, sweep schedule, local solve,
  occupied SVD, one-body recurrence, convergence, and observer calls.
- `src/backend_api.jl`: the five-operation interaction seam.
- `src/backends/`: backend-specific block and window state.

## Current deviations and non-goals

- The sweep advances by one physical chunk; staggered or stride-two schedules
  are not implemented.
- HFDMRG does not search all Hartree--Fock basins.
- It does not construct an MPS or a many-body reduced density matrix.
- Private research features may restart the solver between sweep segments, but
  that is not the ordinary in-sweep environment lifecycle.
