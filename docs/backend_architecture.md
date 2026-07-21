# HFDMRG Backend Architecture

HFDMRG has one sweep engine and several interaction backends. The sweep engine
decides where the moving window is, transforms occupied orbitals, runs the
local SCF, computes energies, and decides when to stop. A backend translates
one particular interaction representation into the block and window objects
needed by that engine.

This is the in-repository extension contract. `solve_hfdmrg` is the only
exported name; the functions in [`src/backend_api.jl`](../src/backend_api.jl)
are an architectural seam for HFDMRG development, not a promise that private
state types or helper names are stable public API.

## The three backend lifetimes

| Object | Lifetime | Responsibility |
|---|---|---|
| Backend | One or more complete solves | Validate and hold the compact global interaction representation. Treat supplied arrays as read-only and make ownership of any copies explicit. |
| Block state | From block construction until that block is replaced during a sweep | Hold the interaction projected into one retained left or right environment and any data needed to absorb the next center chunk. |
| Window state | One left-center-right position | Combine the two block states with the explicit center and own scratch reused by the local SCF Fock builds. |

The core does not require an abstract backend supertype. Multiple dispatch on
the concrete backend, block-state, and window types selects the implementation.
The fallback methods in `backend_api.jl` fail clearly when a required operation
is missing.

## Core/backend boundary

The sweep core owns:

- physical partitioning and left-to-right/right-to-left scheduling;
- projected common or spin-dependent one-body Hamiltonians;
- occupied orbitals, densities, diagonalization, occupation, and damping;
- the SVD that chooses each new environment basis;
- RHF/UHF energy evaluation, convergence, and observer notification; and
- expansion of the final orbitals into the physical basis.

An interaction backend owns:

- validation and interpretation of its solve-wide interaction data;
- interaction projections and caches for retained blocks;
- construction and scratch for one moving interaction window; and
- in-place addition of RHF or UHF interaction terms to the Fock matrices.

A backend does not choose occupations, compute a separate total energy, stop a
sweep, write checkpoints, or transform the one-body Hamiltonian. If a feature
needs only interaction data at the current window, it belongs behind this seam.

## What happens during a sweep

At one position, the physical basis has the ordered form

```text
left physical range | center physical range | right physical range
        phi_L       |          I            |        phi_R
```

`phi_L` and `phi_R` have orthonormal columns and map physical rows into the
retained environment bases. If their column counts are `kL` and `kR` and the
center contains `c` sites, the local window dimension is
`w = kL + c + kR`.

The core then follows this lifecycle:

1. Build initial left and right block states with `vee_init_block`.
2. At each position, assemble the projected one-body matrix and call
   `vee_window` once.
3. Build one Fock matrix for the initial density, then run at most four local
   SCF micro-iterations. Each iteration builds one post-density Fock matrix for
   both its energy and the next diagonalization, so `n` iterations make `n+1`
   backend calls and one window receives at most five.
4. Expand the updated orbitals to the physical basis, retaining the latest
   expansion for end-of-sweep return and observer/checkpoint use.
5. Obtain a new environment basis from an occupied-orbital SVD, absorb the
   departing center chunk with `vee_absorb_block`, and move the window.
6. Complete both directions, test the end-of-sweep energy change, and notify
   the observer.

Backends do not own this schedule. They must not copy the sweep engine to gain
an interaction feature.

## Required backend operations

Use these names and argument positions exactly:

```julia
vee_init_block(side, ra, raV, phi, backend)
vee_absorb_block(side, old, cra, Phi_old, Phi_C, phi_new, raV_new, backend)
vee_window(left, right, Cra, backend)
vee_add_fock_r!(F, rho, window)
vee_add_fock!(Fup, Fdn, rhoup, rhodn, window)
```

The state returned by one operation is opaque to the core; only later methods
for the same backend need to understand it.

### Initialize a block

`vee_init_block` receives:

- `side`, either `:left` or `:right`;
- `ra`, the contiguous physical range represented by the block;
- `raV`, the physical range outside the block that can still couple to it;
- `phi`, an `length(ra) x k` orthonormal block basis; and
- the solve-wide backend object.

A left block grows from low physical indices and has
`raV = (ra[end]+1):N`. A right block grows from high indices and has
`raV = 1:(ra[1]-1)`. The returned state should contain only
representation-specific information needed by later absorption and window
construction.

### Absorb a center chunk

`vee_absorb_block` constructs a replacement state after the SVD moves a center
chunk into an environment. `Phi_old` has shape `k_old x k_new`; `Phi_C` has
shape `length(cra) x k_new`. Before the core's small reorthonormalization, the
new physical basis has the block form

```text
left growth:  [phi_old * Phi_old; Phi_C]
right growth: [Phi_C; phi_old * Phi_old]
```

The core also supplies the final `phi_new` and new outside range `raV_new`.
A backend may transform an old cache recursively with `Phi_old` and `Phi_C`,
or rebuild from `phi_new` when that is the clearer and measured choice. Tests
for a recurrence should compare both left and right growth with a direct
physical-basis projection.

Do not assume the retained dimension stays constant. It is set by the
occupied-orbital SVD and can change from one block to the next.

### Construct a window

`vee_window(left, right, Cra, backend)` returns the interaction object for the
ordered local basis

```text
[left retained basis; center physical sites; right retained basis].
```

This is the right place to combine persistent block caches and allocate scratch
that will be reused by the local SCF. A sliced backend must handle center ranges
that cut through a slice unless its public constructor explicitly rules that
case out.

### Add Fock contributions

The core initializes each Fock matrix from its projected one-body Hamiltonian.
The backend methods add interaction terms in place; they must not overwrite the
one-body contribution.

- `vee_add_fock_r!` adds the RHF `2J - K` contribution for the one-spin density
  `rho`.
- `vee_add_fock!` adds total-density direct and like-spin exchange contributions
  to `Fup` and `Fdn`.

All Fock and density matrices have shape `w x w`. Densities can be damped
mixtures and therefore need not be idempotent. A backend must support that
general symmetric-density input rather than assuming `rho^2 == rho`.

These are hot methods. Reuse window scratch, use `mul!` or explicit auditable
loops where appropriate, and measure steady allocations. Do not retain a Fock
or density argument after the call.

## Numerical contract

For a window density `D`, the core uses

```text
F_RHF = H_B + G_RHF(D)
E_RHF = Tr[D (F_RHF + H_B)]

F_alpha = H_alpha,B + G_alpha(D_alpha, D_beta)
F_beta  = H_beta,B  + G_beta(D_alpha, D_beta)
E_UHF = 1/2 Tr[D_alpha (F_alpha + H_alpha,B)]
      + 1/2 Tr[D_beta  (F_beta  + H_beta,B)].
```

Consequently:

- direct terms use `D_alpha + D_beta` in UHF;
- exchange terms use only the density of the same spin;
- RHF uses a one-spin density and `2J - K`;
- backend output must be consistent with real symmetric Fock matrices; and
- the backend performs no unit conversion and adds no nuclear or
  consumer-defined scalar energy.

The input matrices and orbitals are real. A constructor should reject invalid
dimensions, scalar types, non-finite data, or violated structural symmetry as
early as practical. Floating structural checks should use fixed scale-aware
internal tolerances; do not expose validation-detail tolerance keywords without
a demonstrated user requirement.

Every new tensor convention needs an explicit four-index RHF and UHF oracle.
Index names that look familiar are not proof that direct, exchange, pair
weights, and factors of two agree with the core convention. The public
[numerical conventions](../README.md#numerical-conventions) and the
[target-residual formulas](../notes/target_space_interaction_residual_design_2026-07-17.md#fock-and-energy-formulas)
show the existing contracts.

## Solver overloads and routing

An in-tree backend normally supplies restricted, common-one-body UHF, and
split-one-body UHF overloads:

```julia
function solve_hfdmrg(H, backend::MyBackend, psiup0; kwargs...)
    solve_hfdmrg_core(H, backend, psiup0, psiup0;
        restricted = true, kwargs...)
end

function solve_hfdmrg(H, backend::MyBackend, psiup0, psidn0; kwargs...)
    solve_hfdmrg_core(H, backend, psiup0, psidn0;
        restricted = false, kwargs...)
end

function solve_hfdmrg(Hup, Hdn, backend::MyBackend, psiup0, psidn0; kwargs...)
    Hup === Hdn && return solve_hfdmrg(Hup, backend, psiup0, psidn0; kwargs...)
    solve_hfdmrg_core_split(Hup, Hdn, backend, psiup0, psidn0; kwargs...)
end
```

Add the backend source to `src/HFDMRG.jl`. Backend constructors currently
remain qualified as `HFDMRG.MyBackend`; do not export a new name merely for
convenience. Add a representation-specific convenience overload only when it
is clearer than explicit backend construction and has a demonstrated user.

Preserve exact compatibility routes. For example, an exactly zero target
residual delegates to the original density-density overload before block or
window state is built. Do not silently drop a merely small correction.

## Ownership and mutation rules

- Large supplied interactions should normally be referenced, not converted or
  copied. If a backend needs owned normalized data, make that constructor cost
  explicit. The target-residual backend references `V` but owns copies of `Q`
  and `residual_pair`.
- Treat all arrays supplied by the caller as read-only.
- Block states own caches that survive across window positions. Replacing a
  block returns a new valid state; it must not leave aliases to temporary SVD
  workspaces.
- Window states may own mutable scratch because the sweep is serial. Scratch
  is private to that window and must be reinitialized correctly on each Fock
  call.
- Do not add mutable global numerical state. Instrumentation must be optional
  and must not change results.
- Consumer checkpoint formats, diagnostics, and stop policy belong in an
  observer or wrapper, not in an interaction backend.

## Storage and performance model

Let `N` be the full basis size, `k` a retained block dimension, `c` a center
size, and `w` the window dimension. For sliced data, let `ns` be the number of
slices, `nj` a uniform slice size, `d_s` each ragged slice size, and
`S2 = sum_s d_s^2`. For a target space, let `t` be its dimension and
`P = t(t+1)/2`.

The leading storage terms are:

| Backend | Global interaction | Persistent block state | Window/Fock behavior |
|---|---:|---:|---|
| Density-density | `O(N^2)` | `O(length(ra)*k^2 + length(raV)*k^2 + k^4)` | Dense retained-block and cross-block tensors; no global `N^4` object. |
| Density plus target residual | Base plus `O(N*t + P^2)` | Base plus `O(k*t)` | Projection/lift `O(w^2*t + w*t^2)` and exact pair contraction `O(t^4)`. |
| Sliced projection | Fixed `O(ns^2*nj^4)` or ragged `O(S2^2)` | Primarily the block basis | Builds full-`N` density and Fock intermediates, requiring `O(N^2 + N*w)` temporary storage per Fock route. |
| Cached sliced | Same sliced input | About `O(k^4 + k^2*S2 + N*k)` | Whole-slice windows retain `O(kL^2*kR^2)` cross tensors and reference block caches. Their local Fock work is independent of `ns` at fixed ranks and center; split windows delegate entirely to projection. |

These terms describe interaction-specific storage, not the complete solver.
Actual time depends on block ranks, center placement, slice sizes, BLAS shape,
and the number of local SCF builds. The cached sliced route deliberately spends
more persistent memory than the projection route to avoid repeated full-basis
work. Neither route should create a global four-index interaction.

The accepted cached-sliced design and frozen-Be evidence are preserved in the
[final performance account](../notes/cached_sliced_performance_final_2026-07-20.md).

A design proposal for a new backend must state:

- global, per-block, and per-window storage;
- initialization, absorption, window, and per-Fock work;
- expected steady allocation in the hot Fock method;
- a representative problem size and predicted memory; and
- why the existing backends cannot express the requirement.

Do not introduce a generic additive backend, factorization interface, or
global four-index container for a single consumer. A second concrete use case
and separate design review are required before generalizing the seam.

## Validation ladder

Keep committed tests compact and numerical. A new or materially changed
backend should normally pass this ladder:

1. Constructor rejection tests for important shape, type, finiteness, and
   structural-invariant failures.
2. A tiny signed explicit four-index oracle for RHF and UHF Fock matrices and
   energies.
3. Direct physical-basis comparison of both left and right block absorption.
4. Window-projection parity for random orthonormal block bases and symmetric,
   possibly non-idempotent densities.
5. Exact routing parity for any claimed zero or compatibility case.
6. A small end-to-end sweep against an independent established route, including
   common- and split-one-body cases when both are supported.
7. The complete package test suite.
8. A same-process, warmed representative benchmark reporting wall time,
   allocations, and retained memory. Use a real consumer fixture when the
   feature exists for one.
9. Consumer workflow acceptance before retiring an external oracle or copied
   implementation.

Heavy fixtures, restart files, PID files, and logs belong in machine-local
scratch, not the package. Tests should protect formulas and public behavior;
do not test private helper names, scratch-buffer shapes, or internal
vocabulary.

## Implementation checklist

Before editing source:

1. Write the interaction convention and RHF/UHF formulas.
2. Choose backend, block-state, and window-state representations.
3. State scaling, performance gates, validation fixtures, and a line budget.
4. Confirm that the five backend operations can express the feature. If they
   cannot, stop for architecture review rather than copying or widening the
   sweep engine.

During implementation:

1. Add one bounded file under `src/backends/` and include it from
   `src/HFDMRG.jl`.
2. Validate the solve-wide representation in its constructor.
3. Implement the five backend operations and the necessary solver overloads.
4. Keep hot scratch in the window and persistent projections in block state.
5. Add the smallest numerical tests and one user-facing example that protect
   the live contract.

At review, report the formulas proved, scaling, representative timing and
memory, tests actually run, exact changed files and line counts, any code made
obsolete, and the final git status. Do not modify `src/core.jl` when the
backend seam already expresses the operation.
