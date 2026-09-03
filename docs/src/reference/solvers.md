# Solvers and controls

## Supported solve routes

| Problem | Supported form |
|---|---|
| Compact unit-cell RHF/UHF | `solve_hfdmrg(one_body::BandedOneBody, interaction::UnitCellInteraction; spin=:rhf or :uhf, kwargs...)` |
| Density-density RHF | `solve_hfdmrg(H, V, psiup0; kwargs...)` |
| Density-density common-H UHF | `solve_hfdmrg(H, V, psiup0, psidn0; kwargs...)` |
| Density-density split-H UHF | `solve_hfdmrg(Hup, Hdn, V, psiup0, psidn0; kwargs...)` |
| Explicit backend RHF/UHF | Replace `V` with a supported backend object |
| Sliced convenience RHF/UHF | Insert `layout, V6` or `layout, Vblocks` after the one-body argument(s) |

Passing the same matrix object for `Hup` and `Hdn` routes to the common-H
implementation. Distinct objects select split-H one-body caches and scale-aware
energy comparisons, even if their entries are numerically equal.

```@docs
HFDMRG.solve_hfdmrg
```

## Return values

Historical matrix and sliced routes return `(psiup, psidn, energy)`. Orbital
columns are in the full physical basis and correspond to the last completed
sweep. For RHF, `psiup == psidn`. The energy is electronic and uses the
conventions in [Numerical conventions](@ref).

The compact typed route returns `HFDMRGResult` for RHF or `UHFDMRGResult` for
UHF. These records retain compact block/root state and a `represented_energy`;
they do not materialize global occupied coefficients.

## Producer-energy audit

`producer_energy(result, producer)` accepts only a compact RHF/UHF result and
a `ProducerHamiltonian`. It streams the final determinant through the supplied
uncompressed one-body and density-density interaction, returning decomposed
one-body, Hartree, exchange, constant, total, error-envelope, and block-read
fields in `ProducerEnergy`.

The operation is observational and does not assert producer-Hamiltonian
stationarity. Historical tuple results do not dispatch to it, and neither
solver implementation falls back to the other.

## Keyword contract

The complete user-facing table and interpretation are in [Solver controls](@ref).
The keywords are:

```text
blocksize, block_partition, nblockcenter,
maxiter, cutoff, scf_cutoff, environment_cutoff,
frozen_occupied, observer, verbose
```

Not every specialized mode is supported by every backend. In particular,
`frozen_occupied=:roundoff_exact` is density-oriented and sliced routes reject
it. Private history policy is not a `solve_hfdmrg` keyword.

## Sweep observer information

```@docs
HFDMRG.SweepInfo
```

The observer contract is complete-sweep-only. The orbital arrays must be
treated as read-only. Returning `true` requests a clean stop after that sweep;
returning `false` or `nothing` continues. Any other return value is rejected.
