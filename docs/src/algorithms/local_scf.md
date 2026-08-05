# Local RHF and UHF

## Status

The local SCF is the public numerical kernel used by all ordinary backends.
The same energy and density conventions apply to common- and split-one-body
UHF routes.

## Spaces and dimensions

At one moving window, ``H_B=B^THB`` and each density and Fock matrix is
``w\times w``. RHF has ``n`` occupied spatial orbitals. UHF has independent
``n_\alpha`` and ``n_\beta`` occupied subspaces, with those occupations fixed
for the solve.

## Inputs and outputs

The local solve receives projected one-body matrices, the current occupied
orbitals in window coordinates, one backend window object, and a persistent
mixing weight for that window position. It returns updated occupied orbitals,
the local energy, the number of microiterations, and the local convergence
result. The sweep core expands the orbitals before absorption and observation.

## Equations

For RHF, with one-spin density ``D``:

```math
F=H_B+G_{\mathrm{RHF}}(D),\qquad
E=\operatorname{Tr}[D(F+H_B)].
```

For UHF:

```math
F_\alpha=H_{\alpha,B}+G_\alpha(D_\alpha,D_\beta),\qquad
F_\beta=H_{\beta,B}+G_\beta(D_\alpha,D_\beta),
```

```math
E=\tfrac12\operatorname{Tr}[D_\alpha(F_\alpha+H_{\alpha,B})]
 +\tfrac12\operatorname{Tr}[D_\beta(F_\beta+H_{\beta,B})].
```

The backend adds only the interaction contribution. It must accept symmetric,
non-idempotent damped densities.

## Pseudocode and damping

```text
form determinant densities and build the initial Fock matrix/matrices
repeat at most four times
  diagonalize each Fock and occupy its lowest requested orbitals
  mix old and new determinant densities with the window's weight lambda
  build the post-density Fock matrix/matrices
  evaluate the trace-convention energy
  if energy rose, halve lambda for later updates at this window position
  stop locally when the energy change is below scf_cutoff
end
```

The post-density Fock is reused as the next diagonalization input. Therefore
``n`` microiterations require exactly ``n+1`` backend Fock builds, not ``2n``.
The mixing weight is associated with a window index and persists across
sweeps.

## Occupied eigensolves

The baseline operation diagonalizes a real symmetric Fock matrix. A private
guarded helper requests only the lowest occupied eigenpairs for a dense
`Matrix{Float64}` when

```math
n_{\mathrm{occ}}\le\left\lfloor w/8\right\rfloor.
```

Every other type or dimension uses the original full-spectrum operation. The
subset route changes neither the occupied rank nor the mathematical
eigenspace. At an exact degeneracy across the occupation boundary, individual
rank-``n`` projectors need not be unique; containment in the full degenerate
eigenspace is the meaningful invariant.

## Invariants

- Alpha and beta occupations remain fixed.
- RHF returns identical alpha and beta orbitals.
- Focks use the current complete damped density and the backend's documented
  direct/exchange convention.
- A post-update Fock is not discarded before the next microiteration.
- Returned full-basis orbitals are expanded after the final local SCF update.

## Convergence and stationarity

`scf_cutoff` tests the change in local damped-density energy, and each window
has a hard four-microiteration cap. The sweep can therefore move on without a
fully converged local SCF.

The outer `cutoff` compares complete-sweep energies. Neither test is the same
as the physical occupied--virtual residual

```math
R_\sigma=(I-C_\sigma C_\sigma^T)F_\sigma C_\sigma.
```

A soft response direction can have a small energy change and a material
projector error. High-consequence runs should audit full-space stationarity and
occupation ordering independently.

## Cost and scaling

For a full dense window, diagonalization is ``O(w^3)`` and stores ``O(w^2)``.
Forming densities costs ``O(w^2 n_\sigma)``. The guarded subset operation
reduces eigensolver work when the occupied fraction is small, but backend Fock
costs can still dominate. At most five Fock builds occur per window.

## Code map

- `src/core.jl`: ordinary RHF/UHF loops, damping, subset guard, energy, and
  convergence.
- `src/backend_api.jl`: in-place Fock addition contract.
- `src/backends/*.jl`: representation-specific direct and exchange terms.

## Current deviations and non-goals

- Local iterations use energy stopping rather than a direct residual target.
- HFDMRG does not perform a stability analysis or prove the global HF minimum.
- The correctness-oriented frozen-occupied route is not presented as a general
  local-solve optimization.
