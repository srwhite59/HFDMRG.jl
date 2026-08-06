# Target-Space Interaction Residual: Final Design And Evidence

Date: 2026-08-05
Implementation commit: `64def84`
Status: production-supported specialized backend

## Outcome And Boundary

`HFDMRG.DensityDensityTargetResidualBackend(V, Q, residual_pair)` adds one
small signed four-index correction to the density-density interaction. The
correction lives in a fixed orthonormal target space and is propagated through
the ordinary backend seam. It does not copy the sweep engine, form a physical
`N^4` tensor, or introduce a generic additive-backend or factorization API.

The backend is qualified through `HFDMRG`; `solve_hfdmrg` remains the only
exported name. RHF, common-one-body UHF, and split-one-body UHF are supported.
An exactly zero residual delegates to the original density route, preserving
its methods, operation order, orbitals, and energy exactly.

Static screened-reference subtraction is outside this backend. A caller may
put a reference field in `H` and account for a matching constant separately;
HFDMRG owns only the state-dependent residual and returns its normal electronic
energy.

## Representation And Validation

Let `N` be the physical dimension, `m` the target dimension, and
`P = m(m+1)/2`. `Q` is `N x m` and maps target coefficients into the physical
basis. Its columns must be orthonormal. The residual uses chemists' notation

```text
R[p,q,r,s] = residual_pair[pair(p,q), pair(r,s)]
pair(p,q) = max(p,q)(max(p,q)-1)/2 + min(p,q).
```

Pairs are unnormalized: diagonal pair-density entries have weight one and
off-diagonal entries have weight two. `residual_pair` must be symmetric but
may be indefinite because it is a signed difference, not a Coulomb metric.

`V`, `Q`, and `residual_pair` must share one element type
`T <: AbstractFloat`. `V` is referenced without conversion or copying. The
backend owns dense copies of `Q` and `residual_pair`; after validating the
latter, it removes roundoff antisymmetry by averaging it with its transpose.
Constructor checks are fixed internal rules, not public controls:

```text
norm(Q'Q-I, Inf) <= 128 eps(T) max(N,m)
maxabs(R-R')     <= 128 eps(T) max(1,maxabs(R)).
```

Nonfinite values, incorrect shapes, integer or mixed element types, material
asymmetry, and nonorthonormal `Q` are rejected.

## Block And Window Recurrence

For a retained block on physical range `ra` with basis `phi`, store only

```text
Q_B = phi' Q[ra,:].
```

When a center chunk `cra` is absorbed, the backend API supplies maps into the
final stored, post-reorthogonalization basis. Both directions therefore use

```text
Q_B,new = Phi_old' Q_B,old + Phi_C' Q[cra,:].
```

The window map in `[left; center; right]` order is

```text
Q_w = vcat(Q_L, Q[Cra,:], Q_R).
```

The recurrence was checked against direct `phi_new'Q` after repeated left and
right absorptions. It is important that `Phi_old` and `Phi_C` describe the
final reorthogonalized block basis; preliminary SVD maps are insufficient.

## Fock And Energy Convention

For a target density `D`, define

```text
J[D]pq = sum(rs) R[p,q,r,s] D[r,s]
K[D]pq = sum(rs) R[p,r,q,s] D[r,s].
```

For RHF, `rho` is the one-spin occupied projector:

```text
D       = Q_w' rho Q_w
Delta F = Q_w (2J[D]-K[D]) Q_w'
Delta E = 2<D,J[D]> - <D,K[D]>.
```

For UHF,

```text
Dtot          = Dalpha + Dbeta
Delta Falpha  = Q_w (J[Dtot]-K[Dalpha]) Q_w'
Delta Fbeta   = Q_w (J[Dtot]-K[Dbeta]) Q_w'
Delta E       = 1/2<Dtot,J[Dtot]>
              - 1/2<Dalpha,K[Dalpha]>
              - 1/2<Dbeta,K[Dbeta]>.
```

These factors make the existing core trace expressions reproduce the residual
energy. No backend energy callback is needed.

## Cost And Factorization Decision

The added global storage is `O(Nm + P^2)` and each block adds `O(km)`. A Fock
call projects and lifts in `O(w^2m + wm^2)` and performs an exact `O(m^4)`
target contraction. Window scratch is reused, and the accepted residual kernel
has zero steady allocation.

For the accepted Cr2 case, `N=7069`, `m=46`, `P=1081`. The original design
estimate was about `12.7 MiB` including target maps, block projections, and
window work; the measured resident correction storage was `12.548 MiB`.
The expanded target tensor alone would occupy `34.16 MiB` and two convenient
orientations about `68.3 MiB`.

The weighted symmetric-pair spectrum is signed, spanning approximately
`[-1.619e-2, 8.473e-3] Ha`, with 977 eigenvalues above `1e-12` in magnitude.
Even rank 512 left a spin-Fock error near `1.05e-6 Ha`. An exact dense pair
matrix was therefore both simpler and more suitable than a truncated factor
representation. Any future factor path needs explicit energy, Fock, and sweep
error gates.

## Accepted Evidence

Compact package tests cover a literal signed four-index RHF/UHF oracle,
energy factors, left/right recurrence, constructor rejection, and exact
zero-residual routing. The saved Cr2 acceptance established:

| Quantity | Accepted result |
|---|---:|
| Direct residual energy | `+30.961807 mHa` |
| Exchange residual energy | `-4.646784 mHa` |
| First-window Fock disagreement | `< 7e-15 Ha` |
| Expanded-oracle sweep energy | `-2086.524896228966554 Ha` |
| Compact-backend sweep energy | `-2086.524896228962007 Ha` |
| Sweep-energy difference | `4.55e-12 Ha` |
| Minimum alpha/beta occupied-subspace singular value | `> 0.9999999999999` |
| Returned/recomputed energy mismatch | `-5.75e-10 Ha` |

In the matched same-process comparison, the expanded local tensor route took
`73.1 s` and allocated `182.1 GiB`; the compact backend took `51.6 s` and
allocated `179.4 GiB`, a `21.5 s` (`29.4%`) improvement. Persistent correction
storage fell from `34.2 MiB` to `11.4 MiB` for `Q`, the pair residual, and pair
indices. An earlier isolated measurement found `0.504337 s/sweep` overhead and
`12.548 MiB` resident storage. An extra `2.21 GiB` of cumulative allocation in
that earlier run remained unexplained but was nonblocking.

The consumer-side expanded tensor, direct/exchange insertion, target-window
propagation, and residual sweep fork were consequently retired. Consumer
reference constants, checkpoint policy, and scientific diagnostics remain
outside this backend.

## Limits

This backend does not claim that a target space containing only occupied
directions resolves full occupied-to-virtual stationarity. It supplies an
exact compact interaction correction for the target it is given. Complex
orbitals, a generic interaction-composition layer, factor channels, and a
global four-index container remain outside the supported contract.

-- hfdmrg-manager@rh310l.ps.uci.edu
