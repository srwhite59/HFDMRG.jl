# HFDMRG Target-Space Interaction Residual Design

Date: 2026-07-17

Status: approved for bounded implementation

## Plain-language recommendation

Keep the existing density-density interaction exactly as it is and add one
small, signed correction that acts only in a fixed target orbital space.  The
caller supplies the target orbitals as columns of `Q` and the corrected
two-electron integrals as a dense matrix over unique target-orbital pairs.
For Cr2 this is a `7069 x 46` map plus a `1081 x 1081` pair matrix, not a
four-index object in the 7069-orbital working basis.

Each HFDMRG block should carry only the target orbitals projected into that
block's retained basis.  The moving window then concatenates the left-block,
bare-center, and right-block pieces of that projection.  The residual density
is formed in the 46-orbital target, the direct and exchange corrections are
contracted there, and the resulting Fock correction is lifted back into the
window.  The existing sweep schedule, SCF loop, truncation, and one-body caches
do not change.

For RHF, where `rho` is the one-spin occupied density used by HFDMRG, the
correction is `2J - K`.  For UHF, both spins see the direct field from the
total density, while each spin sees exchange only from its own density:
`F_alpha += J[D_alpha + D_beta] - K[D_alpha]` and likewise for beta.  The
energy factors below are chosen so that the existing core energy expressions
produce the residual two-body energy without any special energy hook.

The first implementation should use the full signed pair matrix, without a
factorization.  The saved Cr2 spectrum is not sufficiently low rank in the
accuracy regime of interest, and the compact pair representation already has
modest memory and a fast allocation-free kernel.

## Concrete public and backend boundary

The first feature should be the specific combined backend

```julia
backend = HFDMRG.DensityDensityTargetResidualBackend(V, Q, residual_pair)

psiup, psidn, energy = HFDMRG.solve_hfdmrg(H, backend, psiup0, psidn0;
    kwargs...)
```

with matching RHF and split-one-body UHF forms:

```julia
solve_hfdmrg(H, backend, psiup0; kwargs...)
solve_hfdmrg(Hup, Hdn, backend, psiup0, psidn0; kwargs...)
```

`DensityDensityTargetResidualBackend` should be qualified through `HFDMRG`, as
the existing specialized backend types are, rather than added to the export
list.  `solve_hfdmrg` remains the only exported solver name.  There should be
no generic additive-backend protocol in this change: there is one live
consumer and one concrete sum, density-density plus a target residual.

The backend owns a normal `DensityDensityBackend(V)` and delegates all of its
base block, window, and Fock work to that backend before adding the target
correction.  The core and `backend_api.jl` already provide every required
operation and should not change.

Screened-Hartree reference subtraction is deliberately not part of this
contract.  For the Cr2 calculation, `-Q * J_R0 * Q'` remains a static one-body
term in the caller's `H`, and the matching reference constant remains caller
bookkeeping.  HFDMRG owns only the state-dependent two-body residual and
continues to return its normal electronic energy.

## Data contract and conventions

Let the physical working-basis size be `N`, the target rank be `m`, and
`P = m * (m + 1) / 2`.

### Target map

`Q` is an owned dense `N x m` real floating-point matrix.  Its columns are the
orthonormal target orbitals in HFDMRG's physical site ordering.  Thus `Q`
maps target coefficients into the physical basis and `Q'` maps a physical
vector into the target.  The constructor should check dimensions, finiteness,
and `Q'Q ~= I`; it should not silently orthogonalize or rotate the supplied
target.

Validation tolerances are fixed implementation details, not public keywords.
For scalar type `T`, use

```text
Q_orthogonality_tolerance = 128 * eps(T) * max(N, m)
```

for `norm(Q'Q - I, Inf)`.  This accounts for the physical contraction length
while remaining about `2.0e-10` for the Cr2 Float64 dimensions.

The first implementation should support real `T <: AbstractFloat`, with `Q`,
`V`, the residual, `H`, and the orbitals required to have compatible scalar
types.  `Float64` is the acceptance type.  Complex orbitals and Hermitian
integrals are outside this bounded feature because the current solver and
oracle use real transpose conventions.

### Canonical pair storage

The target residual uses chemists' notation

```text
R[p,q,r,s] = (p q | r s)_exact - (p q | r s)_density-density
```

and the real-orbital symmetries

```text
R[p,q,r,s] = R[q,p,r,s] = R[p,q,s,r] = R[r,s,p,q].
```

The stored matrix is `residual_pair[P,P]`, with the canonical lower-triangle
ordering

```text
for p = 1:m, q = 1:p
    pair(p,q) = p * (p - 1) / 2 + q
end
```

and

```text
R[p,q,r,s] = residual_pair[pair(p,q), pair(r,s)].
```

There is no `sqrt(2)` normalization in the stored pair matrix.  When a
symmetric density is contracted for the direct field, an off-diagonal pair
therefore has weight two.  If `w[pair(p,q)]` is one for `p == q` and two
otherwise, then

```text
jpair[a] = sum_b residual_pair[a,b] * w[b] * Dpair[b].
```

The constructor should own and symmetrize a copy only after confirming that
the input is symmetric to roundoff.  The fixed internal elementwise symmetry
tolerance is

```text
pair_symmetry_tolerance = 128 * eps(T) * max(1, maximum(abs, residual_pair)).
```

There is no public tolerance keyword.  The constructor should reject the
wrong pair size, nonfinite data, or material asymmetry.  It must not require
positive semidefiniteness: this residual is a signed difference.

The dense pair matrix is the canonical representation.  A full target tensor
is only a tiny-test oracle, and signed factor channels are a possible later
optimization with a separately measured truncation contract; neither belongs
in the first public interface.

## Block and moving-window representation

For a retained block with physical range `ra` and basis `phi_B`, the target
state is

```text
Q_B = phi_B' * Q[ra,:]                    # block_dimension x m.
```

The combined block state is conceptually

```text
DDTargetBlockState(base_density_density_state, Q_B).
```

On absorption, the existing backend API supplies the old-block and center
parts of the new isometry.  Under that API's isometry invariant, both sweep
directions use

```text
Q_B_new = Phi_old' * Q_B_old + Phi_C' * Q[cra,:].
```

This costs only small block/target matrix multiplies.  It must not recompute
`phi_new' * Q[new_ra,:]` over a growing physical range.  A validation test
should compare the recurrence with that direct expression after left and
right absorptions, so any roundoff from the core's protective reorthogonalize
step is measured rather than assumed away.

For the core's `[left block; bare center; right block]` window ordering,

```text
Qwindow = vcat(Q_L, Q[Cra,:], Q_R)        # window_dimension x m.
```

`Qwindow` need not itself be exactly orthonormal: the retained window may not
contain all of `Q`.  Its singular values are useful truncation diagnostics,
but they do not change the contraction.  This is the same projected operation
used by the CR2 oracle.

The window object should contain the delegated density-density window,
`Qwindow`, and reusable target/work buffers.  Fock calls in one window can
then have zero steady heap allocation.  A window is private to the serial
sweep; this feature adds no new thread-sharing guarantee.

## Fock and energy formulas

All matrices here are real and symmetric.  Define

```text
J_R[D]_{pq} = sum_{r,s} R[p,q,r,s] D[r,s]
K_R[D]_{pq} = sum_{r,s} R[p,r,q,s] D[r,s]
```

or, equivalently, `K_R[D]_{pr} = sum_{q,s} R[p,q,r,s] D[q,s]`.
This is the same direct/exchange index convention used in
`density_density.jl` and the sliced-backend oracle tests.  Write
`<A,B> = sum(A .* B) = tr(A'B)`.

### RHF

HFDMRG's restricted `rho = C_occ * C_occ'` is a one-spin density; the total
spatial density is `2rho`.  In one moving window,

```text
D = Qwindow' * rho * Qwindow
G = 2J_R[D] - K_R[D]
DeltaF = Qwindow * G * Qwindow'.
```

The residual interaction energy is

```text
DeltaE_RHF = 2<D, J_R[D]> - <D, K_R[D]>.
```

The core already evaluates `tr(rho * (F + H1B))`.  Because `F` contains
`DeltaF`, its residual contribution is exactly the expression above; there is
no extra factor of one half and no backend energy callback.

### UHF, including split one-body UHF

Let

```text
D_alpha = Qwindow' * rho_alpha * Qwindow
D_beta  = Qwindow' * rho_beta  * Qwindow
D_total = D_alpha + D_beta.
```

Then

```text
DeltaF_alpha = Qwindow * (J_R[D_total] - K_R[D_alpha]) * Qwindow'
DeltaF_beta  = Qwindow * (J_R[D_total] - K_R[D_beta])  * Qwindow'.
```

The corresponding residual interaction components are

```text
DeltaE_direct   =  1/2 <D_total, J_R[D_total]>
DeltaE_exchange = -1/2 <D_alpha, K_R[D_alpha]>
                  -1/2 <D_beta,  K_R[D_beta]>
DeltaE_UHF      = DeltaE_direct + DeltaE_exchange.
```

These factors match the core expression

```text
1/2 tr(rho_alpha * (F_alpha + H_alpha))
+ 1/2 tr(rho_beta  * (F_beta  + H_beta)).
```

The interaction formulas are identical for common-`H` and split-`H` UHF.
The existing `Hup === Hdn` fast-path dispatch must remain intact.

## Exact zero-residual compatibility

The backend constructor should record whether every stored residual entry is
exactly zero.  Each public `solve_hfdmrg` overload should, in that case,
delegate directly to the existing `solve_hfdmrg(H, V, ...)` overload (and its
existing split/common-`H` routing), rather than instantiate combined block and
window states.  This gives the same methods, operation order, arrays, and
energy path as the established density-density solver, so the regression can
require `==`, not merely `isapprox`.

A merely small residual is not zero and must not be dropped.  Sliced backends
are unaffected because no method involving them is added or changed.

## Scaling and Cr2 cost

Let `b` be a retained block dimension, `c` a center size, `w` the moving-window
dimension, and `r` a possible signed factor rank.

- Owned target map: `O(Nm)` storage.
- Dense canonical pair residual: `O(P^2) = O(m^4/4)` storage.
- Target block state: `O(bm)` per stored block.
- Target block absorption: `O((b + c) b_new m)` work using the supplied
  isometry, with no growing-`N` projection.
- Density projection or Fock lift per spin: `O(w^2 m + w m^2)`.
- Dense-pair direct contraction: `O(P^2)`; exchange is `O(m^4)` in total,
  with about half the output elements needed because the result is symmetric.
- A factor representation would store `O(Pr)` and contract exchange in
  `O(rm^3)`.  It is computationally attractive only for a measured effective
  rank roughly below `m`, not merely below `P`.

For Cr2, `N=7069`, `m=46`, `P=1081`, the actual block partition has 56 blocks,
center sizes 129--136, retained rank at most 48, first-window dimension 200,
and maximum window dimension at most about 225.  Float64 incremental storage
is approximately:

| Item | Size |
|---|---:|
| `Q`, `7069 x 46` | 2.48 MiB |
| residual pair matrix, `1081^2` | 8.92 MiB |
| all 56 projected block maps, upper bound `56 x 48 x 46` | 0.94 MiB |
| one window map and reusable work buffers | about 0.4 MiB |
| total target-residual overhead | about 12.7 MiB |

The pair matrix has 1,168,561 values.  The expanded `46^4` tensor has
4,477,456 values, or 34.16 MiB by itself; keeping contiguous direct and
exchange orderings would take about 68.3 MiB.  A global `7069^4` tensor is not
part of the design.

For one UHF residual Fock addition at `w=225`, the projected-density and lift
multiplies plus a symmetric dense-pair direct/exchange kernel amount to roughly
34 million floating-point operations.  There are 106 windows per Cr2 sweep
and at most eight Fock additions per window (two builds in each of four SCF
microiterations).

### Measured design microbenchmark

Using the saved Cr2 residual, `m=46`, `w=225`, two rank-24 spin densities, and
eight BLAS threads on `macmini`, a preallocated dense-pair UHF kernel gave:

```text
pair versus expanded-tensor Fock Frobenius error    2.8e-17 per spin
pair end-to-end residual Fock addition              2.406 ms
pair target contraction alone                       2.203 ms
steady allocation                                   0 bytes/call
```

A benchmark-only expanded representation with both contiguous tensor
orderings took `1.635 ms/call`, but required the roughly 68.3 MiB pair of
four-index orderings.  The compact pair kernel is the better first
time/memory contract.

At the maximum `106 x 8` calls, `2.406 ms` implies about `2.04 s/sweep` in the
residual Fock kernel.  Allowing for window and block-state work, the design
estimate is `2--4 s/sweep` incremental cost.  This must be measured in the
production workflow; it is not yet a full-solver benchmark.  For context, the
CR2 copied-sweep/full-tensor run took `64.9 s/sweep`, versus `44.3 s/sweep` for
the old density-density replay.  An acceptance target is no steady per-Fock
allocation and no more than `5 s/sweep` (about 12% of that base replay) added
on the same restart and process.

## Saved residual spectrum and factorization decision

Because off-diagonal pair coordinates are unnormalized in storage, the
meaningful symmetric-pair spectrum is that of

```text
A = sqrt(W) * residual_pair * sqrt(W),
W = Diagonal(1 for diagonal pairs, 2 for off-diagonal pairs).
```

The saved Cr2 matrix has eigenvalue range
`[-1.61902477e-2, +8.47305474e-3] Ha`.  At threshold `1e-14` it has 546
positive and 491 negative eigenvalues; 977 eigenvalues exceed `1e-12` in
absolute value.  It is emphatically signed, not a Cholesky object.

Ranks needed to retain squared Frobenius weight are:

| Retained weight | Rank |
|---|---:|
| 90% | 31 |
| 99% | 90 |
| 99.9% | 171 |
| 99.99% | 254 |
| 99.9999% | 398 |

More importantly, truncation on the saved imported determinant is not
monotone in energy because the channels are signed:

| Rank | residual-energy error | spin Fock Frobenius error |
|---:|---:|---:|
| 46 | +1130.45 microHa | `3.805e-3 Ha` |
| 128 | +338.44 microHa | `6.084e-4 Ha` |
| 256 | +82.09 microHa | `8.073e-5 Ha` |
| 384 | +17.61 microHa | `1.244e-5 Ha` |
| 512 | +0.786 microHa | `1.053e-6 Ha` |

Rank 512 is already about eleven times `m`, so an `O(rm^3)` exchange is not a
speed optimization over the exact dense-pair `O(m^4)` contraction.  No
factorization or truncation should be present in the first implementation.
Any later factor path needs explicit Fock, energy, and full-sweep error gates,
not just a spectral cutoff.

## Smallest validation ladder

1. **Tiny signed integral oracle.**  For `m=3` or 4, generate a real signed
   symmetric pair matrix and reconstruct the explicit `R[p,q,r,s]`.  Verify
   canonical pair ordering, weighted direct contraction, exchange, and
   symmetry against literal four-index loops.
2. **Window projection and RHF/UHF parity.**  Use a small orthonormal `Q`,
   random retained left/right bases, and symmetric spin densities.  Compare
   the backend window Fock matrices and the formulas above with explicit
   physical integrals transformed through `Q`.  Check RHF and UHF energies,
   including every factor of two, to `1e-11 Ha` in Float64.  Compare both
   left- and right-absorption recurrences with direct `phi'Q`.
3. **Exact zero regression.**  Run RHF, common-`H` UHF, and split-`H` UHF with
   an exactly zero residual and require orbitals and energy to be `==` the
   existing density-density results.  Existing sliced tests must remain
   unchanged and pass.
4. **Saved Cr2 fixed-state parity.**  Outside the package test suite, load the
   existing packet and compare the pair kernel with the executable full-tensor
   oracle.  Required anchors are direct `+30.961807176 mHa`, exchange
   `-4.646784384 mHa`, total `+26.315022792 mHa`, with target Fock parity at
   roughly `1e-11 Ha` or better.  Do not copy the heavy artifact into HFDMRG or
   create a new schema.
5. **One production sweep.**  Run the accepted screened Cr2 restart through
   the new backend and the copied oracle with identical controls.  Require
   electronic-energy parity within `1e-9 Ha`, occupied-subspace parity, and
   matching minimum window-`Q` capture.  Only after this passes should the
   copied residual sweep be retired.
6. **Representative performance gate.**  Measure the same Cr2 sweep in one
   process after warmup.  Report wall time, allocation count/bytes, and
   incremental resident/storage estimates.  Gate on zero steady allocation in
   the residual Fock kernel, about 12.7 MiB backend storage, and no more than
   `5 s/sweep` added relative to the same density-density restart.
7. **Accepted continuation.**  If the user wants a production science result,
   continue from the accepted restart for the agreed sweep count and compare
   final energy and diagnostics with the July 17 result.  This is a workflow
   acceptance run, not a committed unit test.

The committed tests should remain small: one numerical oracle test set plus
one zero/end-to-end workflow test.  The saved Cr2 and timing gates are
acceptance checks, not package fixtures.

## Proposed implementation slices and line budget

| File | Proposed change | Added-line budget |
|---|---|---:|
| `src/HFDMRG.jl` | include the specific backend file | 1 |
| `src/backends/density_density_target_residual.jl` | constructor/validation, combined block and window states, pair kernel, RHF/UHF hooks, three solver overloads | target 185 |
| `test/runtests.jl` | compact numerical tensor, projection, zero, and workflow checks | target 105 |
| `README.md` | public constructor, pair convention, bounded-use example | target 20 |
| **Preferred total** | no core, API, manifest, or script change | **at most 320** |

The hard maximum is 380 added lines across these four files.  If the numerical
contract cannot be protected within that limit, implementation stops for
review rather than widening the framework or line budget.  Tests protect
numerical behavior only; they do not assert private helper names, scratch
shapes, or internal vocabulary.

The public addition is exactly
`HFDMRG.DensityDensityTargetResidualBackend`; the existing
`HFDMRG.solve_hfdmrg` gains backend overloads matching the sliced-backend
shape.  Internal state and scratch names are not a compatibility contract.
There is no new package dependency, artifact format, callback framework,
factor API, or arbitrary interaction-composition layer.

Implementation should be staged as:

1. data validation and allocation-free pair oracle;
2. projected block/window state plus RHF/UHF backend hooks;
3. zero-path and solver overloads;
4. tiny tests and documentation;
5. saved Cr2 numerical and performance acceptance.

## Code retired only after parity

Once the production backend passes fixed-state and one-sweep Cr2 parity, the
CR2-local `pair_to_full`, `jk`, `add_residual_fock!`, local construction of
`Qwindow`, and the residual-specific fork of the HFDMRG sweep become
unnecessary.  Consumer-level routines that add nuclear/reference constants,
screened-Hartree bookkeeping, and final scientific diagnostics remain outside
HFDMRG.

The local runner also has an `on_iteration` checkpoint/stop callback, which is
not supplied by the current HFDMRG public API.  This residual feature should
not silently grow into a callback redesign.  Before deleting the entire CR2
runner, the user must choose whether final-state production diagnostics are
sufficient or whether a separate, independently designed sweep-observer task
is warranted.

## Audit evidence, limits, and decisions still requiring review

This audit proved algebraic agreement of the proposed pair contraction with
the expanded target tensor on a dimension-matched kernel, reproduced the
saved residual energy from signed spectral channels, measured the saved
spectrum/truncation errors, and measured a zero-allocation kernel timing.  It
also established from the live code that the existing backend API carries all
needed block/window/Fock operations.

It did not implement or prove the projected-block recurrence in a live sweep,
bitwise zero-residual routing, full RHF/UHF solver parity, or the estimated
full-sweep speedup.  Those are the implementation gates above.  It also does
not claim that an occupied-only target residual fixes the known
occupied-to-virtual stationarity error; the July 17 Cr2 result already shows
that it does not.

The approved implementation boundary is therefore narrow:

- accept the specific qualified backend name and canonical unnormalized pair
  contract;
- accept an exact dense signed pair matrix for version one;
- keep screened-reference constants and any sweep-observer work outside this
  feature;
- implement within the fixed validation rules and 380-line hard maximum, then
  return the uncommitted source change for review.

Current design-task file effect: this memo only.  No source, test, script,
README, manifest, CR2, or GaussletBases file was changed.

-- hfdmrg-manager@macmini
