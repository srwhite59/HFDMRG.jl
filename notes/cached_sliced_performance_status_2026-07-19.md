# Cached Sliced-Backend Performance Status

Date: 2026-07-19
Owner: `hfdmrg-manager`
Branch: `perf/cached-sliced-20260719`
Base commit: `f097ce49f661a9a6881131a4fc77af8d00d3b9f7`
M3 implementation commit: `ce20a73`
M4 implementation commit: `1c2a837`
M5 runner commit: `6054a2d`
Status: **Milestone 5 executed but not accepted on macmini; held for paper-manager review after frozen convergence and branch gates failed**

Solver architecture, implementation, and line-budget ownership remain with
`hfdmrg-manager`. `hfdmrg-paper-manager` reviews whether the partition,
numerical gates, and workflow evidence are sufficient for scientific
acceptance.

## Scope And Decision

Milestone 0 makes no production source change. The proposed implementation
keeps the corrected Hamiltonian and HFDMRG optimization unchanged.

The design decision is:

1. keep the existing compact block caches `P`, `W`, `Wswap`, and `Vijkl`;
2. align production block boundaries with complete physical slices through an
   explicit opt-in partition, leaving the ordinary integer-`blocksize` route
   unchanged;
3. on an aligned absorption, transform the aggregate `W` and `Wswap`, contract
   raw integrals only for newly absorbed slices, and rebuild the small
   `Vijkl` tensor from the updated `W` and `P`;
4. use a direct pair-map build for boundary initialization, numerical
   comparison, and every partial-slice fallback; and
5. use those same retained tensors for both Coulomb and exchange in an aligned
   moving window, eliminating the global sliced Fock reconstruction and scan.

No new global interaction, generic cache framework, public backend type,
factorization API, or sixth backend lifecycle operation is proposed.
`src/backend_api.jl` does not need to change.

## Governing Numerical Convention

For physical slices `s` and `t`, define

```text
Vst[a,b,c,d] = V[p,q,r,s4]
p=(s,a), r=(s,b), q=(t,c), s4=(t,d)

H2 = 1/2 sum(p,q,r,s4) V[p,q,r,s4] cdag[p] cdag[q] c[s4] c[r].
```

The required contractions are

```text
J[p,r]       = sum(q,s4) D[q,s4] V[p,q,r,s4]
Kspin[p,s4]  = sum(r,q) Dspin[r,q] V[p,q,r,s4]

G_RHF        = 2J(D) - K(D)
G_spin,UHF   = J(Dalpha + Dbeta) - K(Dspin).
```

Coulomb therefore has first/third output indices, while exchange has
first/fourth output indices. The implementation must preserve that distinction
and must not infer either raw slice-pair orientation from the other.

## Retained Block Representation

Let a retained block have dimension `k`, and let

```text
P_s[a,i] = phi[p(s,a),i]
```

be the restriction of its physical basis to slice `s`. Entries are exactly
zero outside the block's physical range.

The existing caches have the following precise meaning:

```text
W_u[i,k,a,b]
  = sum(s in block,c,d) P_s[c,i] P_s[d,k] Vsu[c,d,a,b]

Wswap_u[i,k,a,b]
  = sum(s in block,c,d) P_s[c,i] P_s[d,k] Vus[a,b,c,d]

Vijkl[i,j,k,l]
  = sum(s,t,a,b,c,d)
      P_s[a,i] P_t[c,j] P_s[b,k] P_t[d,l] Vst[a,b,c,d].
```

Equivalently,

```text
Vijkl[i,j,k,l]
  = sum(u,a,b) W_u[i,k,a,b] P_u[a,j] P_u[b,l].
```

`W` and `Wswap` are distinct ordered orientations. Keeping both is necessary
because the backend contract does not assume an additional slice-pair
symmetry. `Vijkl` is small but useful for environment-environment Coulomb and
exchange. No extra persistent Coulomb cache is needed.

The current `SlicedCachedBlockState` fields are therefore already complete:

```text
ra, phi, P, W, Wswap, Vijkl.
```

Keeping `W/Wswap` for all target slices is intentional in the first bounded
implementation. Restricting them to the outside range would roughly halve the
stack of block-cache memory, but would require extra slice maps and a fully
recursive `Vijkl`. That optimization is deferred unless measured memory makes
it necessary.

## Exact Aligned Absorption

Suppose an old block absorbs a center range made of complete slices. Let
`k_old` and `k_new` be the old and new retained ranks. The core forms an SVD
map and then reorthogonalizes the resulting physical basis. Consequently, the
authoritative old-to-new map is derived from the final `phi_new`, not copied
from the earlier `Phi_old` argument.

For left growth,

```text
phi_new = [old physical rows
           center physical rows],
```

and for right growth the two row groups are reversed. In either direction,

```text
A = phi_old' * phi_new[old physical rows, :]
C_s = phi_new[rows belonging to newly absorbed slice s, :].
```

The implementation uses a fixed, scale-aware internal check that
`phi_old*A` reconstructs the final old physical rows. Failure routes to the
direct builder. This is an internal safety check, not a public tolerance.

Define a pair map

```text
R(X)[(a,b),(i,k)] = X[a,i] X[b,k].
M(Vst)[(a,b),(c,d)] = Vst[a,b,c,d].
```

Both definitions use Julia's column-major pair ordering: the first member of
each pair varies fastest. `M` uses the same column-major pair ordering as
`_pair_map_mat!`, which builds `R(P)`. With retained pairs as rows, the ordered
swapped update uses exactly `M(Vus)'`, not `M(Vsu)` or an inferred integral
symmetry.

With `W_u` viewed as a matrix whose rows are retained pairs `(i,k)` and whose
columns are physical pairs `(a,b)`, the exact update is

```text
Wnew_u
  = R(A)' * Wold_u
  + sum(s in newly absorbed slices) R(C_s)' * M(Vsu)

Wswapnew_u
  = R(A)' * Wswapold_u
  + sum(s in newly absorbed slices) R(C_s)' * M(Vus)'.
```

Build `Pnew` directly from final `phi_new`. Then form

```text
Igrouped[(i,k),(j,l)]
  = sum(u in new block) Wnew_u * R(Pnew_u)
```

and explicitly permute pair order `(i,k,j,l)` into the existing storage order
`Vijkl[i,j,k,l]`.

The same recurrence applies to left and right growth; only the placement of
the old and center physical rows changes. It supports changing retained ranks.

### Eligibility And Fallback

The recurrence is exact only when the old block and absorbed center occupy
disjoint complete slices. If one slice is split between them, the pair map of
their sum contains old-center cross products that are absent from the formulas
above.

Incremental eligibility will therefore be determined only from exact
`SliceLayout` boundaries, correct side adjacency, and complete range coverage.
It will never be inferred from coefficient magnitudes. Boundary initialization
and any partial-slice absorption use the exact direct builder.

The old block must also be anchored to its declared side: a left block starts
at orbital `1`, and a right block ends at orbital `N`. A complete-slice internal
range that omits a prefix or suffix is not eligible for the recurrence.

Here complete coverage means `newra` is exactly the disjoint union of `oldra`
and `cra`, with no gap or overlap, in old-then-center order for left growth and
center-then-old order for right growth. Both input ranges must be unions of
complete layout slices. Any failed condition routes to direct rebuilding.

The direct builder stops scanning known-zero global source slices:
it will obtain active slice indices from the physical range, form `W/Wswap`
from those sources, and derive `Vijkl` from `W/P`. It remains the independent
production fallback; committed numerical tests will compare downstream Fock
behavior against the projection and explicit four-index routes rather than
test private cache fields.

## Alternatives Considered

### Rebuild From Active Sources After Every Absorption

Skipping zero slices makes direct initialization much cheaper, but rebuilding
all active-source contributions after every absorption still revisits every
old source. For uniform slices, a chain remains cubic in slice count.

With the proposed pair-map direct builder, the leading contraction count per
state is

```text
2*k^2*A2*G + k^4*A2,
```

plus pair-map and basis work. At full Be support this is approximately
`568,788,480` terms: about `9x` below the current builder but still about `40x`
above the hybrid absorption.

This is retained only as initialization and fallback, not as the normal
absorption path.

### Fully Recursive `Vijkl`

It is possible to transform old/old `Vijkl` and add separate old/center,
center/old, and center/center sectors from `W/Wswap` and raw center integrals.
That saves only the final small `Vijkl` reconstruction.

At the representative Be dimensions the reconstruction being removed is
approximately `1.078e6` of `14.152e6` leading terms, so `7.6%` is only an upper
bound before counting the replacement old/old and mixed-sector work. A direct
pair-space estimate gives about `6%` net arithmetic saving. The recursive form
also adds four sensitive index sectors and allows `Vijkl` to drift independently
from `W/P`. It is rejected for the first implementation and may be reconsidered
only if measured Milestone 2 cache timing fails its gate.

### Component Cache Per Absorbed Slice

Keeping one transformed cache component per absorbed source would contract raw
integrals only once, but storage and repeated basis transformations would grow
with both active and outside slice counts. Aggregating those components is
exactly the chosen `W/Wswap` recurrence and is smaller.

## Window-Local Coulomb And Exchange

For an aligned left/center/right partition, both Coulomb density indices
`(q,s4)` occupy one physical slice. A physical slice belongs wholly to one
region, so Coulomb uses only within-region density blocks. Exchange may connect
two different regions.

Define `T_XY[x,z;y,w]` as the interaction with first/third indices in region
`X` and second/fourth indices in region `Y`. The existing data provide every
ordered region pair:

| Ordered pair | Stored representation |
|---|---|
| `LL` | left `Vijkl` |
| `LC` | left `W` |
| `CL` | left `Wswap` |
| `LR` | `VLR` |
| `RL` | `VRL` |
| `CC` | raw center slice tensors |
| `CR` | right `Wswap` |
| `RC` | right `W` |
| `RR` | right `Vijkl` |

The local formulas are then

```text
J_X[x,z]    += sum(Y,y,w) D_YY[y,w] T_XY[x,z;y,w]
K_XY[x,w]  += sum(z,y)   D_XY[z,y] T_XY[x,z;y,w].
```

RHF applies `2J-K`; UHF forms total-density `J` once and spin-specific `K`.
The formulas accept arbitrary symmetric, non-idempotent densities produced by
damping.

The aligned hot path needs no full-basis lift, global slice-density scan, or
global `S/Srho/rho_slice/J/tmp` scratch. A split window will delegate to the
existing projection window instead of retaining a second copied global
implementation inside the cached backend. Specifically, it will construct the
projection-window representation from the cached states' physical ranges and
block bases and delegate the entire RHF or UHF Fock addition, not exchange
alone.

At fixed retained ranks and center slices, aligned Fock work is

```text
O(kL^4 + kR^4 + kL^2*kR^2
  + (kL^2+kR^2)*Gc + Gc^2),

Gc = sum(center slices s) d_s^2,
```

which is independent of the total number of slices. Window construction still
forms `VLR/VRL` in `O(kL^2*kR^2*G)` work, where
`G=sum(all s)d_s^2`; it will be measured separately rather than hidden in the
per-Fock result.

## Slice-Aligned Partition Seam

Automatically snapping every sliced call that supplies an integer
`blocksize` would change existing numeric behavior. Milestone 1 instead adds
an explicit opt-in core keyword accepting a `SliceLayout` as the physical
block partition:

```text
block_partition = nothing    # exact legacy getblocksizes path
block_partition = layout     # one complete physical slice per block
```

The core validates contiguous `1:N` coverage, complete layout boundaries,
nonempty blocks, more than four blocks, and
`nblocks >= nblockcenter + 4`. It will not require an end slice to contain the
global occupied rank; the occupied restriction may have a smaller exact rank.
Any backend may use the partition; when a backend itself owns a `SliceLayout`,
that layout must match exactly. Both common-H and split-H cores consume the
same validated ranges.
Projection and cached backends can therefore be compared with identical
windows. No new exported name or backend lifecycle method is needed.

The accepted keyword semantics are exact: `block_partition !== nothing`
overrides the numeric `blocksize`, and the structured route requires
`nblockcenter >= 1` while the legacy route remains unchanged. A supplied
partition must describe nonempty contiguous ranges covering exactly `1:N`; a
`SliceLayout` must end at `N`. For either sliced backend, its dimensions and
offsets must exactly equal the supplied layout. A mismatch throws rather than
entering a global fallback. The accepted spelling is
`block_partition=layout`; Milestone 1 adds neither another partition type nor
another backend lifecycle operation.

For radial Be, `block_partition=layout` gives 52 nine-orbital blocks and, with
`nblockcenter=1`, 98 moving windows. The ordinary integer-`blocksize` path will
continue to call the unchanged `getblocksizes` method and produces
identical ranges and solver results.

Only `nothing` and `SliceLayout` are approved values. Supporting another
partition type or widening these semantics requires renewed review; integer
behavior must not change silently.

### Milestone 1 Result

The Milestone 1 implementation changes only
`src/core.jl`, `test/runtests.jl`, and `README.md` outside this maintained
status memo, with `67` additions and `2` deletions against `612acd3`. The
original three-argument `getblocksizes` implementation is unchanged, and
`block_partition=nothing` calls it directly.

The structured selector accepts only a live, internally consistent
`SliceLayout`; requires positive dimensions, exact `1:N` coverage, at least
five blocks, `nblockcenter >= 1`, and enough blocks for the sweep schedule;
and requires exact dimensions and offsets when the backend already has a
`SliceLayout` property. It then uses exactly one `orb_range` per slice.

Committed-test evidence exercises:

- all 52 ranges `1:9` through `460:468` and the resulting 98-window schedule;
- exact legacy equality for omitted versus explicit `block_partition=nothing`;
- exact equality when structured calls use `blocksize=0` and `blocksize=999`;
- projection/cached one-sweep energy differences of `1.78e-15` on common-H
  and `1.59e-9` on split-H routes; and
- rejection of zero-center, mutated, wrong-length, schedule-insufficient, and
  backend-mismatched layouts.

Because each accepted block and center is an exact complete-slice range, the
52-slice route has no split-slice window by construction. The package test
script passes all 84 tests. Milestone 1 adds no backend method, state field,
fallback, or Be-specific range list.

## Cost Model And Local Baseline

Let

```text
G  = sum(all slices s) d_s^2
A2 = sum(active block slices s) d_s^2
C2 = sum(newly absorbed slices s) d_s^2.
```

The pre-Milestone-2 direct builder's leading contractions were approximately

```text
(k^4 + 2k^2) G^2
```

multiply-adds per block state, even when most `P_s` are zero. With a state at
each growth position, that produces cubic uniform-slice work.

The proposed hybrid absorption performs approximately

```text
2*k_old^2*k_new^2*G     transform W and Wswap
+ 2*k_new^2*C2*G        add newly active raw sources
+ k_new^4*A2            rebuild Vijkl from W/P
+ O(length(ra_new)*k)   rebuild P from final phi_new.
```

It is linear in total slice count for one single-slice absorption and
quadratic for a complete block-building chain or sweep.

For `S=52`, `d=9`, and `k_old=k_new=4`, the leading contraction terms are:

| Work component | Estimated leading terms |
|---|---:|
| Current rebuild, one state | `5,109,391,872` |
| Active-source direct rebuild, full support | `568,788,480` |
| Transform `W/Wswap` | `2,156,544` |
| Add one newly absorbed slice | `10,917,504` |
| Worst-position `Vijkl` rebuild | `1,078,272` |
| Proposed worst absorption total | `14,152,320` |

The leading-term reduction is about `361x` per worst-position absorption. It
omits pair-map construction, permutations, and product-chain overhead and is
not a wall-time prediction.

The supplied Be-sized fixed `V6` occupies `141,927,552` bytes (`135.353 MiB`)
and is referenced, not copied. At rank four, `W+Wswap` occupy `1,078,272` bytes
(`1.028 MiB`) per state. Including `P`, `phi`, and `Vijkl`, one measured state
is about `1.06 MiB`. Multiplying the full-support state size by roughly 51
states gives a conservative retained-cache estimate of about `54 MiB`, not a
measured live peak. The hybrid adds no persistent payload.

### Measurements Run For This Checkpoint

Environment:

```text
host: rh310l.ps.uci.edu
CPU: apple-m4
Julia: 1.12.6
Julia threads: 1
BLAS threads: 8
```

All measurements used `julia --project=.` and Julia `@timed`, `@elapsed`, or
`@allocated`. They used synthetic fixed slices with `d=9`, rank four, and the
current commit.

The direct-build table used `MersenneTwister(S)` for each row, normally
distributed `V6`, and QR-orthonormalized normally distributed `phi`. A compiled
`S=2` build preceded the table; each row then ran `GC.gc()` and one `@timed`
build. `Allocated bytes` is `@timed.bytes`, and `State bytes` is
`Base.summarysize(state)`. The one-active-slice check used seed `5201`, one
warmup call, `GC.gc()`, and one measured call.

Pre-Milestone-2 direct `_build_cached_state` timing for a full-support block:

| Slices | Wall time | Allocated bytes | State bytes |
|---:|---:|---:|---:|
| 8 | `0.071169 s` | `210,656` | `174,432` |
| 16 | `0.285134 s` | `418,976` | `346,528` |
| 32 | `1.149779 s` | `835,616` | `690,720` |
| 52 | `3.042623 s` | `1,356,416` | `1,120,960` |

A 52-slice build with only one physically active slice still took
`3.047242 s`, confirming that the old builder scanned the global slice
square rather than the active support.

For an aligned UHF window with one nine-orbital center slice,
`kL=kR=4`, and `w=17`, warmed current Fock behavior was:

| Total slices | Time per Fock | Bytes per Fock |
|---:|---:|---:|
| 8 | `0.074547 ms` | `3,744` |
| 16 | `0.260946 ms` | `10,144` |
| 32 | `1.155931 ms` | `35,232` |
| 52 | `2.954043 ms` | `89,632` |

Thus even an aligned cached Fock currently grows approximately with the global
slice square and is not steady-allocation-free. At 52 slices, warmed window
construction itself took `0.126 ms` and allocated `531,088` bytes; this phase
must be reported separately after the Fock rewrite.

The Fock ladder used seed `7000+S`, normally distributed `V6`, QR block bases,
and symmetric normally distributed alpha/beta densities. Each size received
one warmup Fock call followed by `GC.gc()`. Reported time is the mean of a loop
of 40 calls, clearing both output Fock matrices before every call; allocation
is a separate 40-call `@allocated` loop divided by 40. The 52-slice window
measurement used seed `5218`, one warmup construction, `GC.gc()`, and one
`@timed` construction.

These are engineering baselines, not scientific Be acceptance measurements.

## Milestone 2 Result

Milestone 2 is accepted cross-host through `6cb865d`. Outside this
maintained memo, the change against `8450983` touches only
`src/backends/sliced_basis_cached.jl` and `test/runtests.jl`, with `213`
additions and `80` deletions. It does not change the core, backend API, README,
window Fock algorithm, or interaction storage.

The direct builder now forms pair maps only for source slices touched by the
physical block, accumulates ordered `W` and `Wswap` with matrix products, and
derives `Vijkl` from grouped `W*R(P)` followed by the explicit `(1,3,2,4)`
permutation. Aligned absorption derives `A` from final `phi_new`, transforms
the old aggregate caches, adds only newly absorbed complete slices, and uses
the same finalizer. Stored `phi` is owned; the supplied sliced interaction is
still referenced rather than copied.

Eligibility requires whole-slice old and center ranges, exact side adjacency,
and the exact new outside range. The old-row reconstruction condition is the
fixed private test

```text
norm(old_rows - old_phi*A)
    <= 32*eps(Float64)*max(1,norm(old_rows))*max(size(old_rows)...).
```

Partial slices and failed reconstruction route through the active-source
direct builder. With `HFDMRG_BENCH_TIMING` enabled, the existing private
snapshot reports clearly named `cache_init_calls`, `cache_incremental_calls`,
and `cache_fallback_calls`; no public diagnostic or tolerance was added.

### Numerical Evidence

Compact downstream tests use unsymmetrized fixed and ragged interactions,
arbitrary symmetric non-idempotent RHF/UHF densities, two consecutive left
and right absorptions, ranks `2 -> 3 -> 2`, and final-basis rotations while
passing stale pre-rotation maps. Maximum scale-aware Fock disagreements were:

| Fixture | Incremental vs direct | Incremental vs projection |
|---|---:|---:|
| Fixed | `7.91e-16` | `7.91e-16` |
| Ragged | `6.35e-16` | `6.58e-16` |
| 52-by-9 downstream RHF/UHF | `1.51e-16` | `5.67e-15` |

Split-window projection parity, left and right partial-slice absorption with
nonzero old-center products, and aligned left/right off-span final bases all
pass the `1e-11` scale-aware gate. The compact route checks observe `(init,
incremental,fallback)=(8,8,0)` for accepted fixed/ragged chains and `(9,0,6)`
for partial-slice, off-span, and unanchored-block fallbacks.
The public aligned one-sweep solver route independently reports `(2,9,0)`.

At the final 52-slice right state, direct-rebuild relative errors were
`8.01e-16` for `W`, `8.14e-16` for `Wswap`, and `5.13e-16` for `Vijkl`.
No raw-integral symmetry was used.

The original macmini validation of `31a9c76` passed all 89 tests, but the
subsequent faraday run exposed the platform-sensitive aligned trajectory
described in the corrective-pass record below. That original result is not a
cross-host acceptance claim.

### Performance Evidence

The engineering fixture used host `rh310l.ps.uci.edu`, Julia `1.12.6`, one
Julia thread, eight BLAS threads, the working tree based on `8450983`, fixed
slices of dimension nine, and retained rank four. For each size, one fixed-seed
`V6` and one precomputed basis/map plan were reused for three same-process
trials and were allocated outside the timed region. Times below are medians;
allocations are identical across the three trials and exclude persistent `V6`.

| Slices | Complete initial chain | Allocation | Ratio |
|---:|---:|---:|---:|
| 4 | `0.000171 s` | `0.396 MiB` | - |
| 8 | `0.000739 s` | `1.793 MiB` | `4.31` |
| 16 | `0.003262 s` | `7.600 MiB` | `4.42` |
| 32 | `0.014238 s` | `31.264 MiB` | `4.36` |
| 52 | `0.042467 s` | `83.419 MiB` | `2.98` |

Thus the required `8 -> 16` and `16 -> 32` ratios are below `5`. On the
52-by-9 fixture:

| Stage | Median time | Three-trial range | Allocation | Retained output |
|---|---:|---:|---:|---:|
| Two initial states + 49 right absorptions | `0.042467 s` | `0.041034-0.042730 s` | `83.419 MiB` | `54.143 MiB` |
| 49 left + 49 right absorptions | `0.082025 s` | `0.081111-0.082044 s` | `161.708 MiB` | `106.176 MiB` |

Initial counts were `(2,49,0)` and sweep counts were `(0,98,0)`. The maximum
combined initial-plus-sweep allocation was `245.127 MiB`, below the `512 MiB`
gate; the time gates of `16 s` and `32 s` are also passed with wide margin.
Retained sizes are `Base.summarysize` of the benchmark outputs, while
allocation includes all new states, retained `W/Wswap/P/phi/Vijkl`, pair maps,
and temporary absorption work.

The completed log and PID record are:

```text
~/dmrgtmp/hfdmrg_cached_sliced_20260719/m2_validate.log
~/dmrgtmp/hfdmrg_cached_sliced_20260719/m2_validate.pid
```

The temporary repository benchmark script was removed after the run. These
are synthetic engineering measurements, not frozen-Be scientific acceptance.

### Paper-Review Corrective Pass

The faraday failure was diagnosed before replacing its fixture. Projection and
cached Fock matrices still agreed to `2.5e-16` in the first window, and the
smallest local occupied/virtual gap over the trajectory was `0.774 Ha`.
However, the first post-window block SVD retained a second singular value of
`1.947e-9`, only about twenty times the fixed `1e-10` rank cutoff. Nonlinear
amplification first became material in window two, at `1.01e-8 Ha` in energy
and `2.26e-8` in an occupied projector. The failure was therefore a poorly
conditioned trajectory oracle, not a cache or Fock-contraction disagreement.

The replacement uses an analytic sine basis, unit-spaced one-body spectrum,
six two-orbital slices, and a deterministic symmetric `0.005/(1+abs(p-q))`
density interaction represented in sliced form. Four local iterations are
forced. On macmini, projection/cached differences are:

| Route | Energy | Occupied-projector norm | Minimum local gap | Minimum retained sweep singular value |
|---|---:|---:|---:|---:|
| Common-H RHF | `9.93e-13 Ha` | `4.76e-11` | `1.048 Ha` | `0.0769` |
| Split-H UHF | `7.89e-17 Ha` | `1.44e-15` | `1.015 Ha` | `0.311` |

The committed cross-platform gates are `1e-10` scale-aware energy and `1e-9`
projector norm. Existing direct-cache and arbitrary-density Fock oracles retain
their machine-precision gates.

The general three-center correction is checked on a nonuniform
`[2,3,1,5,1,3,2]` partition with analytically self-consistent, delocalized
orbitals and a prescribed `2 Ha` Fock gap. Corrected common- and split-H
projector errors are at most `4.12e-15`, and independently recomputed energies
agree within `6.2e-15 Ha`. Restoring the old offset gives projector errors of
`6.92e-12` and `1.37e-9`, so the numerical regression rejects it.

The same pass adds the right off-span reconstruction oracle, proves that
complete-slice blocks missing their left or right physical anchor take the
direct route, and proves through a public density-density solve that a supplied
`SliceLayout` may reproduce the legacy ranges and result without a sliced
backend. The macmini test script and isolated `Pkg.test()` pass all 95 tests.
Faraday independently passed the same 95 assertions at `6cb865d`; the signed
acceptance record is committed at `c4e954b`.

## Milestone 3 Result

The aligned cached window now contracts the nine ordered left/center/right
sectors directly from `Vijkl`, `W`, `Wswap`, raw center tensors, and the two
cross tensors `VLR/VRL`. Pair-grouped tensors store `T_XY[x,z;y,w]`; conventional
`Vijkl` stores `T[x,y,z,w]`; and `Wswap` reverses the two grouped region pairs.
No interaction symmetry is inferred. For every ordered sector,

```text
J_X[x,z]   += sum(y,w) D_YY[y,w] T_XY[x,z;y,w]
K_XY[x,w] += sum(z,y) D_XY[z,y] T_XY[x,z;y,w].
```

RHF adds `2J(D)-K(D)`. UHF adds `J(Dalpha+Dbeta)-K(Dspin)` independently to
each spin Fock matrix. A complete, contiguous whole-slice L/C/R partition uses
this local kernel. Any split or off-span window constructs the projection
backend's window and delegates its entire interaction; there is no second
global fallback contraction.

The hot aligned kernel performs scalar contractions into caller-owned Fock
matrices and allocates zero bytes after warmup. At fixed left/right ranks and
center slices its work is independent of the total slice count. Window
assembly still forms `VLR/VRL` separately and remains dependent on the active
block spans; M3 does not claim to optimize that stage recursively.

### M3 Numerical And Workflow Evidence

The compact numerical test activates each of `LL`, `LC`, `LR`, `CL`, `CC`,
`CR`, `RL`, `RC`, and `RR` independently with distinct unsymmetrized tensors
and symmetric non-idempotent densities. It compares both cached and projection
RHF/UHF results with an explicit signed four-index oracle, recovers direct and
exchange pieces separately, and checks both-boundary split delegation.

The direct test entry point and isolated `Pkg.test()` each pass all `169`
assertions. The ordered-sector set contributes `79` assertions. The accepted
structured solve reports `(window_local,window_projection)=(6,0)`, while the
deliberate split fixture reports `(0,1)`. The supported fixed/ragged verifier
uses aligned partitions and gives maximum window and sweep disagreements of
`5.11e-15` and `3.55e-15 Ha`. The largest retained-benchmark Fock disagreement is
`4.32e-15`, below the `1e-12` gate.

### M3 Performance Evidence

The retained harness ran committed tree `ce20a73` on `rh310l.ps.uci.edu` with
fixed nine-orbital slices, rank-four left/right blocks, one center slice, one
Julia thread, eight BLAS threads, three batches, and the median of per-call
batch means. Inputs and cache plans are outside the timed Fock region.

| Route | S=8 cached | S=52 cached | S52 projection | S52 speedup | S52/S8 |
|---|---:|---:|---:|---:|---:|
| RHF | `0.011232 ms` | `0.011259 ms` | `36.5246 ms` | `3244x` | `1.002` |
| UHF | `0.018985 ms` | `0.018989 ms` | `61.4505 ms` | `3236x` | `1.000` |

Both cached Fock batches allocate exactly zero bytes. At S=52, the complete
initial cache chain takes `0.038434 s` and allocates `82.602 MiB`; the two-way
absorption sweep takes `0.074425 s` and allocates `160.074 MiB`, for
`242.676 MiB` combined. Counts are `(2,49,0)` for initialization,
`(0,98,0)` for absorption, and `(1,0)` for the measured aligned window.
Cached-window construction takes `0.0943 ms` and allocates `0.030 MiB`; it is
reported separately from Fock work.

The completed record is:

```text
~/dmrgtmp/hfdmrg_cached_sliced_20260719/m3_validate.log
~/dmrgtmp/hfdmrg_cached_sliced_20260719/m3_validate.pid
PID 28197
```

These are synthetic engineering measurements on the recorded manager host,
not frozen-Be or end-to-end paper scaling claims. M5 scientific acceptance
remains unrun.

## Milestone 4 Result

The core now builds one Fock matrix for the initial density in each moving
window. Every local microiteration diagonalizes the current Fock matrix, mixes
the resulting occupied projector into the density exactly as before, builds a
fresh post-density Fock matrix for the energy, and retains that completed
matrix for the next diagonalization. In symbols,

```text
F_s       = H_window + G(D_s)
P_s       = occupied projector from F_s
D_(s+1)   = (1-lambda_b) D_s + lambda_b P_s
F_(s+1)   = H_window + G(D_(s+1))
E_(s+1)   = existing RHF or UHF trace expression.
```

The pre-M4 loop redundantly reconstructed `F_(s+1)` before the next
diagonalization. Removing only that duplicate changes `2n` backend calls into
`n+1`; four forced updates require five calls instead of eight. Every genuinely
new density still starts from a fresh copy of the projected one-body matrix,
and the backend adds its complete interaction contribution exactly once. No
zero-interaction, delta-Fock, backend-purity, or retained-backend-state
assumption was introduced.

Damping, energy evaluation, local and sweep convergence, block transforms,
physical-basis expansion, and observer notification remain in their original
order. The same rearrangement applies to common-H RHF/UHF and genuinely split
one-body UHF.

### M4 Numerical And Call-Count Evidence

The compact committed fixture wraps real nonzero density-density and cached
sliced backends without changing their results. Forced-four common RHF, common
UHF, split UHF, and aligned cached UHF all make exactly five calls per window;
an immediate-convergence RHF window makes two calls for one update. Wrapped
and unwrapped solve results are exactly equal.

A saved pre/post comparison uses the same analytic six-slice fixture for
density-density and cached-sliced RHF, common UHF, and split UHF. All six
one-sweep energies, full orbital matrices, and occupied projectors are
bit-for-bit identical. All twelve snapshots from the corresponding two-sweep
observer histories are likewise bit-for-bit identical. Calls in every matched
forced-four window change from eight to five.

Both ordinary test entry points pass all `182` assertions, including `13` M4
assertions. The aligned fixed/ragged verifier still has maximum window and
sweep disagreements of `5.11e-15` and `3.55e-15 Ha`.

### M4 Performance Evidence

The timing fixture is sliced projection RHF with 16 nine-orbital slices,
rank-four retained blocks, one center slice, and four forced local updates.
Only delegated backend Fock time and allocation are compared; this is not a
complete-sweep speedup claim. One warmup preceded three measured trials.

| Metric | Pre-M4 | M4 `1c2a837` | Change |
|---|---:|---:|---:|
| Calls per trial | `208` | `130` | `-37.5%` |
| Calls per window | `8` | `5` | `-37.5%` |
| Median delegated Fock time | `0.736171 s` | `0.460419 s` | `-37.458%` |
| Delegated allocation | `80.212 MiB` | `50.133 MiB` | `-37.5%` |

The post-M4 three-trial timing range is `0.458700-0.462582 s`. The saved
baseline and comparison records are under:

```text
~/dmrgtmp/hfdmrg_cached_sliced_20260719/pre_m4_baseline_b3362f3.jls
~/dmrgtmp/hfdmrg_cached_sliced_20260719/pre_m4_observer_history_b3362f3.jls
~/dmrgtmp/hfdmrg_cached_sliced_20260719/pre_m4_projection_timing_b3362f3.log
~/dmrgtmp/hfdmrg_cached_sliced_20260719/post_m4_compare_1c2a837.log
~/dmrgtmp/hfdmrg_cached_sliced_20260719/post_m4_compare_1c2a837.pid
~/dmrgtmp/hfdmrg_cached_sliced_20260719/post_m4_projection_timing_1c2a837.log
~/dmrgtmp/hfdmrg_cached_sliced_20260719/post_m4_projection_timing_1c2a837.pid
```

These remain engineering measurements; the separate frozen M5 evidence follows.

## Milestone 5 Result

Milestone 5 was executed at runner commit `6054a2d` without changing
production source, tests, solver policy, the Hamiltonian, damping, iteration
counts, convergence settings, or frozen acceptance gates. The packet SHA256
was verified as
`fba29651c02cde3cdde90850d90ba30a1abd120ebfa8ac8b64cf2e5a92fea629`.
The full durable record is in
[the M5 evidence memo](cached_sliced_performance_m5_evidence_2026-07-19.md),
with all 50 convergence sweeps in
[the reduced TSV](cached_sliced_performance_m5_sweeps_2026-07-19.tsv).

The numerical and performance routes pass. One identically aligned
projection/cached forced-four sweep differs by `1.17e-13 Ha` in energy and at
most `2.05e-10` in occupied-projector norm. Representative explicit physical,
projection, and cached total Fock matrices agree within `2.27e-13`. Each
forced-four cached sweep uses exactly 98 local windows and 490 Fock builds,
with zero fallback or projection windows.

After one warmup, three same-process trials give median initialization, first
sweep, and steady second-sweep times of `0.036709`, `0.171570`, and
`0.132383 s`. The two sweep allocations are `193.927` and `193.159 MiB` by
`Base.gc_bytes`. A separate fresh-process run records `855.0625 MiB` final
peak RSS by `Sys.maxrss()`. These clear the `40 s` and `2 GiB` readiness gates.

The frozen matched-convergence run does not pass within 50 sweeps. At sweep 50
the independently recomputed energy agrees within `1.78e-15 Ha`, lies
`3.73e-10 Ha` from the oracle, and has orthonormality/idempotency error
`3.64e-15`, but the solver remains unconverged: the last energy change is
`9.64e-11 Ha`, the full-Fock residual is `3.60e-5`, the minimum oracle overlap
is `0.999999999632`, and the maximum projector distance is `3.83e-5`.
The `S2` and `Zspin` deviations are `3.19e-5` and `1.47e-4`, also above their
`1e-6` fingerprint gates.

The recent energy tail contracts smoothly by about `0.793` per sweep, which
suggests approximately ten additional sweeps to the energy cutoff. That is a
diagnosis only: the frozen 50-sweep cap was not widened and no extra sweep was
run. M5 is therefore **executed but not accepted** pending paper-manager
review.

## Validation Ladder

### Milestone 1: Explicit Aligned Partition

- exact contiguous coverage and exact `SliceLayout` boundaries;
- 52 blocks and 98 windows for the representative 52-by-9 layout;
- zero split-slice windows on the opt-in route;
- unchanged outputs from the legacy integer `blocksize` partition;
- exact legacy-result parity when a density-density solve receives an
  equivalent `SliceLayout`;
- identical ranges for projection/cached and common/split-H cores; and
- a nonuniform `nblockcenter=3` stationary-state regression in both core paths;
- compact rejection checks for invalid or insufficient partitions.

### Milestone 2: Incremental Absorption

- multiple consecutive left and right absorptions, not only one step;
- fixed and ragged slices, changing ranks up and down, and unsymmetrized raw
  interactions to expose orientation errors;
- a deliberately nontrivial final-basis correction proving that final
  `phi_new`, rather than stale `Phi_old`, governs the update;
- exact partial-slice direct-fallback parity;
- left and right partial-slice absorptions with nonzero old-center cross
  products, proving the direct fallback retains the omitted recurrence terms;
- deliberately off-span left and right final old-row blocks that fail the
  scale-aware reconstruction check and invoke direct rebuilding;
- complete-slice internal blocks that fail the required left/right anchoring
  and invoke direct rebuilding;
- an observable count proving that the accepted aligned route invokes no
  fallback;
- arbitrary symmetric, non-idempotent RHF and UHF Fock parity against direct
  rebuilding and projection, with scale-aware disagreement at most `1e-11`;
- `S=4,8,16,32,52` chain timing and allocation scaling, with final-state direct
  numerical parity; and
- no full raw rebuild after an eligible aligned boundary.

A practical scaling stop gate is a warmed doubling ratio no worse than `5` for
the `8 -> 16` and `16 -> 32` complete-chain pairs, distinguishing intended
quadratic behavior from the pre-Milestone-2 roughly cubic chain. `S=52` is the
representative endpoint rather than a doubling-ratio point.

The absolute provisional engineering gates for the 52-by-9, rank-four fixture
are:

```text
complete initial cache chain                <= 16 s
forward-plus-reverse absorption cache work  <= 32 s
cumulative cache-chain allocation           <= 512 MiB
```

The allocation measurement excludes the persistent input `V6` and includes
all newly retained block states, retained `W` and `Wswap`, transformed/additive
cache arrays, rebuilt `P/phi/Vijkl`, pair maps, and temporary absorption work
created inside the measured chain. A missed absolute gate requires a measured
explanation and renewed performance review, not a solver-policy change.

### Milestone 3: Window-Local Fock

- explicit signed four-index RHF/UHF oracle remains authoritative;
- separately isolate all nine ordered `LL`, `LC`, `LR`, `CL`, `CC`, `CR`,
  `RL`, `RC`, and `RR` sectors using distinct unsymmetrized `Vst/Vts` and
  symmetric non-idempotent RHF/UHF densities; in particular, exercise
  `LC/CL`, `LR/RL`, and `CR/RC` independently so opposite orientation errors
  cannot cancel;
- fixed/ragged aligned arbitrary-density parity and split fallback parity;
- delegate the entire split-window RHF/UHF interaction to projection;
- frozen workflow window error at most `1e-12` where summation order permits;
- warmed aligned Fock at least `20x` faster than projection;
- `time_Fock(S=52) / time_Fock(S=8) <= 1.5` at fixed ranks and center;
- zero steady backend allocation in the aligned Fock method;
- zero fallback windows in the accepted structured route; and
- separate timing/allocation for initialization, absorption, window assembly,
  and Fock evaluation.

Because the expected local kernel is sub-millisecond, the speed and scaling
gates will use warmed batches large enough to suppress single-call noise and
report the batch count and statistic.

### Milestone 4: Duplicate Fock Reuse

- `n+1` backend Fock calls for `n` local microiterations, including exactly
  five calls for four iterations;
- common-H and split-H RHF/UHF coverage;
- density-density and sliced trajectory, energy, and occupied-subspace parity;
- unchanged damping and convergence policy; and
- about `35%` or better backend-time reduction for forced-four local
  iterations when Fock work dominates.

### Milestone 5: End-To-End Acceptance

Use the execution contract's full scientific, timing, allocation, memory, and
trajectory ladder. No acceptance threshold may be met by changing the
Hamiltonian, local iteration count, damping, or convergence policy.
The provisional readiness gates are a warmed forced-four Be sweep below
`40 s` and cumulative allocation below `2 GiB`.

## File And Line Budget

The status memo is maintained by replacement. The following caps cover added
source, tests, scripts, and user/architecture documentation:

| Milestone | Proposed files | Addition accounting | Deletion accounting |
|---|---|---:|---:|
| M1 | `src/core.jl`, `test/runtests.jl`, `README.md` | `67` actual (`70` cap) | `2` actual (`0` target) |
| M2 | `src/backends/sliced_basis_cached.jl`, `test/runtests.jl` | `213` actual (`220` cap) | `80` actual (`45` target) |
| M2 correction | `src/core.jl`, `src/backends/sliced_basis_cached.jl`, `test/runtests.jl`, `README.md` | `97` actual | `45` actual (`0` target) |
| M3 | Cached source/tests, both sliced scripts, cached example, user/architecture docs | `300` actual (`300` cap) | `829` actual (`360` target) |
| M4 | `src/core.jl`, `test/runtests.jl`, `docs/backend_architecture.md` | `86` actual (`90` amended cap) | `20` actual (`15` target) |
| M5 validation runner | `scripts/bench_radial_be.jl` | `275` actual (`220` preferred) | `0` actual |
| **M3 + M4** | | **`386` actual after M4 amendment** | **`849` actual** |
| **Program total through M5 runner** | | **`1038` actual** | **`976` actual** |

Target net growth is at most `260` lines. Stop for renewed design review if a
milestone exceeds its cap, needs `backend_api.jl`, introduces another public
type/framework, or fails to shrink the cached source materially by Milestone
3.

Addition caps and the total net-growth cap were firm through M3. The M3
paper-review note supersedes only the previous three-addition M4 remainder: M4
again has its original `90`-addition cap, retains the `15`-deletion target, and
may not increase the production `src/core.jl` line count. Counts come from each
milestone commit's `git diff --numstat`. Replacement-only churn in this status
memo, compact validation evidence, and required paper-review governance memos
are excluded.
The corrective implementation uses `97` additions and `45` deletions. The ten
unused M1/M2 additions reduce its debit to `87`, leaving a combined `303` for
M3 and M4 without widening their individual ceilings or the program total.
M3 uses `300` additions and deletes `829` lines. M4 adds `86` and deletes `20`;
`src/core.jl` shrinks from `866` to `863` lines. Across milestone commits
through M4, counted files have `763` additions and `976` deletions, for a
substantial net shrinkage of `213` lines.

The retained M5 validation runner adds `275` lines and deletes none, exceeding
its preferred `220`-line cap by `55`. The cap was preferred rather than a hard
stop; the runner adds no source, public API, test, framework, or machine path.
Including it, the counted program has `1038` additions and `976` deletions,
for net growth of `62` lines, still below the `260`-line program target.

Milestone 2 specifically replaces the scalar global raw-cache loops. Milestone
3 deleted:

- duplicated split-slice exchange fallback in favor of the projection oracle;
- aligned global `S`, `Srho`, `rho_slice`, `J`, and `tmp` scratch and scans;
- unused `rho_slice_up` and `rho_slice_dn`;
- unused center metadata; and
- misleading historical direct/exchange timing labels.

Milestone 4 consolidates duplicate pre/post-density Fock construction without
changing solver policy.

M3/M4 added no source file, include wiring, public type, or backend API
operation.

## Risks And Stop Conditions

1. **Pair-axis permutation.** Pair-grouped `(i,k,j,l)` must be explicitly
   permuted to current `(i,j,k,l)` storage.
2. **Post-SVD reorthogonalization.** Using passed `Phi_old/Phi_C` rather than
   final `phi_new` can violate the recurrence at the requested tolerance.
3. **Partial slices.** Applying the recurrence across a split slice omits real
   old-center pair products.
4. **Hidden integral symmetry.** Collapsing `Wswap`, `VLR`, or `VRL` would make
   an unsupported orientation assumption.
5. **Small matrix overhead.** The arithmetic reduction may expose dispatch,
   allocation, or small-BLAS overhead; the scaling ladder decides whether more
   batching is justified.
6. **Window construction.** Once Fock becomes local, `VLR/VRL` assembly may be
   the next measured bottleneck. It is not recursively optimized in the first
   design.
7. **Public partition scope.** Only `block_partition=nothing` and
   `block_partition=layout::SliceLayout` are approved. Any wider accepted type
   or semantics requires renewed review; the integer path must remain exact.

Numerical disagreement, an unanticipated API expansion, persistent unfavorable
scaling, or a solver/scientific-policy choice stops implementation as specified
by the execution contract.

## Acceptance Review And Exact Next Action

The manager design records three bounded decisions:

1. use the hybrid `W/Wswap` recurrence with `Vijkl` reconstruction rather than
   a fully recursive four-sector update;
2. use an explicit `block_partition=layout` opt-in while preserving ordinary
   integer `blocksize` exactly; and
3. delegate partial-slice cached windows to projection and enforce the
   addition/net-growth budget above.

Paper-manager and faraday accepted M3 through `22f7c07`; paper-manager accepted
M4 and authorized the frozen M5 campaign in `6e9febc`. M5 was executed at
`6054a2d`: parity, timing, allocation, RSS, route, and independent-energy gates
pass, but frozen convergence, occupied-subspace, and branch-fingerprint gates
do not pass within 50 sweeps. Solver architecture and implementation ownership
remain with `hfdmrg-manager`.

The exact next action is paper-manager review of the M5 evidence and failed
frozen gates. Do not run the estimated additional sweeps, change production or
tests, alter solver/scientific policy, or soften a gate without renewed
authority.

-- hfdmrg-manager@macmini
