# Cached Sliced-Backend Performance Status

Date: 2026-07-19
Owner: `hfdmrg-manager`
Branch: `perf/cached-sliced-20260719`
Base commit: `f097ce49f661a9a6881131a4fc77af8d00d3b9f7`
Status: **Milestone 0 complete; paused for paper-manager scientific acceptance review**

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

The implementation will use a fixed, scale-aware internal check that
`phi_old*A` reconstructs the final old physical rows. Failure routes to the
direct builder. This is an internal safety check, not a public tolerance.

Define a pair map

```text
R(X)[(a,b),(i,k)] = X[a,i] X[b,k].
```

With `W_u` viewed as a matrix whose rows are retained pairs `(i,k)` and whose
columns are physical pairs `(a,b)`, the exact update is

```text
Wnew_u
  = R(A)' * Wold_u
  + sum(s in newly absorbed slices) R(C_s)' * pairmatrix(Vsu)

Wswapnew_u
  = R(A)' * Wswapold_u
  + sum(s in newly absorbed slices) R(C_s)' * pairmatrix(Vus)'.
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

Here complete coverage means `newra` is exactly the disjoint union of `oldra`
and `cra`, with no gap or overlap, in old-then-center order for left growth and
center-then-old order for right growth. Both input ranges must be unions of
complete layout slices. Any failed condition routes to direct rebuilding.

The direct builder will itself stop scanning known-zero global source slices:
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
`blocksize` would change existing numeric behavior. Milestone 1 will instead
add an explicit opt-in core keyword accepting a `SliceLayout` as the physical
block partition, provisionally:

```text
block_partition = nothing    # exact legacy getblocksizes path
block_partition = layout     # one complete physical slice per block
```

The core will validate contiguous `1:N` coverage, complete layout boundaries,
nonempty blocks, more than four blocks, and
`nblocks >= nblockcenter + 4`. It will not require an end slice to contain the
global occupied rank; the occupied restriction may have a smaller exact rank.
Both common-H and split-H cores will consume the same validated ranges.
Projection and cached backends can therefore be compared with identical
windows. No new exported name or backend lifecycle method is needed.

For radial Be, `block_partition=layout` gives 52 nine-orbital blocks and, with
`nblockcenter=1`, 98 moving windows. The ordinary integer-`blocksize` path will
continue to call the unchanged `getblocksizes` method and must produce
identical ranges and solver results.

The exact keyword spelling and this anticipated, bounded public solver
addition are part of the Milestone 0 review gate. If review rejects a public
keyword, implementation must stop and choose another explicit opt-in seam; it
must not silently change integer behavior.

## Cost Model And Local Baseline

Let

```text
G  = sum(all slices s) d_s^2
A2 = sum(active block slices s) d_s^2
C2 = sum(newly absorbed slices s) d_s^2.
```

The current direct builder's leading contractions are approximately

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

Current direct `_build_cached_state` timing for a full-support block:

| Slices | Wall time | Allocated bytes | State bytes |
|---:|---:|---:|---:|
| 8 | `0.071169 s` | `210,656` | `174,432` |
| 16 | `0.285134 s` | `418,976` | `346,528` |
| 32 | `1.149779 s` | `835,616` | `690,720` |
| 52 | `3.042623 s` | `1,356,416` | `1,120,960` |

A 52-slice build with only one physically active slice still took
`3.047242 s`, confirming that the present builder scans the global slice
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

The unchanged package baseline was also rerun with
`julia --project=. -e 'using Pkg; Pkg.test()'`; all 73 tests passed.

## Numerical Derivation Check

Before any source edit, the hybrid equations were evaluated in memory against
the current direct raw rebuild with unsymmetrized random interactions and
changing retained ranks:

| Case | max `W` error | max `Wswap` error | max `Vijkl` error |
|---|---:|---:|---:|
| Fixed left | `1.78e-15` | `1.11e-15` | `1.11e-15` |
| Fixed right | `1.33e-15` | `1.33e-15` | `8.88e-16` |
| Ragged left | `8.88e-16` | `8.88e-16` | `3.33e-16` |
| Ragged right | `6.66e-16` | `8.88e-16` | `8.88e-16` |

No integral symmetry was used. These are design checks only; committed tests
will exercise numerical Fock behavior rather than private field names or
buffer shapes.

The fixed check used `MersenneTwister(719)`, `dims=[2,2,2,2]`,
`k_old=2`, and `k_new=3`. It drew `V6` first, then QR-orthonormal old and final
bases for left ranges `old=1:2`, `center=3:4` and right ranges `old=7:8`,
`center=5:6`. The ragged check used `MersenneTwister(720)`,
`dims=[3,2,4,3]`, and drew every `Vblocks[n][m]` in nested `n,m` order before
the bases. Its left case was `old=1:3`, `center=4:5`, ranks `2 -> 3`; its right
case was `old=10:12`, `center=6:9`, ranks `2 -> 4`. In each case the recurrence
above was evaluated with explicit pair-map matrix products and compared with
`_build_cached_state(newrange,newphi,backend)`. This records the random draw
order needed to repeat the design check; Milestone 2 will replace it with
downstream numerical regressions.

## Validation Ladder

### Milestone 1: Explicit Aligned Partition

- exact contiguous coverage and exact `SliceLayout` boundaries;
- 52 blocks and 98 windows for the representative 52-by-9 layout;
- zero split-slice windows on the opt-in route;
- unchanged outputs from the legacy integer `blocksize` partition;
- identical ranges for projection/cached and common/split-H cores; and
- compact rejection checks for invalid or insufficient partitions.

### Milestone 2: Incremental Absorption

- multiple consecutive left and right absorptions, not only one step;
- fixed and ragged slices, changing ranks up and down, and unsymmetrized raw
  interactions to expose orientation errors;
- a deliberately nontrivial final-basis correction proving that final
  `phi_new`, rather than stale `Phi_old`, governs the update;
- exact partial-slice direct-fallback parity;
- arbitrary symmetric, non-idempotent RHF and UHF Fock parity against direct
  rebuilding and projection, with scale-aware disagreement at most `1e-11`;
- `S=4,8,16,32,52` chain timing and allocation scaling, with final-state direct
  numerical parity; and
- no full raw rebuild after an eligible aligned boundary.

A practical scaling stop gate is a warmed doubling ratio no worse than `5` for
the `8 -> 16` and `16 -> 32` complete-chain pairs, distinguishing intended
quadratic behavior from the present roughly cubic chain. `S=52` is the
representative endpoint rather than a doubling-ratio point.

### Milestone 3: Window-Local Fock

- explicit signed four-index RHF/UHF oracle remains authoritative;
- fixed/ragged aligned arbitrary-density parity and split fallback parity;
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

| Milestone | Proposed files | Maximum additions | Deletion target |
|---|---|---:|---:|
| M1 | `src/core.jl`, `test/runtests.jl`, `README.md` | `70` | `0` |
| M2 | `src/backends/sliced_basis_cached.jl`, `test/runtests.jl` | `220` | `45` |
| M3 | `src/backends/sliced_basis_cached.jl`, `test/runtests.jl`, `scripts/bench_sliced.jl`, `README.md`, `docs/backend_architecture.md` | `300` | `360` |
| M4 | `src/core.jl`, `test/runtests.jl`, `docs/backend_architecture.md` | `90` | `15` |
| **Total M1-M4** | | **`680`** | **`420`** |

Target net growth is at most `260` lines. Stop for renewed design review if a
milestone exceeds its cap, needs `backend_api.jl`, introduces another public
type/framework, or fails to shrink the cached source materially by Milestone
3.

Addition caps and the total net-growth cap are firm. Deletion figures are
targets, not correctness gates; they must not force over-compression. Counts
come from each milestone commit's `git diff --numstat`. Replacement-only churn
in this status memo is excluded.

Milestone 2 specifically replaces the scalar global raw-cache loops. Milestone
3 deletes:

- duplicated split-slice exchange fallback in favor of the projection oracle;
- aligned global `S`, `Srho`, `rho_slice`, `J`, and `tmp` scratch and scans;
- unused `rho_slice_up` and `rho_slice_dn`;
- unused center metadata; and
- misleading historical direct/exchange timing labels.

Milestone 4 consolidates duplicate pre/post-density Fock construction without
changing solver policy.

No new source file or include wiring is proposed.

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
7. **Public partition scope.** The explicit `block_partition` keyword is
   anticipated by Milestone 1 but still requires checkpoint approval. The
   integer path must not be silently changed.

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

Paper-manager review is requested to confirm that the complete-slice partition,
the projection/explicit oracle ladder, and the performance and workflow gates
are sufficient for scientific acceptance, or to report a concrete scientific
objection. Solver architecture and implementation ownership do not transfer.

After scientific acceptance review, the exact next action is Milestone 1 only:
implement and test the explicit slice-aligned partition on its own narrow
commit, update this status report, and proceed to Milestone 2 without altering
solver policy.

Until that review, no production source, tests, scripts, README, architecture
guide, manifests, or external workspace will be edited.

-- hfdmrg-manager@macmini
