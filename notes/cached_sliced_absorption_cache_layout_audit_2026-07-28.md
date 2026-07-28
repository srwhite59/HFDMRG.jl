# Cached Sliced Absorption And Cache-Layout Audit

Date: 2026-07-28
Owner: `hfdmrg-manager`
Repository baseline: `main@d8633af56685cdb79b54e3870ae54cb9c8d91f96`
Status: **design/profiling checkpoint; production implementation not authorized**

## Outcome

The accepted H10 attribution identifies the right optimization target:
`105.132 s`, or `77.74%` of the `135.240 s` solve, is in block absorption.
The current recurrence is exact, but its representation does substantially
more allocation and data movement than the window contract requires:

1. `W` and `Wswap` are stored as `2S` separately allocated four-dimensional
   arrays, so an old block transformation launches two small GEMMs per slice.
2. Every new state allocates and zeros `W/Wswap` for all `S` slices.
3. Pair-basis matrices `P` are allocated for all slices, including zero
   matrices outside the block.
4. `W/Wswap` for already absorbed slices are retained and transformed. Internal
   `W` is used only to reconstruct `Vijkl`; internal `Wswap` has no post-build
   consumer.
5. `Vijkl` is rebuilt by scanning every active slice and is then permuted.

The smallest useful sequence is:

1. pack `W/Wswap` without changing their information content, giving two large
   old-state GEMMs and far fewer allocations;
2. separately introduce an exact recursive block-block update, retain only
   exterior `W/Wswap`, derive active slice views from final `phi`, and store the
   block-block interaction in pair-grouped order;
3. do not change the all-cut lifecycle in this program. One state per cut is
   the minimum exact lifecycle for the current two-direction schedule without
   recomputation.

The second step is the material time-and-memory opportunity. It remains within
`src/backends/sliced_basis_cached.jl`; no core, backend API, rank, damping,
convergence, ordering, or center-width change is needed.

## Accepted Evidence And Local Boundary

The paper-side accepted attribution records:

```text
slices / width                    401 / 4
retained ranks                    1-10
sweeps / absorptions              7 / 5,970
solve                             135.239535984 s
absorption                        105.131551232 s
window construction               22.846152901 s
Fock contraction                   2.097543025 s
eigensolves                        2.177754158 s
retained cache arrays          3,656,104,064 B
fallback / projection windows      0 / 0
```

The raw transaction path recorded by the paper result is absent on this host.
No H10 timing was reconstructed or inferred from an unavailable packet.
Instead, the authorized local profiling used `S=401`, `d=4`, ranks `1:10`,
one Julia thread, one BLAS thread, and unsymmetrized random tensors. It is an
engineering surrogate, not publication timing.

Scratch sources:

```text
~/dmrgtmp/hfdmrg_cached_layout_audit_20260728/profile_absorption_layout.jl
  05895a129525dd72f255c819ae1cc62c371759710c72509c2d35fd7e7c6fd0c3

~/dmrgtmp/hfdmrg_cached_layout_audit_20260728/profile_packed_zero.jl
  5d53fee7473c4804eb182a88eed032da6eae773fb9a9964c28b87191a062a00b

~/dmrgtmp/hfdmrg_cached_layout_audit_20260728/profile_exterior_recursive.jl
  e6c47a8dd109af7f1f5c733488a31dcaaf62c924f4c7d3aef46d437b0bcc48ce
```

Host/runtime:

```text
host            rh310l
Julia           1.12.6
Julia threads   1
BLAS threads    1
```

## Unsymmetrized Convention

For physical slices `s,t`, define the pair-grouped matrix

```text
M_st[(a,b),(c,d)] = V_st[a,b,c,d].
```

No relation between `M_st` and `M_ts` is assumed. For a slice-basis matrix
`X`, define the column-major pair map

```text
R(X)[(a,b),(i,k)] = X[a,i] X[b,k].
```

For a block with active slices `A` and final stored slice rows `P_s`,

```text
W_u     = sum(s in A) R(P_s)' M_su
Wswap_u = sum(s in A) R(P_s)' M_us'
G_B     = sum(u in A) W_u R(P_u).
```

`G_B[(i,k),(j,l)]` is the pair-grouped form of the current
`Vijkl[i,j,k,l]`. Distinct `W` and `Wswap` are required. None of the proposed
steps uses exchange symmetry, slice-pair symmetry, or a transpose relation
between them.

## Current Representation And Cost

Let

```text
S       number of slices
d_s     width of slice s
N       sum_s d_s
G       sum_s d_s^2
G_A     sum of d_s^2 over active slices after absorption
G_C     sum of d_s^2 over newly absorbed slices
k0, k   old and new retained ranks
```

For uniform slices, `N=S*d` and `G=S*d^2`.

One current `SlicedCachedBlockState` owns:

```text
phi        active physical rows by k
P          S matrices totaling N*k values, including exterior zeros
W/Wswap    2*k^2*G values
Vijkl      k^4 values
```

Ignoring array headers, its persistent numerical payload is

```text
active_orbitals*k + N*k + 2*k^2*G + k^4 Float64 values.
```

The leading work of one aligned absorption is:

```text
old W/Wswap transform       2*k0^2*k^2*G
new-slice contractions      2*k^2*G_C*G
global Vijkl rebuild        k^4*G_A
```

It additionally allocates/zeros `N*k + 2*k^2*G` output values, constructs
`R(A)`, new-slice pair maps, and `G_A*k^2` active pair-map values, and
permutes `k^4` values.

At fixed `d,k`, a chain of `O(S)` absorptions therefore has `O(S^2)` work
from all three leading terms. Retaining one state per cut costs
`O(k^2*d^2*S^2 + k*d*S^2 + k^4*S)` values, with the full-slice `W/Wswap`
coefficient twice the exterior-only coefficient derived below.

## Synthetic Cost Attribution

The stage profiler used a midpoint block with 200 active slices. Times are
medians of five samples and the table is a sum of separately measured stages;
it deliberately excludes core work and is not an end-to-end benchmark.

At rank 10:

| Stage | Time | Share of stage sum | Allocation |
|---|---:|---:|---:|
| Build full `P` | `0.0295 ms` | `0.31%` | `0.286 MB` |
| Allocate/zero vector `W/Wswap` | `0.4049 ms` | `4.28%` | `11.582 MB` |
| Transform old `W/Wswap` | `6.6322 ms` | `70.08%` | `0.077 MB` |
| Contract the new slice | `0.9441 ms` | `9.98%` | `0.077 MB` |
| Active pair maps | `0.1801 ms` | `1.90%` | `2.885 MB` |
| Rebuild grouped block interaction | `1.2509 ms` | `13.22%` | `0.010 MB` |
| Pair maps, permutation, `phi` copy | `0.0215 ms` | `0.23%` | `0.245 MB` |
| **Stage sum** | **`9.4632 ms`** | **100%** | **`15.159 MB`** |

Rank changes which operation dominates. At rank 1 the new-slice contraction
was `46.63%`, old-`W` transformation `20.77%`, full-`P` construction `9.79%`,
and the block-interaction rebuild `9.41%`. The separately sampled stage sum
rose from `0.185 ms / 0.501 MB` at rank 1 to
`9.463 ms / 15.159 MB` at rank 10.

The active-slice scan is directly visible at rank 10:

| Active slices | Active pair maps | Block rebuild |
|---:|---:|---:|
| 1 | `0.0030 ms` | `0.0092 ms` |
| 50 | `0.0462 ms` | `0.3238 ms` |
| 200 | `0.1768 ms` | `1.2899 ms` |
| 400 | `0.3768 ms` | `2.6131 ms` |

These measurements explain why H10 differs sharply from the earlier low-rank
Be fixture, but the unavailable H10 rank-by-cut distribution prevents an
honest conversion of these fractions into accepted H10 seconds.

## Candidate A: Packed `W/Wswap`

Store each orientation as one matrix:

```text
Wpack     :: Matrix{Float64}  # k^2 by G
Wswappack :: Matrix{Float64}  # k^2 by G
```

Slice `s` occupies a column range of length `d_s^2`. The ragged form uses
pair-width prefix offsets; it does not require padding. The exact old-state
updates become:

```text
Wpack_new     = R(A)' * Wpack_old
Wswappack_new = R(A)' * Wswappack_old.
```

This changes two loops containing `S` small GEMMs into two GEMMs. It does not
combine the two physical orientations and does not copy or reformat the global
interaction.

### Measured engineering effect

At the synthetic midpoint, packed and slice-vector results agreed exactly.
For fixed ranks 1 through 10, the old-transform speedup ranged from about
`1.13x` to `2.23x`; rank 10 changed from `6.632 ms` to `5.280 ms`.
Rank-changing controls were also exact:

```text
k0 -> k     vector          packed
1  -> 2     0.0676 ms       0.0311 ms
4  -> 7     0.6034 ms       0.4786 ms
10 -> 6     2.3953 ms       1.8555 ms
```

Packing also changes zeroing from `2S` arrays to two arrays. At rank 10:

```text
vector arrays    0.4080 ms   11.582 MB allocated
packed arrays    0.0445 ms   10.289 MB allocated
```

Logical `W/Wswap` payload is unchanged. The likely benefits are a roughly
`10-22%` reduction in the measured absorption kernel across the sampled ranks,
fewer than one thousand array objects per state, lower allocator/GC pressure,
and modestly lower RSS. The accepted `3.656 GB` array-payload figure should
remain essentially unchanged.

For H10, a reasonable gate is at least a `10%` absorption-time reduction with
no solve-time regression. A larger claim requires the paper-side packet.

### Risk, scope, and budget

Risks are ragged column-offset errors, hidden allocations from slice views,
and roundoff changes from the larger GEMM shape. Center-sector access must
continue to select `W` and `Wswap` independently.

Expected files:

```text
src/backends/sliced_basis_cached.jl
test/runtests.jl
docs/backend_architecture.md
README.md
```

No `src/core.jl`, `src/backend_api.jl`, or public interface change.

Preferred incremental budget:

```text
source additions      <= 60
test additions        <= 60
docs/README additions <= 20
hard total additions  <= 150
source deletions       >= 20
```

The deleted code should include the per-slice `W/Wswap` constructors and
old-state transformation loops it replaces.

## Candidate B: Exterior-Only State And Recursive Block Interaction

Suppose the final new basis is

```text
phi_new = [phi_old*A; C]       # left
phi_new = [C; phi_old*A]       # right
```

where `A` and every new-slice row block `C_n` describe the final
post-reorthogonalization basis. Put `R_A=R(A)` and `R_n=R(C_n)`.

For an exterior slice `u`, retain exactly:

```text
W'_u     = R_A' W_u
           + sum(n in C) R_n' M_nu

Wswap'_u = R_A' Wswap_u
           + sum(n in C) R_n' M_un'.
```

The exact recursive pair-grouped block interaction is:

```text
G'_B = R_A' G_B R_A
     + sum(n in C) R_A' W_n R_n
     + sum(n in C) R_n' Wswap_n' R_A
     + sum(n,m in C) R_n' M_nm R_m.
```

The four terms are respectively old-old, old-new, new-old, and new-new.
`W_n/Wswap_n` for a newly absorbed slice are read from the old exterior state
before being discarded. The formula remains valid for multiple new slices,
rank increases or decreases, left or right growth, and ragged widths.

The proposed persistent state is:

```text
ra
phi
G_B                         # k^2 by k^2, pair-grouped
exterior W/Wswap            # 2*k^2*G_E values
```

`P` is removed. For aligned states, a slice basis is a row view into the final
stored `phi`; pair maps are built from that view when a window needs them.
Keeping `phi` retains the current final-basis eligibility check and avoids
trusting a stale preliminary SVD map.

The current `Vijkl` permutation can also disappear: local `LL/RR` sectors can
consume `G_B` through the already existing pair-grouped convention.

### Work and storage

The exterior update costs:

```text
old exterior transform      2*k0^2*k^2*G_E
new-to-exterior terms       2*k^2*G_C*G_E.
```

The recursive block update costs, using contraction orders that do not form a
larger temporary:

```text
old-old          O(k^2*k0^4 + k^4*k0^2)
old/new cross    O(k^2*k0^2*G_C + k^4*G_C) per orientation
new-new          O(k^2*G_C^2 + k^4*G_C) for a one-chunk update
```

It has no active-range factor. At fixed `d,k`, the block-interaction part of a
full chain falls from `O(k^4*d^2*S^2)` to
`O(S*(k^6 + k^4*d^2 + k^2*d^4))`. Exterior coupling work remains
`O(S^2)`, because an exact environment must still couple to every slice it
has not absorbed.

One proposed state owns:

```text
active_orbitals*k + 2*k^2*G_E + k^4 Float64 values.
```

For all cuts at fixed rank and width:

```text
current full W/Wswap       about 2*k^2*d^2*S^2 values
exterior-only W/Wswap      about   k^2*d^2*S^2 values
```

The proposal also removes `N*k` zero-padded `P` values per state. Peak
temporary memory for a midpoint absorption falls from `O(k^2*G)` output plus
all active pair maps to `O(k^2*G_E + k^4 + k^2*G_C)`. Near the first cut,
`G_E` is still almost `G`; near the last cut it approaches zero.

### Scratch validation and expected effect

The recursive formula was compared with both the current incremental state and
a direct physical rebuild using fully unsymmetrized interactions:

| Representation | Direction | Rank | Maximum disagreement |
|---|---|---:|---:|
| Fixed | left | `3 -> 3` | `3.1e-16` |
| Fixed | right | `3 -> 4` | `1.8e-15` |
| Ragged | left | `2 -> 4` | `8.9e-16` |
| Ragged | right | `3 -> 3` | `6.7e-16` |

These maxima include recursive `G_B`, direct `Vijkl`, exterior `W`, and
exterior `Wswap` comparisons.

At the synthetic midpoint:

```text
rank       current stage sum        exterior/recursive prototype
1          0.185 ms / 0.501 MB      0.045 ms / 0.058 MB
10         9.463 ms / 15.159 MB     3.170 ms / 5.323 MB
```

The sampled kernel reduction was about `60-75%`; the prototype is not a full
backend and excludes fallback construction, window work, GC, and core
bookkeeping. A conservative H10 expectation is a `35-60%` absorption-time
reduction and a `35-55%` retained-cache-payload reduction. Applied to the
accepted totals, those are design targets, not measured new H10 results.

Suggested paper-side promotion gates are:

```text
absorption time       at least 30% below d8633af
retained payload      at least 30% below 3,656,104,064 B
total solve time      no regression
fallback/projection   exactly zero
```

All accepted P1 endpoint, stationarity, Aufbau, gap, spin, energy, allocation,
and route gates remain unchanged.

### Risk, scope, and budget

The main numerical risk is an orientation error in the two cross terms.
`Wswap_n'` is essential; it cannot be replaced from `W_n`. Other risks are:

- using `Phi_old/Phi_C` before the core's final reorthogonalization;
- aliasing an old exterior buffer that is still live at another cut;
- incorrect exterior ranges on right growth;
- partial-slice fallback accidentally entering the aligned recurrence;
- changing the direct boundary/fallback builder while optimizing the aligned
  path; and
- introducing steady Fock allocation through packed slice views.

The first implementation should allocate a new immutable state on every
absorption. It should not introduce a solve-wide mutable buffer pool.
`R_A`, `R_n`, and matrix-product temporaries may be reused within one call,
but must not survive in persistent state.

Expected files are the same four files as Candidate A. Incremental budget after
Candidate A:

```text
source additions      <= 90
test additions        <= 80
docs/README additions <= 25
hard total additions  <= 220
source deletions       >= 35
```

Deletions should include the zero-padded persistent `P`, global
`_cached_vijkl` active scan, internal-slice `W/Wswap` retention, and the
four-index permutation path made obsolete by pair-grouped `G_B`.

## Pair Maps And Buffer Reuse

The current `R(A)` is used once, while each newly absorbed slice map is reused
across all slice contractions. Global `Vijkl` reconstruction then creates a
second set of active pair maps. Candidate B removes that second set entirely.

Further persistent scratch is not justified:

- all cut states coexist, so an output `W` buffer cannot be recycled from a
  state needed by the reverse pass;
- ranks change, making a fixed persistent scratch size wasteful;
- backend-global mutable scratch would add aliasing and reentrancy risk; and
- a stale pair map is a direct stale-basis correctness defect.

Small call-local matrices may be preallocated and reused among contractions
inside one absorption. Any such reuse must have identical no-reuse numerical
controls and may not increase retained state.

## Why All Cuts Are Retained

Initialization builds the right environment for every future left-to-right
window. During that pass, each consumed right state is replaced by the newly
built left state. Those left states are then required in reverse order by the
right-to-left pass, which overwrites them with new right states for the next
sweep.

Thus the existing array holds one environment per cut, not two complete
directional copies. Under the current two-direction schedule, one state per cut
is the minimum exact lifecycle without recomputation.

A rolling `O(1)`-state scheme would have to rebuild missing opposite
environments repeatedly, increasing a sweep from `O(S)` absorptions to
`O(S^2)` absorptions and changing floating evaluation history. Checkpointing
every `q` cuts trades `O(S/q)` retained states for `O(S*q)` rebuilding work.
Neither is attractive while individual states still contain removable
internal data.

No lifecycle implementation is recommended in this program. After Candidates
A and B, a separate profile may revisit the `22.846 s` window-construction
phase, especially `VLR/VRL`, but it must not be coupled to a sweep redesign.

## Required Validation For Any Later Implementation

Committed tests should remain numerical and compact:

1. Unsymmetrized fixed and ragged physical interactions.
2. Repeated left and right absorptions.
3. Fixed, increasing, and decreasing retained ranks.
4. Complete fixed-edge ranks and interior occupied-subspace ranks.
5. Direct physical-state parity for block-block, block-exterior, and both
   `W/Wswap` orientations.
6. RHF and UHF window Fock/energy parity with symmetric non-idempotent
   densities.
7. Partial-slice, off-span, and unanchored direct fallback parity.
8. Cached/projection verifier with zero fallback for aligned fixed/ragged
   chains.
9. Common-H and split-H end-to-end trajectories.
10. Zero steady allocation in the aligned Fock kernel.
11. Full package suite and `git diff --check`.

Tests should not assert private field names, packed buffer shapes, or helper
vocabulary.

The later paper-side H10 gate must use the accepted U1 start, `nblockcenter=1`,
seven sweeps, one Julia/one BLAS thread, and unchanged convergence/local-SCF
settings. It must compare against the accepted endpoint and repeat the
attribution/memory measurement. Center width 3 is outside this audit.

## Recommendation

Authorize Candidate A alone first. It is the narrowest representation change,
has exact packed-versus-vector controls, reduces allocation objects, and
provides a clean baseline for Candidate B.

If Candidate A passes its numerical gates and produces a measurable H10
absorption improvement, authorize Candidate B as a separate commit. Candidate
B has the stronger memory and scaling benefit, but its four-sector recursive
formula deserves independent review and a separate H10 acceptance run.

Do not authorize an all-cut lifecycle change, generic cache framework, core
edit, new backend API, or center-width experiment from this memo.

-- hfdmrg-manager@macmini
