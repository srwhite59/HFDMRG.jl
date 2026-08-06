# Cached-Sliced Tensor-Kernel Repair Design

Date: 2026-08-05
Base: `main@bfebe9a21b40e0e96ee41b356b85e7adfdc675dd`
Branch: `perf/cached-sliced-tensor-kernels-20260805`
Status: approved bounded implementation design; no chain calculation authorized

## Outcome And Boundary

The H20 result exposes avoidable pair-space work in the generic cached-sliced
backend. Environment blocks are already constructed correctly on the fly; no
block-family rebuilding or sweep-core change is needed. The repair is confined
to two operations behind the existing backend seam:

1. transform block tensors by their true matrix indices instead of explicitly
   forming `A kron A`; and
2. contract left--right fields lazily with each current density instead of
   constructing `VLR` and `VRL` tensors for every window.

The generic backend continues to accept unsymmetrized ordered slice tensors.
It retains distinct `W` and `Wswap`. Physical Coulomb symmetry is a separate
representation design and must never be inferred by this route.

This milestone changes no sweep ordering, environment selection, rank cutoff,
damping, convergence rule, local iteration limit, backend API, public solver
surface, or interaction convention. It runs no H12, H16, H20, or frozen-mode
calculation.

## Existing Ordered Representation

For slice pair `s,u`, define

```text
M_su[(a,b),(c,d)] = V_su[a,b,c,d]
R(X)[(a,b),(i,k)] = X[a,i] X[b,k].
```

For a block `B`, with its basis restricted to slice `s` as `P_s`, the retained
ordered exterior channels are

```text
W_u[i,k,a,b]
    = sum(s in B,c,d) P_s[c,i] P_s[d,k] V_su[c,d,a,b]

Wswap_u[i,k,a,b]
    = sum(s in B,c,d) P_s[c,i] P_s[d,k] V_us[a,b,c,d].
```

The complete block interaction is stored in pair-grouped order,

```text
G_B[(i,k),(j,l)]
    = sum(u in B,a,b) W_u[i,k,a,b] P_u[a,j] P_u[b,l].
```

No symmetry between `V_su` and `V_us`, within either physical pair, or within
either retained pair is assumed.

## A. Staged Successor Transformations

For a matrix `X` and basis map `A` define the two-sided operation

```text
T_A(X) = A' X A.
```

It is evaluated as `tmp = X*A; result = A'*tmp`, without constructing
`R(A)=A kron A`. Applied to every exterior physical pair,

```text
W'_u[:,:,a,b]     = T_A(W_u[:,:,a,b])
Wswap'_u[:,:,a,b] = T_A(Wswap_u[:,:,a,b])
```

before adding newly absorbed slices. If `C_n` is the final stored successor
basis on absorbed slice `n`, those additions are

```text
W'_u[:,:,a,b] += C_n' V_nu[:,:,a,b] C_n
Wswap'_u[:,:,a,b] += C_n' V_un[a,b,:,:] C_n.
```

The second expression preserves the swapped ordered orientation explicitly.

For a four-index tensor represented as a pair matrix `X` of size
`p^2 x q^2`, define

```text
T_{L,R}(X)[i,k,j,l]
    = sum(a,b,c,d) L[a,i] L[b,k] X[a,b,c,d] R[c,j] R[d,l].
```

Implement this without pair maps in two stages:

1. apply `T_L` to each fixed `(c,d)` matrix, producing an
   `r^2 x q^2` intermediate;
2. for each fixed new pair `(i,k)`, apply `T_R` to its `q x q` matrix.

The successor block interaction is exactly

```text
G' = T_{A,A}(G_old)
   + sum(n new) T_{A,C_n}(W_old,n)
   + sum(n new) pairtranspose(T_{A,C_n}(Wswap_old,n))
   + sum(n,m new) T_{C_n,C_m}(M_nm).
```

`pairtranspose` exchanges `(i,k)` with `(j,l)`; it does not transpose either
pair internally. These are the same four sectors as the accepted pair-map
recurrence.

For an input `p^2 x q^2` transformed to rank `r`, staged work is

```text
q^2 (p^2 r + p r^2) + r^2 (q^2 r + q r^2).
```

At `p=q=r=k`, this is `O(k^5)`, versus `O(k^6)` for two generic pair-space
GEMMs. Transforming one channel physical pair costs
`O(k_old^2 k_new + k_old k_new^2)`, versus
`O(k_old^2 k_new^2)`. Across exterior slice pairs it reduces the leading
channel transformation from `O(S d^2 k^4)` to `O(S d^2 k^3)`.

One reusable intermediate of at most `k_new^2 * max(k_old^2,d_max^2)` and
matrix workspaces of at most `max(k_old,d_max) * k_new` are sufficient during
one successor construction. They are temporary and do not enlarge retained
state.

## B. Lazy Left--Right Fields

For ordered regions `X` and `Y`, the eager cross tensor is

```text
V_XY[i,k,j,l]
    = sum(s in Y,a,b) W^X_s[i,k,a,b] P^Y_s[a,j] P^Y_s[b,l].
```

It need not exist. For a density block `D_YY`, first form the physical density
on each `Y` slice,

```text
d_s[a,b] = sum(j,l) P_s[a,j] D_YY[j,l] P_s[b,l].
```

The direct field on `X` is then

```text
J_X[i,k] += sum(s,a,b) W^X_s[i,k,a,b] d_s[a,b].
```

For one spin's cross density `D_XY`, define

```text
T_s[k,a] = sum(j) D_XY[k,j] P_s[a,j]
U_s[i,b] = sum(k,a) W^X_s[i,k,a,b] T_s[k,a]
K_XY[i,l] = sum(s,b) U_s[i,b] P_s[b,l].
```

The backend adds `-K_XY` to the `X,Y` Fock block. Calling this contraction for
`(X,Y)=(L,R)` and `(R,L)` reproduces both previously eager orientations.
RHF adds `2J_X-K_XY`; UHF adds `J_X[Dalpha+Dbeta]` to both spins and subtracts
the like-spin exchange separately.

At equal ranks, one ordered cross direction per Fock call costs
`O(S_Y d^2 k^2 + S_Y d k^2)` and uses
`O(k^2 + kd + d^2)` scratch. The current path constructs
`O(k_X^2 k_Y^2)` tensors in
`O(S_Y d^2 k_X^2 k_Y^2)` work before the first Fock call, then scans them on
every call. The lazy window stores only references, dimensions, and reusable
scratch; it owns no `VLR` or `VRL`.

Steady aligned Fock allocation must remain zero. Damped, symmetric but
non-idempotent densities remain part of the contract.

## C. Separate Physical-Coulomb Representation

A later specialized backend may be admitted only after validating, for every
slice pair and at fixed scale-aware tolerances,

```text
V_st[a,b,c,d] = V_st[b,a,c,d]
V_st[a,b,c,d] = V_st[a,b,d,c]
V_st[a,b,c,d] = V_ts[c,d,a,b].
```

Under these physical Coulomb symmetries, `Wswap_u == W_u`, and both the
retained pair `(i,k)` and exterior pair `(a,b)` are symmetric. A compact state
can therefore use canonical pairs

```text
K = k(k+1)/2,  D_s = d_s(d_s+1)/2
Wsym: K x sum(exterior D_s)
Gsym: K x K
```

with explicit weights of two for off-diagonal pair contractions. Persistent
channel storage falls from

```text
2 k^2 sum(exterior d_s^2)
```

to

```text
K sum(exterior D_s).
```

This path requires its own constructor or backend type, symmetry rejection
tests, weighted direct/exchange oracle, recurrence proof, and storage report.
It will not be implemented in the generic-kernel commit. Separation-dependent
H-chain factors such as `P_n K_|n-m| P_m'` are a still later representation
and are outside this design.

## H20 Motivation And Expected Effect

The accepted H20 attribution reports `143.08 s` process time:

| Phase | Time | Fraction |
|---|---:|---:|
| Window construction | `79.54 s` | `55.6%` |
| Absorption, including initial family construction | `42.16 s` | `29.5%` |
| Local Fock plus eigensolves | `9.50 s` | `6.6%` |

The separate `19.26 s` initialization timer overlaps absorption and is not
added. There were 1,790 absorption calls: 598 initial right-family builds and
1,192 sweep absorptions. Retained arrays grew from approximately `6.111 GiB`
at projected initialization, including `5.598 GiB` of `W/Wswap`, to
`12.032 GB` after the sweep.

This milestone directly removes the dominant eager window tensors and one
power of retained rank from recursive transformations. It does not yet reduce
the persistent two-orientation channel payload. Consequently it should improve
time substantially but is not expected to close the H20 memory problem by
itself.

## Validation And Performance Gates

Committed tests remain compact and numerical:

1. Fixed and ragged unsymmetrized tensors, both absorption directions, multiple
   absorbed slices, and changing old/new ranks: compare `G`, `W`, and `Wswap`
   with direct initialization.
2. Aligned fixed and ragged windows: compare RHF and UHF Focks and interaction
   energies with the projection backend for arbitrary symmetric damped
   densities. Existing ordered-sector tests continue to isolate all nine
   orientations.
3. Require zero steady allocation for aligned RHF and UHF Fock additions.
4. Preserve split-window projection fallback, common/split-H workflows,
   aligned routing, package tests, and the cached/projection verifier.

After correctness review, one bounded synthetic H20-shaped kernel check may
use approximately 600 slices, `d=4`, and rank 20. It compares old formulas kept
only in scratch with the new operations for numerical parity, construction
time, allocation, and temporary memory. It is not a solver run or publication
benchmark. No full H20 replay occurs until separate review.

## Closed File Set And Budget

Implementation may change only:

- `src/backends/sliced_basis_cached.jl`;
- `test/runtests.jl`;
- `docs/src/developer/backend_architecture.md`; and
- this maintained design/evidence note.

Preferred implementation addition is at most 180 lines across source, tests,
and architecture documentation, with deletions expected from eager-window and
pair-map code. Hard maximum is 260 added lines. If exact generic orientation,
zero-allocation scratch, or the numerical ladder cannot fit, stop for review
instead of widening the backend API or copying the sweep core.

-- hfdmrg-manager@rh310l.ps.uci.edu
