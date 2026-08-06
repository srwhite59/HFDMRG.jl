# Cached-Sliced Physical-Coulomb Representation Design

Date: 2026-08-05
Base: `main@422b0e3`
Branch: `perf/cached-sliced-coulomb-symmetry-20260805`
Status: bounded private implementation complete; awaiting review; no chain run

## Outcome And Boundary

The generic cached-sliced backend must continue to represent arbitrary ordered,
unsymmetrized slice tensors with distinct `W` and `Wswap`. A separate private
physical-Coulomb backend will admit only interactions certified to satisfy

```text
V_st[a,b,c,d] = V_st[b,a,c,d]
V_st[a,b,c,d] = V_st[a,b,d,c]
V_st[a,b,c,d] = V_ts[c,d,a,b].
```

It stores symmetric retained and physical pairs once and retains only one
exterior orientation. The sweep core, backend API, truncation, damping,
convergence, window ordering, and generic backend are unchanged. History,
target-residual, frozen occupied, H12/H16/H20, and factorized H-chain work are
outside this milestone.

## Symmetry Certification

Construction checks every fixed or ragged slice block against all three
identities. For a compared block pair, use

```text
scale = max(1, maximum(abs,V_st), maximum(abs,V_ts))
tol   = 256 eps(Float64) max(d_s,d_t)^2 scale.
```

The maximum absolute defect for each identity must not exceed `tol`. These are
fixed internal tolerances; no public tolerance keyword is added. Rejection
occurs before block state exists. The compact route uses a single certified
orientation; roundoff-sized input asymmetry is therefore within the declared
constructor uncertainty rather than silently assumed to be exact.

## Weighted Canonical Pair Coordinates

For dimension `q`, enumerate canonical pairs `I=(i,k)` with `i <= k`, and set

```text
eta_I = 1  if i == k
        2  otherwise.
```

For a symmetric matrix `X`, define the orthonormal packed coordinate

```text
svec(X)_I = sqrt(eta_I) X[i,k].
```

For a pair-symmetric tensor, store

```text
Xhat[I,J] = sqrt(eta_I eta_J) X[i,k,j,l].
```

Thus `svec(J) = Xhat*svec(D)` is the ordered direct contraction without
additional pair weights. An arbitrary ordered value needed by exchange is
recovered as

```text
X[i,k,j,l] = Xhat[canon(i,k),canon(j,l)] /
             sqrt(eta_ik eta_jl).
```

For rank `k`, `K=k(k+1)/2`; for slice width `d_s`,
`D_s=d_s(d_s+1)/2`.

## State And Recurrence

The compact block state contains

```text
Ghat :: K x K
What :: K x sum(exterior D_s)
```

beside the ordinary block range and basis. There is no `Wswap`. Relative to
the generic state, exterior-channel storage changes from

```text
2 k^2 sum(exterior d_s^2)
```

to

```text
K sum(exterior D_s).
```

At `k=20,d=4`, the channel coefficient count is 16.4% of the generic count,
an 83.6% reduction before common basis and index metadata. `G` storage changes
from `k^4` to `K^2`, 27.6% of the generic count.

Define a weighted symmetric two-sided transform by unpacking one packed matrix,
applying `A'XA`, and repacking. Applying it to every exterior physical pair
updates `What` in `O(S d^2 k^3)` work with `O(k^2)` scratch.

For a compact tensor `Xhat` with left dimension `p` and right dimension `q`,
transform one matrix index pair at a time:

1. unpack each right-pair column as a symmetric `p x p` matrix;
2. apply the left two-sided transform and pack into a `K_r x K_q` stage;
3. unpack each stage row as a symmetric `q x q` matrix;
4. apply the right transform and accumulate the packed `K_r x K_r` output.

This retains the accepted `O(k^5)` equal-rank recurrence. Absorbing slice `n`
adds the old/new tensor once and its pair transpose once because the stored
orientation represents both `V_old,n` and `V_n,old`. Newly internal slice pairs
are accumulated from the certified physical tensor. No global pair map or full
four-index retained tensor is formed.

## Window Fock Contractions

Internal block direct fields use packed multiplication; exchange loops over
ordered orbital indices and reconstructs the required value from `Ghat`.
For a block `X` and physical slice `s`, one `What_s` supplies both ordered
sectors. The Fock update adds:

```text
J_X[i,k] = sum(a,b) W[i,k,a,b] D_s[a,b]
J_s[a,b] = sum(i,k) W[i,k,a,b] D_X[i,k]
K_Xs[i,b] = sum(k,a) W[i,k,a,b] D_Xs[k,a]
K_sX[a,k] = sum(b,i) W[i,k,a,b] D_sX[b,i].
```

The same operation applies to block--center and left--right sectors. RHF uses
`2J-K`; UHF uses total-density direct fields and like-spin exchange fields.
All densities may be symmetric, damped, and non-idempotent. Aligned windows
retain only reusable `O(k_max^2+k_max*d_max+d_max^2)` scratch. Split-slice
windows continue to delegate to the projection backend.

## Validation Gates

Committed tests will remain compact and numerical:

1. Fixed and ragged constructors accept exact Coulomb tensors and separately
   reject first-pair, second-pair, and pair-exchange defects.
2. Left and right repeated absorptions with changing ranks compare compact
   windows against projection and generic cached windows.
3. Fixed/ragged RHF and spin-asymmetric UHF Focks and interaction energies
   agree for arbitrary symmetric damped densities, including unequal left and
   right ranks and multi-slice centers.
4. Aligned routes have zero fallback/projection windows and zero steady Fock
   allocation; split windows preserve projection parity.
5. One compact RHF and common-/split-H UHF workflow agrees with the generic
   backend in returned energy and occupied projectors.
6. The package suite, cached/projection verifier, docs build, and diff checks
   pass.

After correctness, a synthetic `S=600,d=4,k=20` construction may report state
coefficient counts and bounded kernel timing. It is not a solver or H-chain
run. No H20 replay occurs before separate review.

## Closed Files And Budget

Implementation is limited to:

- one new private backend file;
- one include line in `src/HFDMRG.jl`;
- compact tests in `test/runtests.jl`;
- developer architecture documentation; and
- this maintained design/evidence record.

Preferred additions are at most 300 source lines and 120 test lines; hard
maximum is 500 implementation lines excluding this design checkpoint. If exact
exchange, bidirectional recurrence, fixed/ragged support, or zero-allocation
execution cannot fit, stop rather than weaken symmetry validation or modify the
generic backend.

## Implementation Evidence

The private sibling backend is implemented in one new file and is not exported.
The generic ordered cached backend file is unchanged. Constructor certification,
weighted canonical storage, one-orientation recurrence, lazy bidirectional
fields, projection fallback, and common-/split-H solve dispatch fit within the
500-line implementation cap: 376 source lines, 118 test lines, one
include line, and three developer-documentation lines.

All 650 package tests pass, including 17 compact Coulomb tests. Fixed and
ragged arbitrary damped-density RHF/UHF Focks agree with generic/projection
oracles within `2e-11`; the explicit development probes gave maximum Fock error
`8.89e-15` and repeated bidirectional recurrence error `2.66e-15`. Aligned
RHF/UHF additions allocate zero bytes steadily, and the compact workflows cover
RHF, common-H UHF, and split-H UHF. The existing cached/projection verifier is
unchanged and passes with maximum Fock difference `5.10702591327572e-15` and
maximum sweep-energy difference `1.0658141036401503e-14 Ha`. Documenter,
doctests, export checks, and links pass.

A disposable one-thread synthetic check at `S=600,d=4,k=20` expanded every
compact coefficient into the generic ordered representation. Persistent
interaction coefficients fell from `7,840,000` to `1,304,100`, an `83.366%`
reduction. Transformed channel and internal-block coefficients agreed within
`1.78e-15` and `1.55e-15`; both channel kernels allocated zero. Timings were
`0.00956/0.00874 s` for generic/compact channels and `0.000961/0.000596 s` for
generic/compact block transforms. These are bounded kernel checks, not H20 or
publication performance. No H12, H16, H20, history, target-residual, frozen,
or factorized-chain calculation was run.

-- hfdmrg-manager@rh310l.ps.uci.edu
