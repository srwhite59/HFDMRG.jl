# Packed-native Coulomb exchange kernel design

Date: 2026-08-06

Status: design checkpoint; no production implementation in this commit.

## Decision

Replace the scalar ordered-index exchange loop in the private physical-Coulomb
cached-sliced backend with a canonical-pair scatter. The replacement will read
and decode each packed exterior coefficient once per slice and apply all of its
symmetry-related contributions to the alpha and beta exchange fields together.

The change is confined to the compact Coulomb backend. It does not alter the
generic unsymmetrized cached backend, the sweep core, backend API, cache
lifecycle, stored block bases, convergence, or any history policy.

The authenticated H20 saved window establishes the priority:

| Cross-window phase | Time per call | Fraction of cross time | Allocation |
|---|---:|---:|---:|
| UHF direct | 0.752939575 ms | 16.94% | 0 bytes |
| alpha plus beta exchange | 3.546315625 ms | 79.77% | 0 bytes |
| complete cross kernel | 4.4458948 ms | 100% | 0 bytes |

The audited ranks were `k=20` on the left, center dimension 12, and right rank
15, with 300 active exterior slices of width four.

## Stored convention

For a block rank `k` and physical slice width `d`, define canonical pairs

```text
I = pair(i,k),  1 <= i <= k <= block rank
A = pair(a,b),  1 <= a <= b <= slice width.
```

Let `s(x,y)` be one on the diagonal and `sqrt(2)` off the diagonal. The stored
channel is

```text
W_hat[I,A] = s(i,k) s(a,b) W[i,k,a,b].
```

The physical coefficient is recovered once as

```text
w = W_hat[I,A] / (s(i,k) s(a,b)).
```

No ordered `k^2*d^2` tensor is reconstructed or retained.

## Existing exchange equation

For spin `sigma`, first project the cross density through the exterior slice
basis `P_s`:

```text
M_sigma = rho_sigma[X,Y] P_s'              # k by d
```

The exchange field before the final multiplication by `P_s` is

```text
E_sigma[i,b] = sum(k,a) W[i,k,a,b] M_sigma[k,a].
```

The window Fock updates are

```text
F_sigma[X,Y] -= E_sigma P_s
F_sigma[Y,X] -= P_s' E_sigma'.
```

The current implementation evaluates the sum over ordered `(i,k,a,b)` and
calls the packed lookup and scale calculation for every term and spin.

## Canonical scatter

For each canonical coefficient with `i <= k` and `a <= b`, decode `w` once and
scatter:

```text
E_sigma[i,b] += w M_sigma[k,a]

if a != b
    E_sigma[i,a] += w M_sigma[k,b]
end

if i != k
    E_sigma[k,b] += w M_sigma[i,a]
    if a != b
        E_sigma[k,a] += w M_sigma[i,b]
    end
end
```

These are exactly the one, two, or four ordered contributions represented by
the canonical coefficient. The same decoded `w` updates both spin fields. RHF
executes only the alpha/restricted field; UHF constructs both `M_alpha` and
`M_beta`, then traverses `W_hat` once.

Pair indices advance directly in the canonical loops. Diagonal/off-diagonal
scales are selected outside the scatter statements, and the coefficient is
decoded outside both spin applications. No pair-map arrays are required in the
hot loop.

The scatter changes floating-point summation order but not the mathematical
terms. It should continue to use fused multiply-add accumulation. Numerical
acceptance is therefore projector/Fock/energy parity at fixed scale-aware
tolerances, not bitwise equality of intermediate fields.

## Work and memory

The density projection and final matrix multiplications are unchanged. The
coefficient traversal changes from

```text
UHF ordered lookup work: 2 k^2 d^2
canonical coefficient reads and decodes: K_k K_d
```

where `K_n=n(n+1)/2`. The number of physical multiply-add contributions remains
the required ordered count; the saving is repeated index normalization, scale
division, and cache reads.

For `k=20,d=4`, each slice has:

```text
current: 6,400 ordered coefficient lookups per spin
new:     2,100 unique packed coefficients total
```

The persistent `G` and `W` representation is unchanged. A fused UHF call needs
four reusable `k*d` matrices: two projected cross densities and two exchange
fields. Extending the current scratch by the beta pair adds only `O(kd)` window
storage. No call-sized allocation is permitted.

The direct kernel is explicitly deferred. Its two hand-written packed
matrix-vector contractions may later be compared with `mul!`, but that work
must not obscure exchange attribution or enlarge this change.

## Implementation boundary

Allowed implementation files:

- `src/backends/sliced_basis_cached_coulomb.jl`
- `test/runtests.jl`

Expected source changes:

1. Add alpha/beta mixed and field scratch, still `O(kd)`.
2. Replace the single-spin ordered exchange loop with restricted and fused-UHF
   canonical scatter paths.
3. Keep the final `E*P` and transpose updates unchanged.
4. Delete the superseded ordered-loop helper once all gates pass.

Preferred budget: at most 60 added source lines and 80 added test lines; hard
maximum 180 total additions before deletions. Stop for review rather than
adding a generic pair framework or persistent lookup tables.

## Validation ladder

1. **Explicit algebra oracle.** Compare compact full-window RHF and UHF Focks
   with direct physical sliced contraction for deterministic symmetric tensors,
   including non-idempotent damped densities.
2. **Layout and direction coverage.** Exercise fixed and ragged layouts, unequal
   left/right ranks, both growth directions, and sequential rank-changing
   absorption.
3. **Spin coverage.** Cover RHF, common-one-body UHF, and split-one-body UHF.
   Require the fused UHF result to match separately evaluated alpha and beta
   exchange oracles.
4. **Routing.** Require aligned compact routing with zero fallback/projection
   windows. Existing misaligned failover behavior remains unchanged.
5. **Workflow parity.** Compare compact and projection energies and occupied
   projectors through complete small sweeps. Record any changed damping decision;
   no unexplained trajectory divergence is acceptable.
6. **Regression.** Run the complete package suite, cached/projection verifier,
   and `git diff --check`.
7. **Allocation.** After warmup, the full compact UHF Fock application must
   allocate zero bytes in steady state.
8. **Authenticated saved window.** Before any HFDMRG solve, repeat the H20
   middle-window audit at one Julia and one BLAS thread. Require exact finite
   output, roundoff symmetry/composition parity, and no exchange timing
   regression. A reduction of at least 25% from the accepted 3.546315625 ms
   two-spin exchange baseline is the useful-performance gate.

The saved-window input SHA-256 is
`2a64d04123330b4e4420ba2e59434ff3842c04c305d7ab44997ab3565141fd20`.
The accepted audit runner and log SHA-256 values are respectively
`afa8c711fd56962a48905acd3c5e2a8a2726ad5fbfa36b491bc00bd32582a6bf`
and `b978c862c8182a576d9f9ab7586322dab718683075a62910e6042cc27383fbda`.

## Stop conditions and non-goals

Stop if correctness requires ordered persistent tensors, nonzero steady
allocation, a generic cache framework, core/backend-API changes, or relaxed
Coulomb certification. Do not run H12, H16, H20, frozen-occupied calculations,
or history calculations during this kernel milestone.

This design does not claim a full-sweep speedup. A complete H20 replay remains
separately gated after kernel correctness and saved-window performance review.
