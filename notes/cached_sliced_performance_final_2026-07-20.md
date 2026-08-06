# Cached Sliced Backend: Final Design And Evidence

Original program date: 2026-07-20
Compacted current account: 2026-08-05
Production commits: `9914ad3`, hardened by `2ae4cce`
Status: production supported for fixed and ragged sliced interactions

## Outcome

The cached sliced backend performs exact local RHF/UHF Fock contractions for
windows aligned to complete slices. It retains the combined interaction of
each block and two ordered channel orientations only to slices still outside
that block. Newly internal channels are folded into the block interaction and
discarded during absorption. Split-slice or otherwise unaligned windows
delegate completely to the projection backend.

The implementation changes representation and cost, not the Hamiltonian,
sweep ordering, damping, local iteration limit, convergence test, or backend
API. It does not construct a global physical four-index tensor or introduce a
generic cache framework.

## Ordered Sliced Convention

For slices `s,t` with local dimensions `d_s,d_t`, define the unsymmetrized
pair matrix

```text
M_st[(a,b),(c,d)] = V_st[a,b,c,d]
R(X)[(a,b),(i,k)] = X[a,i] X[b,k].
```

No relation between `M_st` and `M_ts` is assumed. For a block basis restricted
to slice `s` as `P_s`, its exterior channels and combined interaction are

```text
W_u     = sum(s in block) R(P_s)' M_su
Wswap_u = sum(s in block) R(P_s)' M_us'
G_B     = sum(u in block) W_u R(P_u).
```

`W` and `Wswap` remain distinct because direct and exchange contractions need
both ordered orientations even when the physical input has additional
symmetry.

The RHF and UHF fields follow

```text
G_RHF       = 2J(D)-K(D)
G_alpha,UHF = J(Dalpha+Dbeta)-K(Dalpha)
G_beta,UHF  = J(Dalpha+Dbeta)-K(Dbeta),
```

and use the ordinary core energy traces. Tests explicitly compare interaction
energies with `sum(rho .* F)` for RHF and half the two spin traces for UHF on
symmetric non-idempotent densities.

## Retained State And Exact Absorption

`SlicedCachedBlockState` owns

```text
ra, phi, G, exterior_slices, pair_offs, W, Wswap.
```

The ragged pair offsets pack all exterior slice pairs without padding. There
is no retained per-slice block basis, global `Vijkl` reconstruction, or channel
data for already absorbed slices.

Let `A` map the old retained block basis into the final stored new basis and
let `C_n` be the final new basis on each absorbed slice. With
`R_A=R(A)` and `R_n=R(C_n)`, exact absorption is

```text
W'_u     = R_A' W_u     + sum(n new) R_n' M_nu
Wswap'_u = R_A' Wswap_u + sum(n new) R_n' M_un'

G'_B = R_A' G_B R_A
     + sum(n) R_A' W_n R_n
     + sum(n) R_n' Wswap_n' R_A
     + sum(n,m) R_n' M_nm R_m.
```

The maps must describe final post-SVD reorthogonalized `phi_new`. The cached
path verifies that the final old-region rows lie in `old.phi`; failure uses the
exact direct builder rather than accepting stale coordinates. Hardening also
requires a complete one-sided exterior range.

Whole-slice windows assemble their nine ordered left/center/right sectors
from `G`, the two exterior channels, raw center tensors, and the left-right
cross tensors. A window that is not fully aligned constructs the projection
representation and delegates the entire interaction; it does not mix cached
and projected sectors.

## Cost And Lifecycle

Let `k` be block rank and `S2ext=sum_exterior d_s^2`. One block retains

```text
O(length(ra) k + k^4 + k^2 S2ext)
```

values in its basis, combined interaction, and two exterior channel matrices.
At a fixed edge of physical dimension `d`, completeness makes `k=d`, so the
combined interaction can require `d^4` values. Grown interior environments
continue to use occupied-subspace SVD compression.

At fixed left/right ranks and center width, aligned Fock work is independent
of total slice count. Window construction still forms the left-right cross
terms and therefore depends on the exterior slice-pair extent. The current
two-direction sweep retains one block state at every cut; this is the minimum
exact lifecycle without recomputing the opposite family during reversal.

For local SCF, a window begins with one Fock build and retains each
post-density Fock for the next diagonalization. Thus `n` updates make `n+1`
backend calls; four forced updates make five, not eight. Steady aligned Fock
allocation is zero.

## Engineering Evidence

Before optimization, the accepted H10 attribution had 401 four-function
slices, ranks 1--10, seven sweeps, and 5,970 absorptions. Of `135.240 s` solve
time, absorption consumed `105.132 s` (`77.74%`) and window construction
`22.846 s` (`16.89%`). Retained cache arrays occupied `3.656 GB`; every window
used the aligned cached route.

The packed-channel precursor changed one-small-GEMM-per-slice transformations
into large packed GEMMs. At rank 10, a representative old-channel transform
fell from `6.632 ms` to `5.280 ms`; associated zeroing fell from `0.408 ms` and
`11.582 MB` to `0.0445 ms` and `10.289 MB`.

The exterior-only recurrence then removed global interaction reconstruction
and internal-slice channels. Synthetic midpoint absorption results were:

| Rank | Previous time/allocation | Exterior-only time/allocation |
|---:|---:|---:|
| 1 | `0.185 ms / 0.501 MB` | `0.045 ms / 0.058 MB` |
| 10 | `9.463 ms / 15.159 MB` | `3.170 ms / 5.323 MB` |

The corresponding numerical parity error was about `1.8e-15`. These are
warmed engineering measurements, not a physical scaling law.

Fixed and ragged verifier coverage includes both directions, rank changes,
complete fixed edges, multi-slice right growth, asymmetric explicit
partitions, ordered unsymmetrized sectors, non-idempotent RHF/UHF densities,
and aligned routing. Accepted verifier disagreements are at roundoff, with
representative maxima near `5.1e-15` in Focks and `1.1e-14 Ha` in a sweep.

## Frozen Radial-Be Evidence

The retained physical acceptance fixture has packet SHA-256
`fba29651c02cde3cdde90850d90ba30a1abd120ebfa8ac8b64cf2e5a92fea629`,
dimension `468 = 52 x 9`, and `Nalpha=Nbeta=2`. Its saved oracle energy is
`-14.573351154440187 Ha`; direct recomputation from saved orbitals gives
`-14.573351154439951 Ha`.

One aligned forced-four sweep established:

| Quantity | Result |
|---|---:|
| Projection/cached energy difference | `1.17e-13 Ha` |
| Maximum occupied-projector distance | `2.05e-10` |
| Maximum physical/projection/cached Fock difference | `2.27e-13 Ha` |
| Windows / Fock calls / fallback windows | `98 / 490 / 0` |

After warmup, initialization took `0.036709 s` and allocated `72.157 MiB`;
the first forced-four cached sweep took `0.171570 s` and `193.927 MiB`; a
steady second sweep took `0.132383 s` and `193.159 MiB`. Fresh-process peak
RSS was `855.0625 MiB`. The matched aligned projection route took
`30.47646 s` versus `0.18147 s` cached, about `168x` on this fixture. These are
fixture-specific engineering results.

The 140-sweep characterization also separated cache correctness from solver
stationarity. Residuals at sweeps 60, 80, 100, 120, and 140 were
`1.1839560e-5`, `1.3887388e-6`, `6.4627856e-7`, `6.3292824e-7`, and
`6.3275583e-7`; the minimum was `6.3260414e-7` at sweep 134. At sweep 140 the
energy was within `1.62e-13 Ha` of the oracle and projector distance was
`6.88e-7`. All 140 rows used 98 windows and 98 incremental absorptions, with
zero fallback or projection windows. The largest returned-energy mismatch was
`2.20e-13 Ha`, and the largest Gram or idempotency error `4.56e-15`.

The complete historical per-sweep table was removed during documentation
compaction; Git history preserves it, and `scripts/bench_radial_be.jl`
regenerates the maintained measurements from the authenticated packet.

## Limits

The cached backend proves exact interaction recurrence and a large benefit
when full projection is expensive. It does not prove that end-of-sweep energy
convergence implies a stationary determinant. Stationarity, polishing,
history acceleration, frozen-orbital policies, wider centers, and alternate
cache lifecycles are independent solver questions. The backend remains
deliberately specific to sliced interactions and the existing five-operation
seam.

-- hfdmrg-manager@rh310l.ps.uci.edu
