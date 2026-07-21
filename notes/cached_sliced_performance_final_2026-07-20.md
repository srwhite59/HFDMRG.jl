# Cached Sliced-Backend Final Design And Evidence

Date: 2026-07-20
Owner: `hfdmrg-manager`
Program base: `f097ce49f661a9a6881131a4fc77af8d00d3b9f7`
Accepted through: `03529d72a9dcf9fa7a6401e074fc74966a8c8d39`
Status: **cache/performance program complete; accepted for promotion review**

## Outcome And Scope

The corrected sliced HFDMRG route now uses exact incremental block caches and a window-local
Coulomb/exchange kernel. On the frozen 468-orbital radial-Be fixture, an identically aligned
projected sweep taking about `30.5 s` becomes a cached sweep taking about `0.18 s`; energies,
occupied projectors, and independent full-Fock calculations agree far inside their gates. Every
accepted frozen-Be aligned cached sweep uses 98 local windows and no fallback or projection window.

This changed software organization and cost, not the Hamiltonian or optimization policy. It added
no global four-index interaction, generic cache framework, public backend type, factorization API,
or lifecycle operation, and changed no damping, iteration count, or convergence criterion. It did
not materially change the frozen occupied branch.

The original frozen 50-sweep compound gate remains failed. An authorized 140-sweep study showed
that the cache stays on the saved occupied branch while the full-Fock residual reaches a floor near
`6.33e-7`. Cache acceptance and solver stationarity are separate conclusions.

Key commits are M1 partition `8450983`; M2 absorption `31a9c76` and corrections
`a3dccf0`/`6cb865d`; M3 local Fock `ce20a73`; M4 Fock reuse `1c2a837`; M5 runner/evidence
`6054a2d`/`ead4a4b`; and M5b runner/evidence `b06388e`/`d6b1f2a`.

## Interaction And Energy Convention

For physical slices `s` and `t`, with local indices `a,b,c,d`,

```text
Vst[a,b,c,d] = V[p,q,r,s4]
p=(s,a), r=(s,b), q=(t,c), s4=(t,d)
H2 = 1/2 sum(p,q,r,s4) V[p,q,r,s4] cdag[p] cdag[q] c[s4] c[r].
```

For a one-spin density `D` and spin densities `Dalpha,Dbeta`,

```text
J[p,r]       = sum(q,s4) D[q,s4] V[p,q,r,s4]
Kspin[p,s4]  = sum(r,q) Dspin[r,q] V[p,q,r,s4]
G_RHF        = 2J(D) - K(D)
G_spin,UHF   = J(Dalpha + Dbeta) - K(Dspin).
```

Coulomb has first/third output indices and exchange has first/fourth output indices. Neither
slice-pair orientation is inferred from the other. For validation,
`D_sigma=C_sigma*C_sigma'`, `F_sigma=H_sigma+G_sigma`, and the independent energies are

```text
E_UHF = 1/2 Tr[Dalpha(Falpha + Halpha)] + 1/2 Tr[Dbeta(Fbeta + Hbeta)]
E_RHF = Tr[D(F + H)].
```

The frozen-Be route has `Halpha=Hbeta=H`.
The occupied-virtual residual is the combined norm of
`F_sigma*C_sigma - C_sigma*(C_sigma'*F_sigma*C_sigma)`.

## Partition And Retained Block State

The opt-in `block_partition=layout`, with `layout::SliceLayout`, overrides numeric `blocksize`,
requires `nblockcenter >= 1`, and supplies one nonempty orbital range per complete layout slice.
Those ranges must cover `1:N` contiguously, with at least five blocks and
`nblocks >= nblockcenter + 4`. A sliced backend's dimensions and offsets must match exactly.
`block_partition=nothing` retains legacy integer behavior. The ranges serve both core paths;
non-sliced backends may use an equivalent layout.

For a retained block basis `phi` of rank `k`, define

```text
P_s[a,i] = phi[p(s,a),i]
W_u[i,k,a,b] = sum(s in block,c,d) P_s[c,i] P_s[d,k] Vsu[c,d,a,b]
Wswap_u[i,k,a,b] = sum(s in block,c,d) P_s[c,i] P_s[d,k] Vus[a,b,c,d]
Vijkl[i,j,k,l] = sum(u,a,b) W_u[i,k,a,b] P_u[a,j] P_u[b,l].
```

The retained state is exactly `ra, phi, P, W, Wswap, Vijkl`; `W` and `Wswap` are distinct ordered
orientations, and the large sliced interaction `V6` is referenced rather than copied.

## Exact Aligned Absorption

After the core's SVD and final reorthogonalization, the final `phi_new` is
authoritative. For either sweep direction,

```text
A   = phi_old' * phi_new[old physical rows,:]
C_s = phi_new[newly absorbed slice-s rows,:].
```

Define column-major pair maps, with the first member of each pair varying
fastest:

```text
R(X)[(a,b),(i,k)]    = X[a,i] X[b,k]
M(Vst)[(a,b),(c,d)] = Vst[a,b,c,d].
```

With retained pairs as matrix rows, the exact updates are

```text
Wnew_u = R(A)'*Wold_u + sum(s new) R(C_s)'*M(Vsu)
Wswapnew_u = R(A)'*Wswapold_u + sum(s new) R(C_s)'*M(Vus)'.
```

The transpose on `M(Vus)'` is required. `Pnew` comes from final `phi_new`; grouped
`Wnew*R(Pnew)` products are permuted from `(i,k,j,l)` to stored `Vijkl[i,j,k,l]`.

Eligibility requires complete slices, exact adjacency/union, the requested new outside range to
equal the exact complement, and a side-anchored old block (`1:l` left or `r:N` right). A fixed
scale-aware private check requires `phi_old*A` to reconstruct the final old rows. Partial slices,
unanchored blocks, off-span rows, or any failure use the exact active-source direct builder, as does
boundary initialization.

Aligned complete-slice windows use the local kernel. Split or off-span windows construct the
projection representation and delegate the entire RHF/UHF interaction; no duplicate global cached
fallback algebra remains.

For one absorption, with `G=sum_s d_s^2`, `A2=sum(active slices s) d_s^2`, and
`C2=sum(new slices s) d_s^2`, leading work is

```text
2*kold^2*knew^2*G + 2*knew^2*C2*G + knew^4*A2.
```

Thus one absorption is linear in total slice count and a full chain is quadratic, not cubic. At
rank four, `W+Wswap` occupy `1.028 MiB` per 52-slice state; retained output was `54.143 MiB` for
initialization states and `106.176 MiB` for a two-way chain.

## Window-Local Coulomb And Exchange

For complete left/center/right regions, the nine ordered sectors use `Vijkl` for `LL/RR`, `W` for
`LC/RC`, `Wswap` for `CL/CR`, `VLR/VRL` for `LR/RL`, and the raw center tensor for `CC`.

Writing each grouped sector as `T_XY[x,z;y,w]`,

```text
J_X[x,z]   += sum(Y,y,w) D_YY[y,w] T_XY[x,z;y,w]
K_XY[x,w] += sum(z,y)   D_XY[z,y] T_XY[x,z;y,w].
```

All nine orientations were tested independently with unsymmetrized tensors
and symmetric non-idempotent RHF/UHF densities. At fixed retained ranks and
center size, with `Gc=sum(center slices s) d_s^2`, aligned Fock work is
independent of total slice count:

```text
O(kL^4 + kR^4 + kL^2*kR^2
  + (kL^2+kR^2)*Gc + Gc^2).
```

Window construction still forms `VLR/VRL` in `O(kL^2*kR^2*G)` work, where
`G=sum(all slices s) d_s^2`. This program did not recursively optimize that
separate stage.

## Shared-Core Fock Reuse

Each window builds an initial Fock, then after each density update builds the completed
post-density Fock for both energy evaluation and the next diagonalization. Thus `n` updates use
`n+1`, not `2n`, delegated builds; four forced updates use five calls, not eight. Every density
still gets a fresh one-body copy and one complete interaction. Damping, energy, observer ordering,
and convergence policy are unchanged.

## M1 And M2 Evidence

The 52-by-9 layout produces 52 blocks and 98 windows. Tests cover fixed/ragged slices, rank changes,
both directions, final-basis rotations, partial/off-span/unanchored fallbacks, and arbitrary
symmetric non-idempotent RHF/UHF densities. The exposed `nblockcenter >= 3` right-offset defect was
corrected in both core paths. A deterministic, well-gapped trajectory replaced a
platform-sensitive near-rank-cutoff fixture.

Maximum fixed/ragged incremental-versus-direct Fock error was `7.91e-16`;
incremental-versus-projection error was `7.91e-16`. On the 52-by-9 route the
largest downstream error was `5.67e-15`. Accepted fixed/ragged chains record
zero fallback. At pinned commit `ce20a73`, the retained synthetic cache-chain
engineering checkpoint gave:

| Stage | Median time | Allocation |
|---|---:|---:|
| Complete initial chain | `0.038434 s` | `82.602 MiB` |
| 49 left + 49 right absorptions | `0.074425 s` | `160.074 MiB` |

Combined allocation was `242.676 MiB`, excluding persistent `V6`. The earlier
accepted M2 ladder gave 8-to-16 and 16-to-32 chain-time ratios of `4.416` and
`4.365`. These are synthetic engineering results, not physical or publication
scaling evidence.
Faraday independently passed all 95 assertions with Julia 1.12.6 and ten-thread
OpenBLAS. Its one-sweep split-H projector mismatch of `6.85e-10` motivated the
documented cross-platform `1e-9` gate; further sweeps reduced it to machine
precision without exposing a cache defect.

## M3 And M4 Evidence

The M3 ordered-sector oracle recovers direct and exchange separately. The
retained fixed/ragged verifier's maximum window and sweep disagreements are
`5.11e-15` and `3.55e-15 Ha`. On macmini, warmed aligned kernels gave:

| Route | S=52 cached | S=52 projection | Speedup | S52/S8 |
|---|---:|---:|---:|---:|
| RHF | `0.011259 ms` | `36.5246 ms` | `3244x` | `1.002` |
| UHF | `0.018989 ms` | `61.4505 ms` | `3236x` | `1.000` |

Steady aligned Fock allocation is zero. Faraday independently passed all 169
assertions and the verifier, reporting `3044x/3796x` RHF/UHF speedups and
`0.979/0.991` S52/S8 ratios. These are warmed synthetic-kernel engineering
measurements, not complete-solver or physical scaling claims. M3 removed the
global density lift, slice-pair scans, duplicate split fallback, and stale
scratch; cached source shrank from 886 to 422 lines.

M4 proves five calls per forced-four window and two for one accepted update.
Faraday passed all 182 assertions and found six density-density/cached RHF,
common-UHF, and split-UHF trajectories—including twelve two-sweep observer
snapshots—bit-for-bit identical before and after M4. Delegated projection-Fock
time changed from `0.736171` to `0.460419 s`, and allocation from `80.212` to
`50.133 MiB`, both about `37.5%`; these are engineering-only delegated-work
figures. M4 added `86` and deleted `20` counted lines, while `src/core.jl`
shrank by three lines.

## Frozen-Be Acceptance

The frozen packet is identified by:

```text
SHA256: fba29651c02cde3cdde90850d90ba30a1abd120ebfa8ac8b64cf2e5a92fea629
dimension: 468 = 52 radial slices * 9 angular functions
occupations: Nalpha = Nbeta = 2
```

The macmini runs used host `rh310l.ps.uci.edu` (Apple M4 Pro), Julia 1.12.6,
one Julia thread, eight-thread ILP64 OpenBLAS, and the normal home depot. The
independent full six-index contraction used the J/K and determinant-energy
conventions above; it did not call either sliced backend.

The distant seed energy is `-13.489740341026078 Ha`. The packet records the
oracle reference energy as `-14.573351154440187 Ha`; independently recomputing
the energy of its saved oracle orbitals gives `-14.573351154439951 Ha`, a
`2.36e-13 Ha` consistency difference. The oracle residual is `1.32e-8`, with
`S2=0.13165158762762275` and `Zspin=-1.1194338511924447`.

One identically aligned forced-four sweep gives:

| Quantity | Result |
|---|---:|
| Projection/cached energy difference | `1.17e-13 Ha` |
| Maximum occupied-projector distance | `2.05e-10` |
| Maximum physical/projection/cached total-Fock disagreement | `2.27e-13 Ha` |
| Windows / Fock calls / fallback windows | `98 / 490 / 0` |

After one warmup, three same-process cached trials give:

| Stage | Median time | Cumulative Julia allocation |
|---|---:|---:|
| Initialization | `0.036709 s` | `72.157 MiB` |
| First forced-four sweep | `0.171570 s` | `193.927 MiB` |
| Steady second sweep | `0.132383 s` | `193.159 MiB` |

Fresh-process peak RSS was `855.0625 MiB`. Aligned projection took
`30.47646 s` versus `0.18147 s` for aligned cached execution, about `168x` on
this exact fixture. Numeric `blocksize=9` projection took `31.20375 s` but has
a different window geometry, so it is a control rather than a trajectory
comparison. All 182 package assertions and the fixed/ragged verifier passed.
These end-to-end numbers are paper-usable with this provenance, but they do
not establish a scaling law.

## Frozen Convergence And M5b Characterization

With `cutoff=1e-11`, `scf_cutoff=1e-12`, and the original cap of 50, the
solver did not converge. At sweep 50:

| Quantity | Result |
|---|---:|
| Last energy change | `9.64e-11 Ha` |
| Energy distance to oracle | `3.73e-10 Ha` |
| Independent returned-energy mismatch | `1.78e-15 Ha` |
| Full-Fock residual | `3.60e-5` |
| Minimum oracle singular value | `0.999999999632` |
| Maximum projector distance | `3.83e-5` |
| S2 / Z-spin deviations | `3.19e-5 / 1.47e-4` |

This immutable failure is part of the accepted record. The later ordinary
M5b run, from the original seed with an 80-sweep cap, stopped under the
solver's energy-change-only production policy at sweep 60. The diagnostic run
disabled the outer energy stop and completed 140 sweeps with the same local
`scf_cutoff=1e-12`. Its first 50 reduced rows reproduce the original M5
trajectory byte-for-byte.

| Gate | First satisfying sweep | Behavior through 140 |
|---|---:|---|
| Energy distance `<=1e-8 Ha` | 36 | remains satisfied |
| Minimum singular value `>=0.9999999999` | 56 | remains satisfied |
| Energy change `<1e-11 Ha` | 60 | remains satisfied |
| Projector distance `<=1e-5` | 62 | remains satisfied |
| S2 deviation `<=1e-6` | 76 | remains satisfied |
| Z-spin deviation `<=1e-6` | 83 | only sweeps 83--88 |

Residuals at sweeps 60, 80, 100, 120, and 140 are `1.1839560e-5`,
`1.3887388e-6`, `6.4627856e-7`, `6.3292824e-7`, and `6.3275583e-7`.
The minimum is `6.3260414e-7` at sweep 134, about `47.9x` the oracle residual.
At sweep 140, energy is within `1.62e-13 Ha` of the oracle, projector distance
is `6.88e-7`, minimum singular value is `0.999999999999881`, and the Z-spin
deviation is `2.81e-6`. The largest upward energy step after the plateau is
`2.56e-13 Ha`; roundoff-scale fluctuations begin after sweep 80, and this step
is below the `128*eps(Float64)` scale allowance of `4.14e-13 Ha`.

Each of the 140 rows records 98 sweep windows, 98 incremental absorptions, and
98 local cached windows, with zero cache fallbacks and zero projection windows.
The largest independent returned-energy mismatch is `2.20e-13 Ha`; the largest
orthonormality or idempotency error is `4.56e-15`. These facts establish branch
continuity and cache integrity, despite the transient historical Z-spin gate.
With `||Z||=22.87`, the Z-spin miss is also below the rough operator-norm bound
`2*22.87*6.88e-7 = 3.15e-5`; it is not evidence of branch loss.

The residual floor is consistent with, but does not prove, limitation by the
current local tolerance and damped mixed-density microiteration. No
stationarity-based production stopping criterion, tighter local policy,
determinant-variational update, polishing step, or acceleration method was
validated by this program.

## Accounting, Retained Evidence, And Next Boundary

Across maintained source, tests, scripts, examples, README, and architecture
documentation from the program base through M5b, the branch has 1,026 added
and 965 deleted lines, for net growth of 61. The runner is 274 lines. No M5 or
M5b production source or committed test change was made.

Durable reproduction artifacts are:

- [frozen-Be runner](../scripts/bench_radial_be.jl), parameterized by packet,
  output directory, modes, trials, maximum sweeps, and energy cutoff;
- [reduced 140-sweep trajectory](cached_sliced_performance_m5b_sweeps_2026-07-20.tsv),
  which contains the original 50 rows as an exact prefix; and
- the ordinary package tests and
  `scripts/verify_cached_vs_projection.jl` for maintained numerical coverage.

Raw logs, PID files, and serialized checkpoints were machine-local under
`~/dmrgtmp/hfdmrg_cached_sliced_20260719/`; they are provenance, not repository
artifacts.

After review, promote this cache branch to `main`. Begin any stationarity,
convergence, polishing, or acceleration study on a fresh branch. The old March
work predates the sliced correction and cache program: use it as reference and
port only independently reviewed determinant-variational or history pieces;
do not merge it wholesale.

The first convergence diagnostic should compare the sweep-140 state under a
tighter local SCF tolerance and with one full-space Fock/polishing step before
selecting any production algorithm change.

-- hfdmrg-manager@macmini
