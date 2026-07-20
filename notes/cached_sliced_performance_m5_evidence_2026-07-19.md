# Cached Sliced-Backend Milestone 5 Evidence

Date: 2026-07-19
Owner: `hfdmrg-manager`
Runner commit: `6054a2d44161a70bac8b96019fa22d33e333ecd6`
Status: **EXECUTED BUT NOT ACCEPTED**

## Plain-Language Outcome

The cached sliced backend is numerically faithful to both the projection
backend and an independent physical full-Fock calculation, and it clears the
frozen one-sweep time, allocation, memory, route-count, and alignment gates by
wide margins. A warmed forced-four sweep takes about `0.13-0.17 s`, allocates
about `193 MiB`, uses exactly 98 local windows and 490 Fock builds, and invokes
no projection fallback. Fresh-process peak resident memory is `855.0625 MiB`.

Milestone 5 is nevertheless not accepted. With the frozen solver settings and
50-sweep cap, the Be determinant is still approaching the saved oracle but has
not met the convergence, occupied-subspace, or branch-fingerprint gates. This
is a slow, smooth convergence tail rather than a cache disagreement. No source,
test, Hamiltonian, damping, local-iteration, convergence, or gate change was
made. Rough extrapolation suggests about ten additional sweeps, but those
sweeps are outside the frozen contract and were not run.

## Frozen Fixture And Execution Environment

The packet checksum was verified before execution:

```text
packet: /Users/srw/dmrgtmp/hfdmrg_radial_be_20260718/radial_be_lmax2_v1.ser
SHA256: fba29651c02cde3cdde90850d90ba30a1abd120ebfa8ac8b64cf2e5a92fea629
name: Be-radial-Ylm-s015-lmax2-v1
dimension: 468 = 52 radial slices * 9 angular functions
occupations: Nalpha = Nbeta = 2
```

The exact packet fields are
`metadata,H,V6,dims,seed_up,seed_dn,oracle_up,oracle_dn,Zop,radial_centers,radial_weights,angular_offsets,oracle_energy,oracle_residual,oracle_fock_angle,oracle_s2,oracle_zspin,oracle_parity,random_parity,current_backend_known_mismatch`.
The runner verified `H` as `468 x 468`, `V6` as
`9 x 9 x 9 x 9 x 52 x 52`, `dims == fill(9,52)`, angular offsets
`[1,2,5,10]`, the recorded `V6` permutation `(2,3,5,6,1,4)`, finite ordered
radial data, symmetric one-body operators, orbital shapes, and orthonormality.

All six runs recorded the same execution environment:

```text
host: rh310l.ps.uci.edu
CPU: Apple M4 Pro
Julia: 1.12.6
Julia threads: 1
BLAS threads: 8
BLAS: LBTConfig([ILP64] libopenblas64_.dylib)
project: repository Project.toml
depot: normal home depot; JULIA_DEPOT_PATH unset
HFDMRG_BENCH_TIMING: 1
allocation method: Base.gc_bytes differences
memory method: Sys.maxrss() in a fresh process
```

The timing mode alone used one warmup followed by three same-process raw
trials and reports their median and range.

## Independent Numerical Convention

For each spin, the runner forms the physical density

```text
D_sigma = C_sigma C_sigma'.
```

It contracts the frozen six-index interaction directly, without calling either
sliced backend:

```text
G_sigma[p,r] += sum(q,s) (D_alpha[q,s] + D_beta[q,s]) V[p,q,r,s]
G_sigma[p,s] -= sum(r,q) D_sigma[r,q] V[p,q,r,s]
F_sigma       = H + G_sigma.
```

Thus direct and exchange retain their signed first/third and first/fourth
output conventions. The independently recomputed determinant energy is

```text
E = 1/2 Tr[D_alpha(F_alpha + H)]
  + 1/2 Tr[D_beta (F_beta  + H)].
```

Occupied-virtual stationarity is the combined norm of
`F_sigma C_sigma - C_sigma(C_sigma' F_sigma C_sigma)`. Occupied-subspace
agreement uses minimum singular values and stable two-sided projection
distances. The branch fingerprints are

```text
S^2 = Sz(Sz+1) + Nbeta - ||C_alpha' C_beta||_F^2
Zspin = Tr[(D_alpha-D_beta) Zop].
```

Orthonormality and idempotency are evaluated from the small occupied Gram
matrices. Explicit window checks compare total local Fock matrices, including
the projected one-body term, for physical, projection, and cached routes.

## Fixture And Oracle Check

| State | Independent energy | Reported-energy mismatch | Residual | S2 | Zspin |
|---|---:|---:|---:|---:|---:|
| Distant seed | `-13.489740341026078` | n/a | `1.736001724185797` | `0.9999986209944836` | `-1.5135347128552081` |
| Frozen oracle | `-14.573351154439951` | `2.36e-13 Ha` | `1.32e-8` | `0.13165158762762275` | `-1.1194338511924447` |

The oracle alpha/beta minimum singular values are at least
`0.9999999999999997`; projector distances are at most `6.67e-16`, and
orthonormality/idempotency errors are `4.48e-16`. Recorded and recomputed
oracle Fock angles are both zero. The frozen oracle branch is recovered by the
independent calculation.

## Same-Range And Explicit Full-Fock Parity

One forced-four aligned sweep from the same distant seed gives:

| Quantity | Projection | Cached | Difference/gate |
|---|---:|---:|---:|
| Energy | `-14.57318656816803` | `-14.573186568167912` | `1.17e-13 Ha` / `1e-9` |
| Independent energy mismatch | `1.42e-14` | `1.51e-13` | `1e-10` |
| Alpha projector distance | | | `8.98e-11` / `1e-8` |
| Beta projector distance | | | `2.05e-10` / `1e-8` |
| Sweep windows / Fock calls | `98 / 490` | `98 / 490` | exact |

Cached initialization used two initial states and 49 incremental absorptions.
The sweep used 98 incremental absorptions, 98 local windows, zero fallback
absorptions, and zero projection windows.

Representative total-Fock maximum absolute disagreements were:

| Center slice | Cached vs projection | Projection vs physical full Fock |
|---:|---:|---:|
| 2 | `1.78e-15` | `2.27e-13` |
| 26 | `1.42e-14` | `2.84e-14` |
| 51 | `8.88e-16` | `3.82e-14` |

All are below the frozen `1e-12` Fock gate.

## Timing, Allocation, And Resident Memory

One warmup preceded three same-process trials. The table reports medians and
raw ranges; allocation is cumulative Julia allocation for the measured stage,
not resident memory.

| Stage | Median wall time | Three-trial range | Allocation | Counts |
|---|---:|---:|---:|---|
| Initialization | `0.036709 s` | `0.036558-0.039782 s` | `72.157 MiB` | 2 init, 49 incremental |
| First forced-four sweep | `0.171570 s` | `0.169527-0.173178 s` | `193.927 MiB` | 98 windows, 490 Focks |
| Steady second sweep | `0.132383 s` | `0.124844-0.173918 s` | `193.159 MiB` | 98 windows, 490 Focks |

Every measured sweep used 98 local windows, 98 incremental absorptions, and
zero cached or projection fallback. Both sweep stages are far below `40 s`
and `2 GiB`.

The fresh-process TSV captured `Sys.maxrss() = 896204800` bytes
(`854.6875 MiB`) before its final evidence write. The process's final logged
peak was `896598016` bytes, or `855.0625 MiB`; its sweep took `0.195355 s`
and allocated `203753120` bytes. The final logged value is the absolute
process peak RSS reported here, distinct from `Base.gc_bytes`.

## Frozen 50-Sweep Convergence Result

The aligned cached route used the fixed `cutoff=1e-11`,
`scf_cutoff=1e-12`, and 50-sweep cap. Every sweep used exactly 98 local
windows, 98 incremental absorptions, and zero fallback/projection windows.
Fock calls declined from 480 to 302 because local SCF steps were allowed to
converge early; the five-call count applies only to forced-four windows.

At sweep 50:

| Gate quantity | Result | Frozen gate | Decision |
|---|---:|---:|---|
| Solver converged | `false` | true by sweep 50 | **fail** |
| Last energy change | `9.64e-11 Ha` | `<1e-11 Ha` | **fail** |
| Independent energy mismatch | `1.78e-15 Ha` | `<=1e-10 Ha` | pass |
| Energy distance to oracle | `3.73e-10 Ha` | `<=1e-8 Ha` | pass |
| Orthonormality/idempotency | `3.64e-15` | `<=1e-11` | pass |
| Full-Fock residual | `3.60e-5` | stationarity required | not reached |
| Minimum alpha/beta singular value | `0.999999999632` | `>=0.9999999999` | **fail** |
| Maximum projector distance | `3.83e-5` | `<=1e-5` | **fail** |
| S2 deviation | `3.19e-5` | `<=1e-6` | **fail** |
| Zspin deviation | `1.47e-4` | `<=1e-6` | **fail** |

The complete reduced trajectory is in
`notes/cached_sliced_performance_m5_sweeps_2026-07-19.tsv`. Recent energy
changes contract smoothly by about `0.793` per sweep. A geometric extrapolation
from sweep 50 estimates roughly ten more sweeps to reach the energy cutoff.
That estimate is diagnostic only: the frozen 50-sweep cap was not widened and
no additional sweep was run.

## Alignment Control

| Route | Energy | Init | One sweep | Allocation | Windows/Focks |
|---|---:|---:|---:|---:|---:|
| Numeric `blocksize=9` projection | `-14.573004632463093` | `0.01696 s` | `31.20375 s` | `4570306976 B` | `100/500` |
| Aligned projection | `-14.57318656816803` | `0.00115 s` | `30.47646 s` | `4478711232 B` | `98/490` |
| Aligned cached | `-14.573186568167912` | `0.03890 s` | `0.18147 s` | `203752928 B` | `98/490` |

Numeric and aligned projection use different window geometries and therefore
different trajectories; this is a geometry/runtime control, not an algorithm
ranking. Identically aligned projection and cached routes satisfy the energy,
projector, returned-energy, and route-count gates.

## Package And Workflow Validation

The full package suite passed all `182` assertions:

```text
julia --project=. -e 'using Pkg; Pkg.test()'
```

The supported fixed/ragged verifier also passed:

```text
julia --project=. scripts/verify_cached_vs_projection.jl
```

Its largest fixed/ragged window disagreement was `5.11e-15`; its largest
sweep-energy disagreement was `3.55e-15 Ha`.

## Commands And Raw Records

The six exact recorded runner commands differed only by mode/output directory;
the timing command also specified three trials:

```text
/Users/srw/.julia/juliaup/julia-1.12.6+0.aarch64.apple.darwin14/Julia-1.12.app/Contents/Resources/julia/bin/julia --project=. /Users/srw/Library/CloudStorage/Dropbox/codexhome/work/hfdmrg/scripts/bench_radial_be.jl --packet=/Users/srw/dmrgtmp/hfdmrg_radial_be_20260718/radial_be_lmax2_v1.ser --outdir=/Users/srw/dmrgtmp/hfdmrg_cached_sliced_20260719/m5_be/fixture --modes=fixture
/Users/srw/.julia/juliaup/julia-1.12.6+0.aarch64.apple.darwin14/Julia-1.12.app/Contents/Resources/julia/bin/julia --project=. /Users/srw/Library/CloudStorage/Dropbox/codexhome/work/hfdmrg/scripts/bench_radial_be.jl --packet=/Users/srw/dmrgtmp/hfdmrg_radial_be_20260718/radial_be_lmax2_v1.ser --outdir=/Users/srw/dmrgtmp/hfdmrg_cached_sliced_20260719/m5_be/parity --modes=parity
/Users/srw/.julia/juliaup/julia-1.12.6+0.aarch64.apple.darwin14/Julia-1.12.app/Contents/Resources/julia/bin/julia --project=. /Users/srw/Library/CloudStorage/Dropbox/codexhome/work/hfdmrg/scripts/bench_radial_be.jl --packet=/Users/srw/dmrgtmp/hfdmrg_radial_be_20260718/radial_be_lmax2_v1.ser --outdir=/Users/srw/dmrgtmp/hfdmrg_cached_sliced_20260719/m5_be/timing --modes=timing --trials=3
/Users/srw/.julia/juliaup/julia-1.12.6+0.aarch64.apple.darwin14/Julia-1.12.app/Contents/Resources/julia/bin/julia --project=. /Users/srw/Library/CloudStorage/Dropbox/codexhome/work/hfdmrg/scripts/bench_radial_be.jl --packet=/Users/srw/dmrgtmp/hfdmrg_radial_be_20260718/radial_be_lmax2_v1.ser --outdir=/Users/srw/dmrgtmp/hfdmrg_cached_sliced_20260719/m5_be/convergence --modes=converge
/Users/srw/.julia/juliaup/julia-1.12.6+0.aarch64.apple.darwin14/Julia-1.12.app/Contents/Resources/julia/bin/julia --project=. /Users/srw/Library/CloudStorage/Dropbox/codexhome/work/hfdmrg/scripts/bench_radial_be.jl --packet=/Users/srw/dmrgtmp/hfdmrg_radial_be_20260718/radial_be_lmax2_v1.ser --outdir=/Users/srw/dmrgtmp/hfdmrg_cached_sliced_20260719/m5_be/alignment --modes=alignment
/Users/srw/.julia/juliaup/julia-1.12.6+0.aarch64.apple.darwin14/Julia-1.12.app/Contents/Resources/julia/bin/julia --project=. /Users/srw/Library/CloudStorage/Dropbox/codexhome/work/hfdmrg/scripts/bench_radial_be.jl --packet=/Users/srw/dmrgtmp/hfdmrg_radial_be_20260718/radial_be_lmax2_v1.ser --outdir=/Users/srw/dmrgtmp/hfdmrg_cached_sliced_20260719/m5_be/rss --modes=rss
```

Raw logs are under
`~/dmrgtmp/hfdmrg_cached_sliced_20260719/m5_be/{fixture,parity,timing,convergence,alignment,rss}.log`.
The corresponding PIDs were `34415`, `34445`, `34477`, `34520`, `34558`,
and `34595`; PID files are in each mode directory as
`bench_radial_be.pid`. Long runs were monitored through their active tool
sessions and mode logs; the recorded process-polling command is

```text
ps -p "$(cat ~/dmrgtmp/hfdmrg_cached_sliced_20260719/m5_be/MODE/bench_radial_be.pid)" -o pid,etime,rss,vsz,pcpu,pmem,command
```

The completed package and verifier logs are
`~/dmrgtmp/hfdmrg_cached_sliced_20260719/m5_be/package_test.log` and
`~/dmrgtmp/hfdmrg_cached_sliced_20260719/m5_be/verify_cached_vs_projection.log`.

## Scope And Line Accounting

The retained runner commit adds only `scripts/bench_radial_be.jl`: `275`
additions and no deletions. It exceeds the preferred `220`-line runner cap by
`55` lines, but adds no public API, framework, production source, or committed
test. M5 execution itself made no production or test change. This evidence
memo, the reduced sweep TSV, and replacement status text are the only pending
documentation/evidence changes.

The next action is paper-manager review of this failed frozen acceptance gate.
No extra sweep, source repair, policy change, or softened gate is authorized.

-- hfdmrg-manager@macmini
