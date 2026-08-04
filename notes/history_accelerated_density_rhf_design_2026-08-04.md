# History-Accelerated Density-Density RHF Design

Date: 2026-08-04

Status: approved for bounded private implementation after normalized
production-main acceleration and external stationary-basin validation.

## Decision

Implement the first package version as a private density-density RHF driver
that composes complete ordinary HFDMRG solves. Do not change the sweep core,
backend API, density-density backend, or any public `solve_hfdmrg` route.

The driver alternates two complementary operations:

1. ordinary HFDMRG sweeps discover structured orbital changes in small moving
   windows; and
2. a small HF solve in the anchored union of recent sweep-end occupied spaces
   applies those changes globally.

Every proposed history determinant is checked in the full physical basis. A
failed, unsafe, or overly large proposal falls back to ordinary sweeps from the
unchanged incoming determinant. The first implementation deliberately rebuilds
HFDMRG environments for every cleanup segment; block reuse is a separate
optimization.

## Authority And Evidence

Production authority for this design is `main@4164cf86b9734bb1ae8d55d2e0885df6895335f3`
with tree `5b89afca40b3a2139cc94f9c608f5aab4833bd49`.

The accepted mechanism evidence used the numerically later experimental source
tree `af073b64d5579424b7b8f3ed60d1d0434ae42b6d`. The history feature did not
depend on its frozen-occupied work, but implementation validation must repeat
the physical gates from production `main`; tree similarity is not assumed.

Authenticated inputs and evidence are:

- WL packet SHA-256
  `8e4a4d8bd32eeeef9be4118a3b3c738466f4f823cec12175f499fc127f6c38d6`;
- PQS packet SHA-256
  `1f3ccacd50d7d42cb7029cf91ee856cb47dffdacae5373c1fd1a2e2a45a5bf17`;
- original review-bundle manifest
  `4d68939a0e5817125ce2d47f1d7f263a464485471b5ca43e36e0ea998924a2e8`;
- hardening-addendum manifest
  `4ab0d00477eb0d16f3a8d20aab4857df5997d296fd736b890652c50c34e8dc4d`.

The hardening proved commutator DIIS, quadratic energy evaluation, clean
fallback, the `1e-6` WL/PQS endpoints, and common WL stationary-basin parity.
It did not implement the requested history-count normalization. Consequently,
the raw response-Gram cutoff `1e-10` remains mechanism evidence, not the package
rank rule. The normalized rule below is a new implementation gate and requires
a fresh WL/PQS validation before acceptance.

### Normalized Production-Main Replay

The fixed-policy replay used the design authority revision with SHA-256
`46055f6ccf32cf29b3777b82ac8c6ca5f6fa8f427d30cfb18cea0090df27b710`
and scratch runner SHA-256
`06f3b1953ee2d2f47941b8df3f9ec60b572480ddcdc498b35cee947c1e3c7b7f`.
It ran from production `main`, with one Julia thread, ten BLAS threads, and no
change to either normalized history cutoff. Two earlier preflights stopped
before loading a physical packet: the first deliberately exposed that its
compact model exceeded the workspace guard, and the second exposed the need
for the unit-scale numerical-rank anchor now specified above.

The successful replay gave:

| Route | Ordinary 100-sweep residual | Accelerated sweeps/cycles | Time to `1e-4` | Final residual | Driver phase time | Resource fallback |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| WL | `7.618813904209707e-4` | 36 / 15 | `16.988396958 s` | `5.758068423206059e-7` | `25.795514666 s` | 0 |
| PQS | `1.556268013985071e-4` | 36 / 15 | `20.064558791 s` | `8.212746622142931e-7` | `25.056859375 s` | 1 |

WL used contracted ranks 18--32. PQS used ranks 18--32 and cleanly rejected
one selected rank-33 model at the predeclared pair-map/model workspace guard;
ordinary cleanup retained the incoming determinant and the route subsequently
reached the target. There were no Pulay rejections. Independent readback
recomputed all four endpoint energies, residuals, and Gram errors and verified
every output manifest entry.

At the same 36-sweep count, ordinary WL and PQS residuals were
`7.993594508671634e-4` and `4.383917209569072e-4`, respectively. Accelerated
outer wall times were `26.245845834 s` and `25.422487959 s`; the 100-sweep
ordinary controls took `75.075037334 s` and `64.166645875 s`. Accelerated
cumulative allocations were `55.285006832 GB` and `54.194120704 GB`, while
peak RSS was `2.608594944 GB` and `2.567487488 GB`. The corresponding
100-sweep control values were `112.797004112 GB`/`2.255077376 GB` and
`108.693339152 GB`/`2.201747456 GB`. These are engineering diagnostics from
one process per route, not publication timing or matched-sweep memory gates.

Against the previously hardened references, the accelerated raw-target
energies differ by `-9.706013770482969e-12 Ha` (WL) and
`1.8260948309034575e-11 Ha` (PQS), and both occupied--virtual gaps are positive
(`0.3391996496223715 Ha` and `0.33465141713347657 Ha`). The raw determinants
at the fixed `1e-6` stopping target differ from those references by projector
distances `4.4694235036968176e-5` and `1.4228783218720863e-4`. Those distances
are diagnostics, not stationary-projector gates: a positive Fock gap does not
bound the self-consistent HF response Hessian, and a soft response direction
can produce this motion at residual `1e-6` without changing basins.

Validation-only floor-aware full-space commutator DIIS was therefore applied
independently to each sealed accelerated endpoint and its matched production-
main ordinary-control endpoint. This solver is not part of the package driver
or acceleration timing. All four endpoints reached `norm(FD-DF) <= 2e-11`:

| Route | Control `K` | Accelerated `K` | Polished energy difference | Polished projector distance | Gate |
| --- | ---: | ---: | ---: | ---: | ---: |
| WL | `1.496909266798822e-11` | `1.2227172087338298e-11` | `-4.362732397567015e-12 Ha` | `2.139100693159639e-9` | pass |
| PQS | `1.8059472742285123e-11` | `1.6970822030310353e-11` | `1.2150280781497713e-11 Ha` | `4.7513499860296286e-8` | pass |

The polished comparisons pass the existing `2e-9 Ha` energy and `2e-6`
projector gates by wide margins. Polishing took `101.570312709 s` and
`23.010359042 s` from the WL control and accelerated endpoints, and
`32.464640666 s` and `27.343971209 s` for PQS. These validation costs remain
separate from the acceleration result. Independent quadratic-energy readback
gave pair differences `-9.023892744153272e-12 Ha` and
`1.2292389328649733e-11 Ha`; the few-picohartree change from the runner values
is contraction-order roundoff and does not approach the gate. Neither history
cutoff was tuned.

The sealed replay root is
`~/dmrgtmp/hfdmrg_history_rhf_normalized_main_20260804`. Its four manifest
SHA-256 values are, respectively, WL control
`0c7566b9a673fd401600e140ad2dbde050ca8cef9649f41bd4a4f54e1fff4c58`,
WL accelerated
`2e0b0e798dfb5afa8c9dd56bbb4c58e33b9da092160f288a2a356916ce9d4658`,
PQS control
`5923255159bcfd5b72e7d6a4f9c181cecba5780522e0bcecb792afdf14d9d9bd`,
and PQS accelerated
`168f35aec1c90e03b9ce2219d3d9258823f92638d0c2423675f8c831c1323f32`.
The external closure used runner SHA-256
`92986368d72fea6b3d033723b315200ada7ab0ec80860354429c161e3c500426`
and produced manifest SHA-256
`ab49830d58e8a155d1273f37e2f8e748a59d784228242f6fd7476f88f288770a`.

## Scope

The bounded implementation includes only:

- real density-density RHF;
- a private fixed development policy;
- complete sweep-end history snapshots;
- an anchored, normalized response basis;
- contracted gauge-invariant RHF with rank-revealing Pulay mixing;
- full-basis energy, residual, orthogonality, and overlap acceptance checks;
- fresh two-sweep HFDMRG cleanup segments; and
- a private result and trace sufficient for numerical and cost review.

It excludes UHF, split-H, sliced interactions, target residuals, frozen
occupied environments, residual-vector enrichment, retained block state,
adaptive cleanup length, a public observer/trace extension, a public keyword,
and a generic subspace-acceleration framework.

The first implementation accepts only `Matrix{Float64}` inputs for `H`, `V`,
and `C0`. The fixed physical, covariance, and DIIS gates are not meaningful for
every `AbstractFloat`; in particular, their nominal scale would be far too
loose for `Float32`. The private driver validates dimensions, finiteness,
physical symmetry, occupation, and initial orthogonality with fixed
scale-aware `Float64` gates. It does not convert or copy `V`. Mixed floating
types, other floating types, abstractly typed storage, and integer inputs
reject before a sweep or history allocation.

## Physical Convention

Let `C` be an `N x nocc` real orthonormal occupied-orbital matrix for one spin
and let

```text
D = C C'
d = diag(D)
F(D) = H + 2 Diagonal(V d) - V .* D.
```

The electronic RHF energy is

```text
E(D) = 2 dot(D,H) + 2 dot(d,V*d) - sum(V .* D .* D).
```

The orbital residual and gauge-invariant projector residual are

```text
R = F C - C (C' F C)
K = F D - D F.
```

For an idempotent real density, `norm(K) = sqrt(2)*norm(R)` up to roundoff.
The driver uses `norm(R)` as its target. It never forms the full physical
commutator `F*D-D*F`, which would require cubic matrix multiplication. Only the
small contracted commutator `Kq` is used as the Pulay error. The full audit
forms `F*C` and uses the orbital residual `R`. It also never forms `D*(F+H)`
merely to compute energy.

All accepted and returned states are pure determinants. Damped local densities
remain internal to ordinary HFDMRG segments and never enter the history queue.

## History State And Lifecycle

The fixed first policy is:

```text
bootstrap sweeps       6
bootstrap captures     2, 4, 6
history FIFO length    6, including the current determinant
cleanup sweeps         2 after every history attempt
maximum cycles         25
full residual target   1e-6
contracted DIIS cap    40 iterations
contracted DIIS memory 7 Fock/error pairs
minimum overlap        0.99
```

The bootstrap is one six-sweep `solve_hfdmrg` call with `cutoff=0`; an internal
observer copies the complete physical determinant after sweeps 2, 4, and 6.
This pays environment initialization once. A history proposal is then followed
by one fresh two-sweep solve, also with `cutoff=0`. The second sweep amortizes
the fresh environment construction and was materially better than one-sweep
cleanup in the scratch ladder.

Only the final determinant of each two-sweep cleanup is appended to the FIFO.
Rejected proposals append the endpoint obtained by cleaning up the unchanged
preproposal determinant. The driver tests its full residual only after a
complete bootstrap or cleanup segment and never returns a raw history proposal.

The driver stops when the full residual target is met or the cycle cap is
exhausted. `target_reached` means only that the independent full residual is at
or below the private target; it is not a claim of global-minimum stability.

## Anchored Normalized History Basis

Let the newest FIFO determinant be `C` and let `C_j`, `j=1,...,h`, be the older
snapshots. Define the normalized response matrix

```text
Z = (I - C C') [C_1 C_2 ... C_h] / sqrt(h).
```

Apply the projector twice to suppress roundoff leakage, then compute a thin SVD
`Z = U*S*W'`. This construction is invariant under independent occupied
rotations or sign changes in every snapshot. Its covariance is

```text
Z Z' = (1/h) sum_j (I-D) D_j (I-D),
```

so `lambda_i = S[i]^2` is comparable across FIFO lengths.

The current occupied space is never rank reduced. The contracted basis is

```text
Q = orth([C, U[:,1:k]]).
```

The response rank uses three bounds:

1. Numerical rank retains only singular values above
   `64*eps(Float64)*max(N,h*nocc)*max(1,S[1])`. The unit-scale anchor is
   essential: when all snapshots have the current projector, an entirely
   relative `S[1]` threshold would promote arbitrary roundoff singular vectors
   instead of classifying a zero response.
2. The normalized absolute covariance floor is `lambda > 2e-11`. This is the
   FIFO-six analogue of the successful raw `1e-10` floor because a full FIFO
   has five older snapshots.
3. The relative rule chooses the smallest `k_tail` satisfying

   ```text
   sqrt(sum(lambda[k_tail+1:end]) / sum(lambda)) <= 1e-4.
   ```

Use `k = min(k_numerical, max(k_absolute, k_tail))`. Thus any resolvable mode
above the absolute floor is retained, while the tail rule may conservatively
retain smaller modes when their combined response remains meaningful. If the
response is numerically zero, skip the contracted solve as a recorded no-op.

The fixed `2e-11` and `1e-4` values are private development choices, not energy
error estimates or public controls. Promotion requires the physical gates and
a rank/cost study; failure does not authorize tuning them on the same attempt.

Before contraction require

```text
norm(C - Q*(Q'*C)) <= 32*eps(Float64)*sqrt(N*size(Q,2))
norm(Q'*Q - I)     <= 32*eps(Float64)*sqrt(N*size(Q,2)).
```

Failure is fail-closed and leaves the incoming determinant unchanged.

## Contracted Density-Density Model

For `r = size(Q,2)`, construct

```text
Hq[a,b] = sum_pq Q[p,a] H[p,q] Q[q,b]
P[p,(a,b)] = Q[p,a] Q[p,b]
G[(a,b),(c,d)] = P[:,(a,b)]' * V * P[:,(c,d)].
```

`G` is the complete four-index interaction only inside the small history
space. No physical `N^4` object is created. Store one `r x r x r x r` tensor;
do not allocate a second exchange-permuted copy.

For contracted density `Dq = c*c'`, use

```text
Fq[a,b] = Hq[a,b]
        + 2 sum_cd G[a,b,c,d] Dq[c,d]
        -   sum_cd G[a,c,d,b] Dq[c,d]
Eq = dot(Dq, Fq + Hq).
```

Symmetrize `Fq` only at the eigensolve boundary. Compact tests must prove that
`Fq`, `Eq`, and the lifted determinant agree with direct physical contractions.
`Eq` is the exact physical energy restricted to `span(Q)`, not a surrogate.

## Contracted RHF And Pulay Solve

Start from `c0 = Q'*C`. At each iteration:

1. form `Dq`, `Fq`, energy, orbital residual, and `Kq=Fq*Dq-Dq*Fq`;
2. append `Fq` and `Kq` to a FIFO of length seven;
3. form the Pulay Gram matrix `B[i,j]=dot(K_i,K_j)`;
4. divide its Gram block by `maximum(abs,B)`;
5. solve the augmented constrained system with SVD, retaining singular values
   above `128*eps(T)*(m+1)*smax`; and
6. diagonalize the mixed symmetric Fock and occupy its lowest `nocc` vectors.

No fixed diagonal ridge or unqualified inverse is permitted. Let `s=sum(c)`.
After the rank-revealing solve, require finite coefficients, coefficient norm
at most `1e6`, and

```text
abs(s-1) <= 1024*eps(Float64)*max(1,norm(c,1)).
```

If that gate passes, replace `c` by `c/s` before mixing Focks, enforcing the
affine constraint to roundoff. A truncated constraint or any failed coefficient
gate records a Pulay failure and uses the current unmixed Fock for that
iteration. Two consecutive Pulay failures reject the complete proposal.
Programming errors are not caught and disguised as fallback.

The contracted convergence test is

```text
norm(Kq) <= 1e-10 * max(1,norm(Fq)).
```

Seed the best valid determinant with the incoming `c0` before the first update,
so the no-change determinant is always retained. Track later candidates with a
scale-aware energy tie

```text
tau_E = 128*eps(Float64)*max(1,abs(E_best),abs(E_trial)).
```

A lower energy by more than `tau_E` wins. Within the tie, the smaller
commutator residual wins. This prevents a floating-point-low but less
stationary iterate from replacing a near-equal stationary iterate. At a hard
iteration cap, return this tie-aware best contracted determinant with explicit
status.

## Full-Basis Acceptance And Fallback

Lift `Ctrial = Q*c` and independently recompute its full physical energy,
residual, and orthogonality. Accept only if:

- all values are finite and the occupation count is unchanged;
- the orbital Gram error passes the same scale-aware matrix gate;
- `Etrial <= Ebefore + tau_E`;
- `minimum(svdvals(Cbefore'*Ctrial)) >= 0.99`; and
- the contracted/direct Fock and energy orientation is certified by tests.

Immediate residual growth is diagnostic, not a rejection condition. The WL
mechanism included a useful downhill proposal whose residual rose by `5.72x`
before cleanup.

Any ordinary numerical rejection runs the two cleanup sweeps from `Cbefore`.
No partially accepted orbital, history mutation, or contracted state survives.
Only after cleanup succeeds is its endpoint appended atomically to the FIFO.

## Pair-Map/Model Workspace Guard

Model construction is allowed only when the dominant pair-map/model workspace
satisfies

```text
2*N*r^2 + r^4 <= N^2
```

in `Float64` elements. The left side covers `P`, `V*P`, and `G`; it does not
claim to bound the history queue, SVD, DIIS arrays, full physical audit, Julia
allocator high-water state, or total process RSS. Evaluate every multiplication
and addition with overflow-safe checked integer arithmetic before allocating
`P`. Construct the model in its own scope so `P` and `V*P` are unreachable
before the full physical audit begins. Exceeding the guard is an ordinary
recorded fallback, not a rank truncation that violates the selected tail rule.

For WL (`N=2551`) the exact inequality admits `r=32` and rejects `r=33`.

This guard is intentionally conservative and private. It prevents a long
history or unexpectedly high response rank from creating hidden large state.

## Return Semantics

The private driver returns a private result containing:

- the RHF orbital matrix (shared by alpha and beta);
- independently recomputed determinant energy and residual;
- `target_reached` and terminal reason;
- sweep and history-cycle counts;
- accepted, rejected, no-op, and resource-fallback counts;
- contracted ranks and DIIS statuses; and
- phase times and allocated bytes needed for review.

The reported energy is always the physical energy of the returned determinant,
not the last damped local-SCF scalar. No public `solve_hfdmrg` return convention
changes.

## Code Organization And Closed File Set

The proposed implementation file set is:

- `src/HFDMRG.jl`: one private include;
- `src/history_accelerated_rhf.jl`: all private history, contraction, DIIS,
  orchestration, and result code; and
- `test/runtests.jl`: compact numerical and workflow tests.

Do not edit `src/core.jl`, `src/backend_api.jl`, density backend files, README,
or manifests. Do not export a symbol or add a public keyword. The provisional
entrypoint is qualified and private, conceptually
`HFDMRG._solve_history_rhf(H,V,C0; _policy, solver_controls...)`.

Only ordinary density-RHF controls needed by the physical gates are forwarded:
block size/partition, center width, local SCF cutoff, environment cutoff, and
verbosity. Outer cutoff and sweep counts belong to the fixed history policy;
an external observer is not supported in this first implementation.

## Cost Model

Let `h` be the number of older snapshots, `n` the occupied count, `r=n+k` the
history rank, and `N` the physical dimension.

| Phase | Work | Additional temporary/retained memory |
| --- | --- | --- |
| History queue | none between captures | `O((h+1)Nn)` retained |
| Anchored response/SVD | `O(N(hn)^2)` | `O(Nhn)` |
| One-body projection | `O(N^2 r)` | `O(Nr + r^2)` |
| Pair map and interaction | `O(N^2r^2 + Nr^4)` | `O(Nr^2 + r^4)` retained, approximately twice `Nr^2` at construction |
| Contracted DIIS iteration | `O(r^4+r^3)` | `O(r^4 + m r^2)` |
| Full acceptance audit | `O(N^2n)` | `O(N^2)` temporary |
| Conventional dense HF iteration | `O(N^3)` | `O(N^2)` |

The history correction is useful when `r << N` and the structured sweeps are
cheap enough that a coherent global mode otherwise dominates sweep count. It
is not promised to beat conventional dense HF on every `N≈2500` fixture.

Observed WL ranks were 6--28 in `N=2551`, and 18 hardened corrections used
`0.349 s` for model setup plus `0.024 s` for contracted HF. The remaining
workflow time was dominated by fresh cleanup sweeps and full diagnostics.

## Validation Ladder

Implementation proceeds only in this order.

### 1. Compact algebra and invariance

- direct physical versus contracted Fock and energy parity `<=2e-12`;
- analytic quadratic energy versus the accepted trace convention `<=2e-12`;
- invariance of history projector, rank, energy, and proposal projector under
  deterministic occupied rotations and sign changes `<=2e-12`;
- normalized rank behavior for FIFO lengths 3, 4, and 6 with deterministic
  covariance spectra straddling both cutoffs;
- exact current-space containment and orthogonality under their scale-aware
  gates; and
- rank-revealing Pulay behavior on full-rank and deficient Gram systems.

Tests assert numerical behavior and returned determinants, not private helper
names, buffer shapes, or internal vocabulary.

### 2. Compact complete workflow

- six-sweep capture followed by at least two history/cleanup cycles;
- one accepted proposal with full energy/projector parity;
- one rejected proposal whose endpoint equals ordinary cleanup from the
  unchanged preproposal determinant;
- one rotated-history run with identical physical result;
- one zero-response no-op; and
- one resource-guard fallback with no large allocation.

The ordinary public density RHF result must remain bitwise identical when the
private driver is not called.

### 3. Package and backend regression

- complete package suite;
- existing cached/projection verifier; and
- `git diff --check`.

The public symbol list and every public default route remain unchanged. Exact
default-parity controls cover density RHF/UHF, split-H UHF, target residual,
projection sliced, and cached sliced calls because merely including the private
file must be arithmetically inert.

### 4. Authenticated WL gate on production main

Use the same projected seed, physical ordering, six-sweep bootstrap, FIFO six,
and two-sweep cleanup policy. Under the new normalized rank rule require:

- residual `<=1e-4` in less than `25 s` on the matched host/thread setup;
- residual `<=1e-6` by the cycle cap;
- every accepted/rejected/fallback event to satisfy its contract;
- raw-target energy and projector distances to the conventional reference,
  reported as diagnostics rather than stationary gates;
- validation-only polishing of the sealed accelerated and matched ordinary
  endpoints to full-space commutator residual `<=2e-11`, followed by energy
  and projector parity within `2e-9 Ha` and `2e-6`; and
- reported ranks, tail ratios, setup/DIIS/cleanup/audit time, allocation, and
  peak RSS.

A clean resource fallback proves safety but does not pass the acceleration
gate. If the normalized WL route misses `1e-6` for any reason, including the
`r=32` boundary, stop for review before package implementation.

Run an ordinary-HFDMRG control from the identical seed on the same production
`main` tree, host, Julia/BLAS thread setup, and physical ordering. Record its
residual and energy trajectory through at least the accelerated route's total
sweep count, plus the established 100-sweep comparison when practical. The
normalized route is compared only with this matched production-main baseline,
not with the later experimental-tree trajectory.

This is engineering timing, not publication timing.
The external full-space closure is a basin test and its time is never included
in the package driver or acceleration timing.

### 5. Matched PQS gate

With the identical fixed policy, PQS must reach `1e-4` and either reach `1e-6`
or continue through clean ordinary fallback with no cost or integrity
pathology. Preserve its stationary basin under the same external `2e-11`
commutator closure and polished energy/projector gates used for WL. No
parameter tuning follows a failed attempt.

## Line Budget And Shrinkage

Preferred additions are:

```text
src/history_accelerated_rhf.jl   <= 230 lines
src/HFDMRG.jl                    1 line
test/runtests.jl                 <= 120 lines
total                            <= 351 lines
hard maximum                     380 lines
```

If the hard limit cannot hold the numerical contract, stop for review rather
than widening the feature. No existing production code is expected to be
deleted in the first implementation. Do not copy scratch serialization,
portable-runner authority machinery, conventional full-space controls, or
paper diagnostics into the package.

## Promotion Boundary

Passing this ladder authorizes only an experimental private density-RHF
implementation checkpoint. Public exposure requires a second design after:

- normalized history ranks are studied across a controlled size family;
- total time is compared with both ordinary HFDMRG and conventional HF at a
  common stationarity target;
- memory remains bounded by the declared guard; and
- at least one additional scientific family confirms the mechanism.

UHF, sliced interactions, cache reuse, adaptive triggering, and residual-vector
enrichment remain separate future milestones.
