# History-Accelerated Density-Density UHF Design

Date: 2026-08-04

Status: proposed for bounded private implementation; source implementation and
sliced-backend work are not yet authorized by this design checkpoint.

## Decision

Port the accepted private density-density RHF history accelerator to
density-density UHF without changing the sweep engine, backend API, public
solver interface, or history policy.  The first port supports both the common
one-body route and genuinely split `Hup`/`Hdn` routes.

Each history entry is one paired determinant `(Cup, Cdn)`.  A single contracted
one-particle basis spans the current alpha and beta occupied spaces and the
rank-revealed response of both spin histories.  Inside that basis, alpha and
beta retain separate densities, Fock matrices, occupations, and eigensolves.
One Pulay coefficient vector is obtained from both spin commutators and mixes
the paired Focks.  A proposal is accepted, rejected, or sent to resource
fallback atomically across both spins.

This is a composition layer around complete ordinary HFDMRG sweep segments.
It neither retains blocks between segments nor changes how an ordinary UHFDMRG
sweep builds, transforms, or discards environments and interaction caches.

## Authority And Existing Evidence

The design authority is production
`main@031da8f566f0992dd07f8232a7b7bbe51ea9d446`, tree
`4a9b24474cee76c14a44f5777a2711ae8d87a743`.

That tree contains:

- the accepted private density-density RHF history driver;
- guarded lowest-occupied subset eigensolves in ordinary RHF, common-H UHF,
  and split-H UHF; and
- no public history keyword or exported history symbol.

The RHF driver has passed compact numerical tests, the complete package suite,
WL/PQS acceleration, and stationary-basin closure.  The UHF port reuses its
fixed sweep/history policy and numerical constants unchanged.  It does not
infer that RHF evidence proves UHF basin preservation.

A scratch-only algebra probe at
`/tmp/history_uhf_algebra_probe.jl`, SHA-256
`33551166272ab932a3646e5982b90dd507869aee1c2e060e48d1df36667a1049`,
gave:

```text
alpha contracted/direct Fock error       5.576870229801381e-15
beta contracted/direct Fock error        4.5444669795149594e-15
contracted/direct total-energy error     1.3322676295501878e-15
rotated-history basis-projector error    2.4218531962557843e-14
```

This establishes the tensor orientation and paired-basis algebra on a random
`N=19`, rank-eight, spin-asymmetric fixture.  It is not implementation or
scientific acceleration evidence.

## Scope

The bounded implementation includes only:

- real `Matrix{Float64}` density-density UHF;
- common `H` and split `Hup`/`Hdn` one-body routes;
- paired complete-sweep history snapshots;
- a common anchored and normalized spin-history basis;
- contracted UHF with separate spin Focks and coupled commutator DIIS;
- full-basis energy, residual, Gram, overlap, and spin diagnostics for each
  proposal, with Aufbau diagnostics only at termination or external closure;
- atomic proposal acceptance and ordinary-sweep fallback; and
- a private result/trace sufficient for correctness and resource review.

It excludes sliced and cached-sliced interactions, target residuals, retained
block state, block/cache reuse, frozen occupied modes, residual-vector
enrichment, adaptive scheduling, public controls, and a generic accelerator
framework.  It does not modify or tune the accepted RHF policy.

The first port accepts one empty spin channel, which is a valid fully polarized
UHF determinant.  Empty-spin residual, Gram, and overlap contributions are
defined as zero, zero, and one.  At least one electron is required, and neither
spin occupation may fill the complete physical basis.  This behavior must be
tested explicitly rather than arising accidentally from empty matrix slices.

## Physical UHF Convention

Let `Ca` and `Cb` be real, separately orthonormal occupied-orbital matrices,
with fixed ranks `na` and `nb`.  Define

```text
Da = Ca*Ca'
Db = Cb*Cb'
dt = diag(Da) + diag(Db)

Fa = Ha + Diagonal(V*dt) - V .* Da
Fb = Hb + Diagonal(V*dt) - V .* Db.
```

For common-H UHF, `Ha === Hb === H`; the two Focks are nevertheless distinct
when the spin densities differ.  For split-H UHF, only the one-body matrices
change.  The interaction algebra is identical.

The determinant energy is

```text
E = 0.5*dot(Da, Fa + Ha) + 0.5*dot(Db, Fb + Hb).
```

This is the existing HFDMRG UHF convention.  The driver must not substitute an
RHF factor of two or form a spin-averaged Fock.

For each nonempty spin,

```text
Rs = Fs*Cs - Cs*(Cs'*Fs*Cs)
Ks = Fs*Ds - Ds*Fs.
```

The full target is

```text
max(norm(Ra), norm(Rb)) <= 1e-6,
```

with the empty-spin norm defined as zero.  Both values and their root-mean-
square are reported.  The max criterion prevents one well-converged spin from
hiding a poorly converged partner and reduces exactly to the RHF target when
the two spin projectors coincide.

The spin diagnostics are

```text
Sz = (na - nb)/2
S2 = Sz^2 + (na + nb)/2 - norm(Ca'*Cb)^2.
```

`Sz` is fixed by occupation.  `S2` may change during legitimate UHF
optimization and is diagnostic, not an acceptance gate.  It becomes a basin
comparison after independent stationary polishing.

The full physical audit uses orbital residuals and quadratic energy
contractions.  It never forms either full commutator `Fs*Ds-Ds*Fs`, which
would introduce unnecessary cubic work.

## Paired History And Common Basis

The queue stores immutable-by-convention pairs

```text
[(Ca_1,Cb_1), ..., (Ca_h,Cb_h), (Ca,Cb)].
```

Only complete-sweep observer snapshots and complete cleanup endpoints enter
the queue.  Damped local densities never do.

First construct a current spin-union anchor.  For

```text
X0 = [Ca Cb],
```

use a thin SVD and retain only its scale-aware numerical rank.  Call the
resulting orthonormal basis `A`.  Require separately

```text
norm(Ca - A*(A'*Ca)) <= 32eps(Float64)*sqrt(N*size(A,2))
norm(Cb - A*(A'*Cb)) <= 32eps(Float64)*sqrt(N*size(A,2)).
```

An empty spin passes by definition.  Failure is fail-closed.  In an RHF-like
state with identical alpha and beta occupied spaces, the anchor rank is `n`,
not `2n`.

Let `q` be the number of nonempty spin channels and let `h` be the number of
older paired snapshots.  Form

```text
Z = (I - A*A') [Ca_1 Cb_1 ... Ca_h Cb_h] / sqrt(q*h),
```

omitting empty channels.  Apply the anchor projector twice before the thin
SVD.  Its covariance is

```text
Z*Z' = (1/(q*h)) sum_j,s (I-A*A') D_s,j (I-A*A').
```

The normalization removes the trivial overall scaling with the number of older
snapshots and populated spin channels.  It does not erase the genuine spectral
differences carried by distinct alpha and beta histories.  When
`Ca_j == Cb_j` for all `j`, the duplicate spin terms cancel the factor of two
and reproduce the accepted RHF covariance exactly.

Use the unchanged RHF rank rules:

1. numerical singular-value threshold
   `64eps(Float64)*max(N,q*h*nmax)*max(1,S[1])`;
2. normalized covariance floor `S[i]^2 > 2e-11`; and
3. smallest relative rank whose discarded covariance tail is at most `1e-4`.

As in RHF, select

```text
k = min(k_numerical, max(k_absolute, k_tail))
Q = orth([A U[:,1:k]]).
```

Require `Q` to be orthonormal and to contain `Ca` and `Cb` separately under
the scale-aware gates above.  Independent occupied rotations and sign changes
within every spin and snapshot must leave the `Q` projector, selected rank,
contracted energy, and proposed spin projectors unchanged to roundoff.

## Contracted UHF Model

For `r=size(Q,2)`, construct one shared pair map and interaction tensor:

```text
Hqa = Q'*Ha*Q
Hqb = Q'*Hb*Q
P[p,(a,b)] = Q[p,a]*Q[p,b]
G[(a,b),(c,d)] = P[:,(a,b)]' * V * P[:,(c,d)].
```

When `Ha === Hb`, compute one projected one-body matrix and retain the same
object for both spin fields.  Split-H computes two matrices.  Do not duplicate
`G`, construct an exchange-permuted tensor, or form any physical four-index
object.

For contracted occupied matrices `ca` and `cb`,

```text
Dqa = ca*ca'
Dqb = cb*cb'

Fqa[a,b] = Hqa[a,b]
           + sum_cd G[a,b,c,d]*(Dqa[c,d] + Dqb[c,d])
           - sum_cd G[a,c,d,b]*Dqa[c,d]

Fqb[a,b] = Hqb[a,b]
           + sum_cd G[a,b,c,d]*(Dqa[c,d] + Dqb[c,d])
           - sum_cd G[a,c,d,b]*Dqb[c,d]

Eq = 0.5*dot(Dqa,Fqa+Hqa) + 0.5*dot(Dqb,Fqb+Hqb).
```

Symmetrize each Fock only at its eigensolve boundary.  Compact oracles must
compare both spin Focks and the total energy against independent physical
contractions for common-H and split-H cases.

The accepted pair-map/model guard remains

```text
2*N*r^2 + r^4 <= N^2
```

in `Float64` elements, evaluated with checked integer arithmetic before
allocating `P`.  UHF adds only one optional `r^2` split-H one-body projection
and a second set of small DIIS arrays.  The inequality continues to describe
the dominant pair-map/model construction, not total heap ownership.  A guard
failure preserves both incoming spin determinants and proceeds to ordinary
cleanup.

## Coupled UHF DIIS

Seed the contracted state with

```text
ca0 = Q'*Ca
cb0 = Q'*Cb.
```

At every iteration compute separate `Fqa`, `Fqb`, `Kqa`, and `Kqb`.  The Pulay
Gram matrix uses their combined Frobenius product:

```text
B[i,j] = dot(Kqa_i,Kqa_j) + dot(Kqb_i,Kqb_j),
```

omitting an empty spin.  Equivalently, the existing rank-revealing Pulay
solver can receive horizontally paired Fock and commutator matrices.  One
coefficient vector mixes both Focks, preserving their coupled UHF history.
Never run independent alpha and beta DIIS coefficient solves.

After mixing, diagonalize the two small symmetric Focks separately and occupy
their respective lowest `na` and `nb` eigenvectors.  The contracted convergence
gate requires each nonempty spin to satisfy

```text
norm(Kqs) <= 1e-10*max(1,norm(Fqs)).
```

The best iterate is seeded with the incoming pair `(ca0,cb0)`.  Lower total
energy wins outside the accepted `128eps(Float64)` energy tie; inside the tie,
the smaller RMS two-spin commutator wins.  The Pulay affine constraint,
rank-revealing solve, coefficient bound, two-consecutive-failure rejection,
DIIS memory, and iteration cap remain exactly those of the RHF driver.

Any nonfinite field, eigensolve result, coefficient set, or spin determinant
rejects the complete paired proposal.  No alpha-only or beta-only state may
escape from the contracted solver.

## Full-Basis Atomic Acceptance And Fallback

Lift both trial determinants:

```text
Catrial = Q*ca
Cbtrial = Q*cb.
```

Independently reconstruct their physical UHF Focks, energy, residuals, Gram
errors, `S2`, and spin-labeled occupied overlaps.  Do not diagonalize the full
physical Focks or compute Aufbau gaps during proposal acceptance.  Accept the
pair only if:

- all spin orbitals and audit quantities are finite;
- both occupation counts are unchanged;
- each nonempty spin passes its orthogonality gate;
- the total energy satisfies `Etrial <= Ebefore + tau_E`;
- each nonempty spin has
  `minimum(svdvals(Cs_before'*Cs_trial)) >= 0.99`; and
- the contracted/direct spin-Fock and energy orientations have passed the
  independent compact oracle.

Residual growth and `S2` movement are recorded but do not reject a downhill
proposal.  A common-H calculation may legitimately break or restore spin
symmetry.

The acceptance routine prepares both spin decisions and diagnostics without
mutating the queue or current pair.  It then commits both or neither.  On
rejection, numerical failure, no response, or resource fallback, ordinary
cleanup begins from the unchanged incoming `(Ca,Cb)`.  Only the complete
cleanup endpoint is appended to the FIFO.

This atomic boundary must be fault-tested with an alpha-valid/beta-invalid
trial, verifying unchanged determinants, queue, counters, and cleanup start.

## Sweep Composition And Route Dispatch

The fixed RHF policy is reused without tuning:

```text
bootstrap sweeps       6
captures               2,4,6
FIFO length            6 paired determinants
cleanup sweeps         2
maximum cycles         25
target                 max spin orbital residual <= 1e-6
DIIS iterations        40
DIIS memory            7 paired Focks/errors
minimum overlap        0.99 per nonempty spin
```

The common-H private route calls existing
`solve_hfdmrg(H,V,Ca,Cb; ...)`.  The split-H route calls existing
`solve_hfdmrg(Ha,Hb,V,Ca,Cb; ...)`.  Both use `cutoff=0`, complete-sweep
observers, and fresh two-sweep cleanup segments.  The existing `Ha === Hb`
fast-path behavior remains authoritative.

The driver forwards only block size/partition, center width, local SCF cutoff,
environment cutoff, and verbosity.  It does not forward a user observer or
alter ordinary convergence semantics.

## Return Semantics

The private result contains:

- `Cup` and `Cdn` from the last complete sweep segment;
- the independently recomputed physical determinant energy;
- alpha, beta, maximum, and RMS orbital residuals;
- Gram, `S2`, and overlap diagnostics; and
- terminal-only Aufbau/gap diagnostics, timed separately from acceleration;
- target/terminal status, sweep/cycle counts, and atomic event counts;
- paired-history ranks and DIIS status; and
- phase timing and allocated-byte summaries.

The returned energy is never a damped local-SCF scalar or a contracted model
surrogate.  No public `solve_hfdmrg` return convention changes.

## Cost And Memory Model

Let `nmax=max(na,nb)`, `q` be the number of nonempty spins, `h` the number of
older history pairs, and `r` the common history rank.

| Phase | Work | Additional memory |
| --- | --- | --- |
| Paired history queue | none between captures | `O(h*N*(na+nb))` retained |
| Current spin-union SVD | `O(N*(na+nb)^2)` | `O(N*(na+nb))` |
| Paired response SVD | `O(N*(q*h*nmax)^2)` | `O(N*q*h*nmax)` |
| One-body projection | common `O(N^2*r)`, split `O(2N^2*r)` | `O(N*r+r^2)` |
| Shared pair map/model | `O(N^2*r^2+N*r^4)` | `O(N*r^2+r^4)` retained |
| Coupled DIIS iteration | `O(r^4+r^3)` | `O(r^4+2m*r^2)` |
| Full paired proposal audit | `O(N^2*(na+nb))` | several `O(N^2)` temporaries |
| Terminal Aufbau audit | `O(N^3)` once, outside acceleration timing | `O(N^2)` temporary |

UHF does not duplicate the dominant `P` or `G`.  Relative to RHF at the same
`N` and `r`, its contracted iterations approximately double small Fock,
commutator, and eigensolve work, while interaction setup remains shared.

The method is attractive only when the common rank remains much smaller than
`N`.  Strongly different alpha/beta histories can increase the anchor and
response ranks enough to trigger the resource guard.  That is a safe ordinary
fallback, not permission to trim one spin more aggressively after observing a
fixture.

## Smallest Implementation Seam

The closed implementation file set is:

- `src/HFDMRG.jl`: one private include;
- new `src/history_accelerated_uhf.jl`: paired basis, physical/model metrics,
  coupled solve, atomic acceptance, segment dispatch, and private driver; and
- `test/runtests.jl`: compact numerical and workflow coverage.

The new file may reuse the private RHF policy values, checked workspace guard,
and rank-revealing Pulay implementation.  It must not turn them into a generic
backend accelerator or change accepted RHF arithmetic.  Archival references to
`_HistoryRHFPolicy` remain valid; public naming cleanup is deferred.

Do not edit `src/core.jl`, `src/backend_api.jl`, density backend files, README,
or manifests.  The conceptual entry points are private methods:

```text
_solve_history_uhf(H,V,Ca,Cb; ...)
_solve_history_uhf(Ha,Hb,V,Ca,Cb; ...).
```

No symbol is exported.

Preferred addition budget:

```text
src/history_accelerated_uhf.jl   <= 215 lines
src/HFDMRG.jl                    1 line
test/runtests.jl                 <= 115 lines
total                            <= 331 lines
hard maximum                     380 lines
```

If the atomic two-spin contract cannot fit within 380 added lines, stop rather
than widening the design.  No accepted RHF code is expected to be deleted in
the first port.  After UHF is independently accepted, exact duplicate scalar
bookkeeping may be considered for shrinkage in a separate review.

## Validation Ladder

### 1. Compact algebra and invariance

- Common-H and split-H alpha/beta contracted Focks versus independent physical
  projections, each `<=2e-12`.
- Contracted versus physical total energy `<=2e-12`.
- RHF-limit parity when `Ca==Cb`, including common basis rank, energy,
  residual, and proposal projector.
- Independent alpha/beta occupied rotations and sign changes leave selected
  rank, basis projector, energy, and both proposed projectors unchanged within
  `2e-12`.
- Deterministic covariance spectra straddle the numerical, absolute, and tail
  rules for paired FIFO lengths three, four, and six.
- Full-rank and deficient coupled Pulay systems enforce one affine coefficient
  vector to roundoff.

Tests protect numerical results and determinant behavior, not private buffer
shapes or helper vocabulary.

### 2. Atomic workflow controls

- Common-H and split-H six-sweep capture followed by at least two complete
  history/cleanup cycles.
- At least one accepted paired proposal with full Fock/energy parity.
- Alpha-valid/beta-invalid rejection with neither spin or queue mutated before
  ordinary cleanup.
- Independent rotated-history workflow with identical physical endpoints.
- Zero-response no-op and pair-model resource fallback.
- One-empty-spin common-H and split-H controls.
- Exact ordinary-cleanup equality after rejection or resource fallback.
- No public symbol or keyword appears.

### 3. Package and route regression

- Complete package suite.
- Existing cached/projection verifier.
- Exact public default parity for representative density RHF, common-H UHF,
  split-H UHF, target-residual, projection-sliced, and cached-sliced routes.
- Existing private RHF history workflow remains numerically unchanged.
- `git diff --check` and exact changed-file/line-count review.

### 4. Broken-symmetry and basin validation

Before promotion beyond an experimental private checkpoint, use a deterministic
density-density family with a genuinely spin-broken endpoint.  Run common-H
and a related split-H control from identical paired seeds.  Require:

- the fixed history policy, with no cutoff or cycle tuning;
- final maximum spin residual `<=1e-6`;
- positive occupied--virtual gaps for both nonempty spins;
- finite Gram, energy, `S2`, and spin-labeled projector diagnostics;
- no atomicity or resource-integrity failure; and
- external full-space UHF commutator polishing of accelerated and matched
  ordinary endpoints until each spin has `norm(Ks)<=2e-11`.

After polishing, require total-energy agreement `<=2e-9 Ha` and each
spin-labeled projector distance `<=2e-6`.  Compare `S2` after closure and report
any alpha/beta-swapped symmetry partner explicitly rather than hiding it in a
spin-averaged projector.

This gate establishes the same attraction basin, not the global UHF minimum.
Its polishing time is separate from driver timing.

### 5. Engineering evidence

On the broken-symmetry fixture report ordinary and accelerated sweeps, cycles,
common ranks, accepted/rejected/fallback events, model/DIIS/audit times,
allocation, and RSS.  The first port has a correctness gate but no universal
speed claim.  A failure to accelerate is reported without policy tuning.

No H10 calculation belongs here: production H10 uses a sliced interaction and
cannot exercise this density-only port.

## Promotion Boundary And Roadmap

Passing compact tests, the package suite, and default parity permits one
private experimental implementation commit.  Promotion for scientific use
also requires the independent broken-symmetry/basin gate.

Only after that review may work proceed to sliced interactions.  The intended
sequence is:

1. density-density UHF history parity;
2. independent broken-symmetry and stationary-basin validation; and
3. projection-sliced and cached-sliced UHF history support under a separate
   backend-aware design.

The density port alone does not accelerate production H10.  Sliced support
must preserve each backend's compact contraction and may not materialize a
physical global four-index interaction merely to reuse this implementation.

## What This Design Does Not Prove

The algebra probe and RHF evidence do not prove that UHF history acceleration
is faster, that it preserves every metastable spin basin, or that a common
spin-history basis remains small for strongly polarized systems.  They also do
not establish sliced-backend support, cache reuse, public policy, or a new
default center geometry.  Those remain explicit later gates.
