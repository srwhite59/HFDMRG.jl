# Frozen-Occupied Environment Design

Date: 2026-07-30

Status: amended F0 design checkpoint for paper-manager review. F1 is not
authorized. This does not authorize source/test edits, an API addition,
threshold freezing, H100, merge, or push.

## Outcome

Exactly completed occupied orbitals can leave the moving environment only if
they remain individually recoverable. The proposed solver therefore:

1. classifies completed alpha and beta combinations separately;
2. contracts their Fock fields and constant energy, but stores each physical
   frozen column once in a directional ledger; and
3. audits every frozen column against the physical Fock matrix after each
   sweep, rebuilding environments if any column must thaw.

Ordinary backend caches contain only the common spatial span of the remaining
fractional occupied directions. This is one opt-in solver policy, not generic
backend composition. F1 is roundoff-exact only. Occupation leakage is not an
energy bound; threshold-controlled freezing is outside this design.

The common spatial basis is not itself an allowed orbital space for either
spin. Every local solve uses a spin-specific null space that excludes that
spin's frozen columns while leaving the other spin free to share them.

## Baseline And Evidence

The inspected repository state was:
```text
branch  main
HEAD    3f5a11da970a68efcf6923f4ef65d4b51546c5a2
tree    6d9b9662d8a8b089c57ff85b2d8944c498ada052
relation to paper baseline: 0 ahead, 0 behind
status:
 M AGENTS.md
?? JuliaStyle.md
```
Those unrelated files remain excluded.

The accepted six-function oracle found maximum errors of `6.66e-16 Ha` in
energy, `2.29e-16` in Fock Frobenius norm, `2.22e-16 Ha` in the recursive
scalar, and `1.36e-16` in the recursive field. Stationary frozen residuals
were at most `2.29e-16`; a forced `2.5e-5` coupling was recovered within
`1.76e-18`. Its tail study showed errors approximately linear in tail
amplitude although leakage is quadratic.

## Algebra And Energy

For spin `sigma`,
```text
D_sigma = D_f,sigma + D_a,sigma
G_alpha(D) = J(D_alpha + D_beta) - K(D_alpha)
G_beta(D)  = J(D_alpha + D_beta) - K(D_beta)
h_eff,sigma = h_sigma + G_sigma(D_f)

E_f = sum_sigma Tr(D_f,sigma h_sigma)
    + 1/2 sum_sigma Tr[D_f,sigma G_sigma(D_f)]

E = E_f + sum_sigma Tr(D_a,sigma h_eff,sigma)
  + 1/2 sum_sigma Tr[D_a,sigma G_sigma(D_a)]

F_a,sigma = h_eff,sigma + G_sigma(D_a).
```
For RHF, alpha and beta projectors are identical and promotion is paired:
```text
h_eff = h + 2J(D_f) - K(D_f)
E_f = 2Tr(D_f h) + 2Tr[D_f J(D_f)] - Tr[D_f K(D_f)].
```
The core must use `h_eff` in both positions of its trace formula:
```text
E_UHF = E_f + 1/2 sum_sigma Tr[D_a,sigma(F_a,sigma+h_eff,sigma)]
E_RHF = E_f + Tr[D_a(F_a+h_eff)].
```
Adding the frozen field only to `F` would lose half the active-frozen cross
energy.

At each window form one effective one-body matrix per spin,
```text
H1B_eff,sigma = H1B_sigma + G_frozen,window,sigma.
```
That same matrix seeds every candidate Fock and occupies the second position
of the energy trace. Frozen fields do not also enter the active interaction
backend. Add the left and right frozen self energies and their cross scalar
once. Even when the bare `H` is common, UHF frozen exchange generally makes
`H1B_eff,alpha != H1B_eff,beta`; exact common-H UHF therefore takes the
existing split-matrix local route internally without changing default mode.

For disjoint left/right frozen densities, use
```text
E_f,LR = E_f,L + E_f,R
       + 1/2 sum_sigma {
           Tr[D_f,L,sigma G_sigma(D_f,R)]
         + Tr[D_f,R,sigma G_sigma(D_f,L)] }.
```
Both orientations remain explicit. Exact-mode sliced input must pass fixed
internal physical-symmetry checks so this scalar and `h_eff` are reciprocal;
`W` and `Wswap` nevertheless stay distinct.

## Spin-Resolved Classification

Let `C_a,sigma` be active occupied columns and `P_E` the completed physical
environment. Diagonalize separately by spin:
```text
M_sigma = C_a,sigma' P_E C_a,sigma
M_sigma V_sigma = V_sigma diag(lambda_sigma).
```
Candidates are `C_a,sigma V_sigma`. A combined `D_alpha+D_beta` spectrum is
never used because it can exceed one and cannot identify spin exchange.

For element type `T`, define
```text
gamma(m) = m eps(T)/(1-m eps(T))
tau_lambda = 128 gamma(length(E)) max(1,n_a,sigma)
tau_tail = 128 eps(T) max(1,n_a,sigma) sqrt(max(1,length(Ebar)))
           max(1,norm(C_a,sigma,Inf)).
```
Exact mode requires `length(E) eps(T) < 1/2`. An eigenvalue outside
`[-tau_lambda,1+tau_lambda]` is resolvably invalid and stops the solve. Values
inside that envelope may be clamped for classification. A unit candidate
must pass both:
```text
1-lambda_clamped <= tau_lambda
norm(P_Ebar C_a,sigma v) <= tau_tail.
```
The explicit tail norm prevents cancellation in `1-lambda`. Store normalized
`P_E C_a,sigma v` only after its change from the candidate passes the same
scale. These are fixed internal checks, not public tolerances, and do not
reuse `environment_cutoff`.

A direction thawed by the residual audit stays protected from promotion for
the rest of the solve. Transport it by projection into the current occupied
subspace; a candidate with roundoff-unit overlap with that protected subspace
remains active.

After removing certified modes, form
```text
X_E = [P_E C_remaining,alpha  P_E C_remaining,beta].
```
Apply the existing SVD and strict `sigma > environment_cutoff` rule to `X_E`.
Thus that public control still compresses the remaining occupied span; it is
not a freezing or energy threshold. If no direction passes, choose the first
deterministic coordinate direction in the orthogonal complement of the
frozen spatial span. If the complement is empty, leave one certified
candidate active instead of constructing a zero-rank block.

Fixed edges perform no promotion and retain today's complete identities.
Promotion begins only after an edge absorbs its first center chunk.

## Spin-Specific Allowed Coordinates

For the physical window basis `B` and frozen columns of spin `sigma`, compute
an orthonormal null-space basis
```text
A_sigma = C_f,sigma' B
Z_sigma = null(A_sigma)
C_f,sigma' B Z_sigma = 0.
```
For `w=size(B,2)` and `f_sigma=size(C_f,sigma,2)`, take a full SVD of
`A_sigma` and define
```text
tau_Z = 128 gamma(max(1,w,f_sigma)) max(1,norm(A_sigma,2)).
```
Singular directions strictly above `tau_Z` are constrained; the remaining
right singular directions form `Z_sigma`. Certify
`norm(A_sigma Z_sigma,2)<=tau_Z` and orthonormal `Z_sigma` at the corresponding
fixed scale. This is not `environment_cutoff`; empty `C_f,sigma` gives `Z=I`.
Require `n_a,sigma <= size(Z_sigma,2)`. The spin solve is
```text
F_allowed,sigma = Z_sigma' F_window,sigma Z_sigma,
C_a,window,sigma = Z_sigma U_occ,sigma.
```
Thus alpha cannot duplicate an alpha-ledger orbital even when beta keeps that
spatial direction in the common basis; beta is constrained only against its
own frozen columns.

All candidate and damped occupied projectors are formed in `Z_sigma`
coordinates. Active densities use `B Z_sigma U_occ,sigma`; energy uses those
same densities. After promotion, structural demotion, or thaw, recompute
`Z_sigma`, reduce `D_sigma-D_f,sigma` into `B Z_sigma`, and reorthogonalize
there before the next update. Reconstruction lifts through `B Z_sigma`.
Post-update orthogonality is an invariant check, not the enforcement method.

`n_a,sigma=0` is valid independently for either spin. Represent its occupied
coefficients as `size(Z_sigma,2) x 0`, set its active density to zero, skip its
eigensolve and damping, and retain its effective Fock for the other spin and
the sweep audit. If both UHF spins, or all RHF pairs, are frozen, skip the
local occupied solve and evaluate the frozen-only scalar; local convergence
is vacuous but outer convergence still requires the energy and sweep audits.
Structural demotion may restore a positive active count at the next cut.
The rank-one environment fallback remains a spatial basis rule and never
creates a rank-one occupied result.

## State And Cut Lifecycle

The opt-in path adds private state:
```text
FrozenLedger
    C_alpha, C_beta       physical columns stored once per side generation
    spin/promotion order
    protected thaw subspaces

FrozenBlockState
    side, physical range
    alpha/beta ledger references and counts
    E_f                  self energy of this block's frozen density
    backend-specific frozen field
```
`LRBlock.phi` and its ordinary one-body/backend caches describe only the
active basis. A cut references immutable ledger columns rather than copying
an `N x f` projector.

At every window reconstruct
```text
C_sigma = [C_f,L,sigma, B_a C_a,window,sigma, C_f,R,sigma]
n_a,sigma = N_sigma-f_L,sigma-f_R,sigma.
```
Require orthonormality, left/right frozen orthogonality, and containment of
both frozen projectors in `D_sigma` at fixed scale-aware tolerances. After a
move, promotion or structural demotion can change `n_a,sigma`; reconstruct
`D_a=D-D_f,L-D_f,R`, choose an orthonormal active gauge, and reduce it into
the new window. This avoids stale column order and restores a column when a
smaller saved block no longer classifies it as frozen.

Each immutable ledger generation belongs to one side and one growth pass.
Cuts created by that pass reference prefixes of its append-only column list;
fixed edge cuts have empty references. At the start of a left-growing pass,
create a new left generation. Each newly written left cut drops its reference
to the old left generation and takes a prefix reference to the new one; right
cuts continue to reference the retained right generation. Release the old
left generation only after its final cut and any in-flight window reference
are gone. Apply the mirror rule at reversal for a new right generation.

Consequently an old same-side, new same-side, and opposite-side generation
may coexist during overwrite; the steady completed direction has two. No
reference migrates without reconstructing that cut, and no cut refers to a
mutable aggregate density. A thaw first materializes the full determinant
while both old directional generations are valid, invalidates all their cuts,
then releases them before fresh initialization. Observer data owns its
expanded determinant and cannot extend a ledger lifetime.

## Promotion And Final Coordinates

For all projectors promoted in one absorption, define
`G_P,sigma=G_sigma(P_alpha,P_beta)`. Both growth directions use:
```text
h_eff,sigma' = h_eff,sigma + G_P,sigma

E_f' = E_f
     + Tr(P_alpha h_eff,alpha) + Tr(P_beta h_eff,beta)
     + 1/2 Tr(P_alpha G_P,alpha)
     + 1/2 Tr(P_beta G_P,beta).
```
Simultaneous alpha/beta promotion is mandatory because it contains their
promoted-promoted Coulomb cross term.

Preserve the current final-basis computation. For preliminary maps `A0,C0`
and reorthogonalization `R`,
```text
left:  phi_new = [phi_old A0; C0]R
right: phi_new = [C0; phi_old A0]R
A=A0R, C=C0R.
```
Every one-body, interaction, and frozen-field recurrence receives final
`A,C`, never `A0,C0`. In particular,
```text
Kbar_f,sigma' =
    A' Kbar_f,sigma A + phi_new' K_sigma(P_sigma) phi_new.
```
Left/right differ only in physical concatenation and exterior range.

## Thaw

After initial classification and every complete sweep, construct physical
densities and Fock matrices through the backend's existing full contraction.
For every frozen column, audit stationarity:
```text
r_f,sigma = (I-D_sigma)F_sigma c_f,sigma
           = F_sigma c_f,sigma
           - C_sigma(C_sigma'F_sigma c_f,sigma).
```
Use the fixed internal gate
```text
gamma_N = N eps(T)/(1-N eps(T))
tau_thaw,sigma = 128 gamma_N max(1,norm(F_sigma,Inf)).
```
Strictly thaw when `norm(r_f,sigma)>tau_thaw,sigma`. A triggered audit
prevents convergence return and forces another sweep.

Residual stationarity is insufficient for occupation ordering. Let
`C_sigma` contain every occupied orbital and let `U_sigma` span its physical
orthogonal complement. Define
```text
epsilon_occ,max = lambda_max(C_sigma' F_sigma C_sigma)
epsilon_vir,min = lambda_min(U_sigma' F_sigma U_sigma)
Delta_Aufbau = epsilon_occ,max-epsilon_vir,min
tau_Aufbau = 128 gamma_N max(1,norm(F_sigma,Inf)).
```
The audit is vacuous when either the occupied or virtual space is empty.
`Delta_Aufbau>tau_Aufbau` is a resolved occupied-virtual inversion even if
the frozen residual is zero. It prevents convergence and deterministically
thaws all frozen orbitals of that spin; those columns are protected from
re-promotion for the rest of the solve. Values within the fixed tolerance
are roundoff-degenerate and pass. If a resolved inversion remains after that
spin has no frozen columns, exact mode reports unresolved non-Aufbau
stationarity rather than claiming convergence or loosening the gate.

Thaw rebuilds rather than mutating all saved cuts:

1. protect the triggered occupied continuation;
2. reconstruct the same full determinant;
3. discard both directional block families and contracted frozen fields;
4. rebuild from that determinant; and
5. continue with the column active.

The ledger and original backend data make this reversible. The inverse scalar
identity is
```text
h_eff,sigma' = h_eff,sigma-G_P,sigma
E_f' = E_f
     - Tr(P_alpha h_eff,alpha) - Tr(P_beta h_eff,beta)
     + 1/2 Tr(P_alpha G_P,alpha)
     + 1/2 Tr(P_beta G_P,beta),
```
but production rebuilds so every cut agrees. Repeated rebuilds are a measured
stop condition, not authority to loosen the threshold.

## Density-Density State

On entry to exact mode, first reject any nonfinite `V`. Exhaustively certify
all stored pairs with the maximum-entry norm:
```text
s_V = max(1,maximum(abs,V))
delta_V = max(p,q) abs(V[p,q]-V[q,p])
tau_sym = 128 eps(T) s_V.
```
Require `delta_V<=tau_sym`; do not symmetrize the input. This certification is
specific to `:roundoff_exact` and leaves default density arithmetic unchanged.

For block range `R`, store
```text
j_f = V[:,R] diag(D_f,alpha+D_f,beta)                 length N
Kbar_f,sigma = phi'(V[R,R].*D_f,sigma)phi
E_f, counts, ledger references.
```
For `P=P_alpha+P_beta`,
```text
j_f' = j_f + V[:,R_new]diag(P)
Kbar_f,sigma' = A'Kbar_f,sigma A
              + phi_new'(V[R_new,R_new].*P_sigma)phi_new.
```
At a window, project diagonal `j_f,L+j_f,R` into left/center/right active
bases and subtract each side's `Kbar_f,sigma` on its active block. Exchange
cross is zero for disjoint support. Evaluate the scalar cross in both direct
orientations:
```text
E_cross = 1/2 {
  diag(D_f,L,total)'j_f,R[L] + diag(D_f,R,total)'j_f,L[R] }.
```
Reconstruct those diagonals from ledger columns instead of copying charge at
every cut. Thaw uses
```text
F_sigma = h_sigma + Diagonal(V*diag(D_alpha+D_beta)) - V.*D_sigma.
```
## Cached-Sliced State

Use the unsymmetrized orientation
`V_nm[a,b,c,d]=V6[a,b,c,d,n,m]`. For
```text
T_f,m[c,d]=(D_f,alpha+D_f,beta)[(m,c),(m,d)],
```
store every output-slice direct field
```text
J_f,n[a,b] = sum(m in block,c,d) V_nm[a,b,c,d]T_f,m[c,d],
```
and spin exchange projected into the active block:
```text
K_f,sigma[(n,a),(m,d)] =
    sum(b,c) V_nm[a,b,c,d]D_f,sigma[(m,c),(n,b)]
Kbar_f,sigma=phi'K_f,sigma phi.
```
Pack `J_f` with existing ragged offsets. Promotion adds its ordered
`V(output,promoted)` channel and transforms/adds exchange with final maps. No
`V_nm=V_mn` equality is used in these recurrences.

Candidate B remains `G` plus distinct exterior `W/Wswap`. Promoted active
coordinates disappear from those arrays; their effects collapse into
`J_f,Kbar_f,alpha,Kbar_f,beta,E_f`. Physical columns, spin, and order remain
in the ledger, and the global sliced interaction remains referenced, so thaw
can rebuild exactly. Never merge `W` and `Wswap`.

Aligned whole-slice windows add `J_f,L+J_f,R` per slice and both projected
exchange blocks. Compute the frozen cross scalar in both ordered direct
orientations. Exact mode requires `block_partition==layout`; split-slice
fallback is rejected. For every valid stored real entry, exact mode rejects
nonfinite data and exhaustively checks
```text
V_nm[a,b,c,d] = V_nm[b,a,c,d]
V_nm[a,b,c,d] = V_nm[a,b,d,c]
V_nm[a,b,c,d] = V_mn[c,d,a,b].
```
With `s_6=max(1,maximum(abs,V6))`, let `delta_6` be the maximum-entry norm of
all three difference tensors and require
`delta_6<=128 eps(T)s_6`. Do not average or repair entries. Default cached
arithmetic continues to accept ordered unsymmetrized data, and certification
does not permit merging `W` with `Wswap`.

## Cost And Lifecycle

Let `S` be slices, `N=Sd`, `q` common fractional rank, `f` frozen count, and
`C=Theta(S)` retained cuts. Cached active all-cut storage remains
```text
O(CNq + Cq^4 + CSd^2q^2).
```
Frozen additions at steady completed directions are
```text
ledgers                   O(Nf)
all-cut direct fields     O(CSd^2)
projected exchanges       O(Cq^2)
counts/scalars/refs       O(C).
```
Total cached retained storage is
```text
O(CNq+Cq^4+CSd^2q^2+Nf+CSd^2+Cq^2).
```
At fixed `d,q`, retained state and supplied sliced interaction are quadratic
in `S`; this is not a linear-scaling claim. A window adds worst-case
`O(Sd^2q^2)` frozen direct-field projection work. Promotion contracts ordered
direct channels and exchange; thaw causes a rare rebuild. The sweep audit
uses two `N x N` Fock temporaries already supported by density or sliced
projection, never a global `N^4` interaction. Window scratch is `O(w^2)` and
steady local Fock addition must allocate zero.

Density-density frozen storage is `O(CN+Cq^2+Nf)`, worst-case window
projection `O(Nq^2)`, and audit `O(N^2)`.

During a direction overwrite, three ledger generations can coexist, so
physical frozen-column storage peaks at at most three directional
generations, `O(3Nf)=O(Nf)`, rather than the steady two-generation constant.
The all-cut field terms are unchanged. Thaw additionally holds one materialized
`N x (N_alpha+N_beta)` determinant while old generations are released.

H100 estimates must use measured cut-by-cut `q,f`. The preflight's
leakage-`1e-6` rank near five per spin is not the roundoff-exact classifier
and is not a production prediction.

## Production Seam

After F1 validation, add one keyword to every public route:
```text
frozen_occupied=:off
```
Accept `:off` and `:roundoff_exact`. Omission and explicit `:off` enter the
existing path before policy allocation or validation, preserving exact
trajectories. Exact mode initially supports density-density and aligned
cached-sliced RHF, common-H UHF, and split-H UHF. Projection sliced and
nonzero target-residual routes reject it clearly. An exactly zero residual
still routes to density-density and may use it.

Threshold freezing will not share this keyword or `environment_cutoff`; it
requires separate scientific design and error reporting.

Private dispatch covers frozen initialization, promotion, window field/scalar,
and physical Fock audit. The five active-interaction operations remain
unchanged. Do not compose generic backends or copy the sweep.

## Milestones, Files, And Budget

F1 implements density-density classification, dynamic active occupations,
determinant reassembly, RHF/UHF promotion, fields, energy, audit, protected
thaw, and rebuild. Allowed files:
```text
src/HFDMRG.jl
src/core.jl
src/backend_api.jl
src/frozen_occupied.jl
src/backends/density_density.jl
test/runtests.jl
README.md
docs/backend_architecture.md
```
Preferred additions: source 330, tests 150, docs/include 60, total 540. Hard
stop: 650 total. Do not touch sliced cache arithmetic in F1.
F2 starts only after F1 review and adds aligned cached-sliced exact mode in:
```text
src/backends/sliced_basis_cached.jl
src/frozen_occupied.jl
test/runtests.jl
README.md
docs/backend_architecture.md
scripts/verify_cached_vs_projection.jl  only if required
```
Preferred addition: 380 lines; hard stop: 480. No accepted code is obsolete
at F0. Share classification, reassembly, and energy code across common/split
paths; do not commit scratch oracle translations or private-shape tests.

## Validation Ladder

Before edits, freeze the current default trajectories and cached fixed/ragged
verifier outputs.

1. Translate the independent six-function density and sliced-pair oracle into
   compact Julia numerical tests.
2. Require explicit four-index alpha/beta Fock and energy parity before/after
   recursive promotion.
3. Require spin-swapped promotion plus initial/promoted thaw controls.
4. Add simultaneous nonzero alpha/beta promotion and its Coulomb cross term.
5. Require omission versus `:off` exact equality for RHF, common/split UHF,
   density, projection/cached sliced, and zero/nonzero target-residual routes.
6. In moving sweeps require roundoff no-promotion, exact left/right promotion,
   structural demotion, and forced thaw/rebuild with invariant determinant.
7. Compare density and cached left/right/rank-changing recurrences with fresh
   physical contractions, including post-reorth maps, fixed/ragged slices,
   and zero steady Fock allocation.
8. Run the full package suite and maintained cached/projection verifier.
9. Require H10 endpoint energy, occupied projectors, physical Fock residual,
   spin diagnostics, and route parity before claiming compression.
10. Only then recompute H20/H40/H80/H100 cut-based time/memory estimates. Do
    not run H100 before review.
11. In UHF, let beta require a spatial direction frozen for alpha and prove
    spin-specific allowed spaces prevent only alpha duplication.
12. Give a frozen exact Fock eigenvector a resolved lower virtual; require
    zero residual but an Aufbau-triggered thaw/rebuild.
13. Exercise common-H UHF with unequal alpha/beta frozen exchange fields.
14. Require one-spin-zero-active and both-spins-zero-active UHF plus
    all-frozen RHF, including frozen-only energy and structural demotion.
15. Hold old/new same-side and opposite-side ledger generations live through
    reversal, then thaw and prove deterministic release and reconstruction.
16. Accept symmetric and reject nonfinite or resolvably nonsymmetric density
    and sliced inputs only in exact mode.
17. Independently decompose energy and prove each frozen self scalar,
    left-right cross scalar, and active-frozen field enters exactly once.

## Risks And Stop Conditions

- A field without recoverable physical frozen columns is irreversible: stop.
- A cut whose frozen projector is not contained in the live determinant is
  stale: stop; do not loosen tolerances.
- Preliminary rather than post-reorth maps repeats the repaired coordinate
  defect: stop.
- Repeated thaw or active rank returning to occupied count is a valid negative
  result, not authority for threshold mode.
- Failed sliced physical symmetry rejects exact mode; it does not authorize
  merged orientations or changed default semantics.
- A frozen residual pass with a resolved Aufbau inversion triggers thaw; it
  cannot be accepted as stationary.
- Spin-specific null-space failure or an active count exceeding its allowed
  dimension stops exact mode.
- Any default-off trajectory difference blocks the feature.
- Exceeding a milestone hard cap requires new review.
- H10 failure blocks long-chain calculations.

## F0 Decision

This amended design needs neither a global four-index interaction nor a copied
sweep engine. It is not approved for implementation. F1 requires
paper-manager review of the spin-specific allowed spaces, residual-plus-Aufbau
audit, zero-active behavior, interaction certification, reversible lifecycle,
energy convention, and bounded private seam.

-- hfdmrg-manager@rh310l
