# Roundoff-Exact Frozen Occupied Environments: Final Contract

Date: 2026-08-05
Accepted implementation commits: `313149b`, `9d8d510`, `4164cf8`
Status: supported density-only correctness mode; not a recommended optimization

## Outcome And Scope

`frozen_occupied=:roundoff_exact` allows an occupied orbital that is already
complete inside a grown environment to leave the active environment basis
while retaining its exact density contribution, Fock field, and energy. The
physical frozen column remains recoverable in a directional ledger. Sweep-end
audits thaw and protect a direction when it is no longer stationary or no
longer occupation ordered.

The mode supports density-density RHF, common-one-body UHF, split-one-body UHF,
and the exact-zero target-residual route that delegates to density-density.
Nonzero target-residual and sliced backends reject it. Omission and explicit
`:off` use the unchanged ordinary path and are required to agree exactly.

This is a correctness-oriented roundoff classification, not a user-selected
leakage threshold. It has not demonstrated a general performance benefit and
is not recommended as the default solver policy.

## Density, Fock, And Energy Algebra

For each spin, decompose the pure occupied determinant density as

```text
D_sigma = D_f,sigma + D_a,sigma.
```

With density-density direct and exchange maps,

```text
G_alpha(D) = J(Dalpha+Dbeta)-K(Dalpha)
G_beta(D)  = J(Dalpha+Dbeta)-K(Dbeta)
h_eff,sigma = h_sigma + G_sigma(D_f).
```

The frozen scalar and active energy are

```text
E_f = sum_sigma Tr(D_f,sigma h_sigma)
    + 1/2 sum_sigma Tr[D_f,sigma G_sigma(D_f)]

E = E_f + sum_sigma Tr(D_a,sigma h_eff,sigma)
  + 1/2 sum_sigma Tr[D_a,sigma G_sigma(D_a)].
```

Equivalently, the existing local trace formulas use `h_eff` in both positions:

```text
E_UHF = E_f + 1/2 sum_sigma Tr[D_a,sigma(F_a,sigma+h_eff,sigma)]
E_RHF = E_f + Tr[D_a(F_a+h_eff)].
```

Adding a frozen field only to `F` would lose half the active-frozen cross
energy. Common-H UHF still needs separate alpha and beta effective one-body
matrices because frozen exchange is spin specific; it remains in the common-H
sweep core rather than copying or switching sweep engines.

The determinant stored and carried between windows is always the pure
idempotent determinant. Damping applies only to the temporary local SCF density.
Frozen bookkeeping must never replace the carried determinant with a damped,
non-idempotent density.

## Classification And Allowed Spaces

For a spin-specific active occupied matrix `C_a` and completed environment
projector `P_E`, diagonalize

```text
M = C_a' P_E C_a.
```

A rotated occupied direction is frozen only when its eigenvalue is
roundoff-equal to one and its explicit norm outside the environment is also at
roundoff scale. The implementation uses fixed, scale-aware tolerances derived
from `eps(T)`, dimensions, and matrix scales. These are not public keywords and
are distinct from `environment_cutoff`, which continues to compress the
remaining occupied environment span.

Alpha and beta are classified independently. In a physical window basis `B`,
each spin solves in its own allowed complement:

```text
A_sigma = C_f,sigma' B
Z_sigma = null(A_sigma)
C_f,sigma' B Z_sigma = 0.
```

This prevents alpha from reoccupying alpha-frozen orbitals and beta from
reoccupying beta-frozen orbitals without incorrectly excluding the other
spin's spatial directions. Zero-active-electron windows are valid, including
one or both empty active spin sectors.

If all ordinary retained directions would disappear, the structural-demotion
rule keeps one certified candidate active or selects an allowed complement
direction. It never constructs a zero-rank block or indexes an empty coupling
range.

## Ledger, Fields, And Recovery

Each cut stores a linked ledger generation, spin counts, the accumulated
direct field, spin exchange fields in the retained block basis, and a frozen
energy scalar. A physical frozen column is stored once per generation and
recovered by walking the ledger. Cross energy between left and right ledgers
is added once when a window is formed.

Block growth updates fields and scalars recursively from the final stored block
basis. Detailed diagnostics can compare that recurrence with direct physical
contractions, but the ordinary path does not rebuild full physical fields on
every window.

At each complete sweep, the solver builds full physical Focks and performs two
audits:

1. Frozen residual coupling checks
   `(I-D_sigma) F_sigma C_f,sigma` at a fixed roundoff scale.
2. On prospective convergence, an Aufbau ordering check compares occupied and
   complement spectra.

The Aufbau audit is dense `O(N^3)` in this first implementation and is one
reason the mode is limited to small density-density validation rather than
long-chain production optimization.

Audit triggers are collected for both spins, thawed once, and followed by one
block-family reconstruction and another forced sweep. Observers are called
only after complete sweeps; incomplete stale-recovery attempts do not produce
observer records.

Containment is certified before ledger columns are trusted. A stale ledger
cannot silently replace a determinant direction. Structural demotion and stale
recovery preserve the determinant projector, and ledger generations are
released when no cut or in-flight successor references them.

## Accepted Evidence

The original six-function algebra oracle reported maximum errors of
`6.66e-16 Ha` in energy, `2.29e-16` in Fock Frobenius norm, `2.22e-16 Ha` in
the recursive scalar, and `1.36e-16` in recursive fields. A forced `2.5e-5`
coupling was recovered within `1.76e-18`.

The accepted F1 validation included exact default/off parity, RHF and UHF,
common- and split-one-body routes, spin-specific exclusion, zero-active-spin
cases, damped/pure separation, stale-ledger rebuilding, structural demotion,
recursive Fock/scalar parity, rotation-invariant residual audits, exact ledger
release, observer ordering, and cached/projection rejection controls. The
independent acceptance run passed 519/519 assertions. A public `N=12` UHF
demotion reproducer returned `39.0 Ha` with zero Gram error, and the split-H
demotion matrix and frozen-subspace rotation oracle passed.

Current package tests retain these numerical contracts while avoiding a public
freeze threshold or backend extension.

## Limits And Related Experimental Work

Roundoff-exact freezing does not cover finite leakage, cached-sliced frozen
fields, cutoff schedules, successor rebasing, or retirement ledgers. Those
were explored on separate experimental branches and are not part of the
production contract summarized here. The mode supplies no memory or timing
claim beyond correctness.

A scalable matrix-free Aufbau audit, threshold semantics, sliced support, and
any simple frozen-local eigensolve would require separate designs. Do not infer
those features from the `:roundoff_exact` keyword.

-- hfdmrg-manager@rh310l.ps.uci.edu
