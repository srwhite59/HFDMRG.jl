# Hosted CI H20 Finite-Test Contract

Date: 2026-09-09
Candidate: `844ad927394e935cd2b04c99035ad51fe93441de`
Hosted run: `34409081278`, test job `102659008293`
Status: test-only correction prepared; unstaged and uncommitted

## Finding

The sole hosted failure was an independently started physical-versus-reflected
finite RHF energy comparison: `2.0797360278379529e-7` against `2e-7`. The bound
entered with `18afb6f` and was explained as the sum of two `1e-7` convergence
envelopes. That explanation is not valid. The compact stopping rule compares
energies two half-sweeps apart within each trajectory and requires two quiet
comparisons. It does not bound either endpoint's absolute error from a common
stationary energy.

No hosted coefficient matrices, projectors, complete energy histories, or
rank histories were retained. The hosted endpoints therefore cannot be
classified as the same or different HF basin. The earlier local same-basin
inference and proposed `2sqrt(energy_tolerance)` projector bound are withdrawn:
neither follows without curvature, gap, and actual-error information.

## Replay semantics

After nonlinear stopping, the RHF facade records `replay_before`, performs two
`run_fixed_half_sweep!` calls, and records `replay_after`.
`replay_energy_change` is their difference. Fixed replay performs state-link
factorization, optional finite recompression, block publication, root
close/open, and turnaround. It does not construct or accept a new nonlinear
Fock eigensolution. The value therefore certifies represented-energy stability
under a full transport/recompression replay, not nonlinear stationarity.

For a finite converged result, `abs(replay_energy_change)` gates acceptance at
the requested energy tolerance plus the existing Float64 energy envelope. An
excess changes the result to `reason=:finite_replay_floor`. The exact path uses
its roundoff envelope and `reason=:exact_replay_failed`.

`replay_projector_change` compares root covariance entries before and after
replay only when root ranks match. Those entries can be expressed in different
outer bases. It is a reported root-coordinate diagnostic and does not gate
RHF or UHF result acceptance. This is consistent with the established policy
that compact convergence has no mandatory projector or residual condition;
no production safeguard defect was found.

## Correct finite H20 contract

Each orientation is now checked independently for:

- `converged` with `reason=:energy_converged`;
- requested/effective tolerance values and their reported reason;
- two terminal same-orientation changes satisfying the explicit `1e-7`
  request for this bounded fixture, without applying the final reported
  threshold retrospectively;
- the actual finite representation-replay energy bound, reconstructed from
  `replay_after`, `replay_energy_change`, and the root rank;
- zero later original-row reads, one 19-entry collection, and ten occupied
  orbitals;
- independently reconstructed dense represented energy within `2e-10`;
- orbital orthonormality and physical-density idempotency within `2e-11`.

The generic finite cross-orientation endpoint comparison is removed rather
than enlarged or replaced by an inferred projector bound. The existing exact
orientation tests and the independent monotone-publication transaction test
remain unchanged.

## Local focused evidence

Julia 1.12.6, macOS, ILP64 OpenBLAS, four BLAS threads:

- ordinary physical/reflected energies:
  `-4.4883418679544125` / `-4.48834175434509`;
- ordinary half-sweeps: 40 / 39; cross difference `1.1360932283821512e-7`;
- checked-bounds energies:
  `-4.488341851731549` / `-4.48834189980209`;
- checked-bounds half-sweeps: 40 / 45; cross difference
  `4.807054043709513e-8`;
- independent dense certificate error was zero at printed precision;
- representation-replay energy changes were at most `1.78e-15`.
- focused ordinary contract: 8/8 assertions passed;
- focused checked-bounds contract: 8/8 assertions passed.

The differing ordinary and checked-bounds trajectories reinforce that the
discarded endpoint assertion was path-sensitive. They do not establish a
basin classification.

## Static inventory of analogous comparisons

Independent nonlinear endpoints, not direct symmetry transforms:

- finite H2/H10 RHF physical-versus-reflected energy, `2e-7`;
- exact H2/H10 RHF physical-versus-reflected energy, `2e-11`;
- exact H20 RHF physical-versus-reflected energy, `1e-7`;
- finite/exact H2/H10 UHF physical-versus-reflected energy, `2e-7`;
- finite/exact H20 UHF physical-versus-reflected energy, `2e-7`.

These remain unchanged pending case-specific evidence. The UHF spin-swap
checks compare the same Hamiltonian under the actual alpha/beta symmetry and
remain true covariance tests. Fixed-state reflection/reopening tests likewise
transport one supplied state and remain distinct from comparing independently
optimized endpoints.

No production, API, solver-control, CI-version, external-repository, or
visibility change is part of this correction.

-- hfdmrg-manager@Mac-Studio
