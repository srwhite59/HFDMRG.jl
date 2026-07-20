# Paper-Manager Review: Cached Sliced Milestone 5b

Date: 2026-07-20
Reviewed through commit: `d6b1f2a`
Decision: **M5b characterization accepted; cache program complete; current solver stationarity remains an algorithm issue**

## Bottom Line

The ordinary solver reaches its sweep-energy criterion at sweep 60. A forced
140-sweep trace remains on the same broken-symmetry occupied branch and
reaches the oracle energy to roundoff, but its full occupied-virtual Fock
residual plateaus near `6.33e-7` rather than reaching the oracle residual of
`1.32e-8`.

This does not reopen the cached-backend result. The first 50 rows reproduce
byte-for-byte, every sweep uses the intended 98 local cached windows with zero
fallback, and all independent energy and determinant-integrity checks pass.
The M1--M5 cache/performance program is complete.

## Scientific Interpretation

The current production stopping rule is an energy-change rule, not a
stationarity rule. At its ordinary sweep-60 stop, the energy change is below
`1e-11 Ha`, but the full Fock residual is still `1.18e-5`. Continuing removes
most of that residual before reaching a floor around `6.33e-7`.

The floor is consistent with, but does not by itself prove, limitation by the
current local `scf_cutoff=1e-12` and damped mixed-density microiteration. The
local loop can stop after a small energy change even though the occupied
orbitals are not eigenvectors of the final rebuilt Fock matrix. A tighter
local tolerance, a determinant-variational update, or an external polishing
step could therefore change the floor. That is a separate algorithm study.

The late Z-spin miss is not evidence that the trajectory leaves the branch.
At sweep 140 the maximum projector distance is only `6.88e-7`, the minimum
singular value is `0.999999999999881`, and the energy is within `1.62e-13 Ha`
of the oracle. The observed `2.81e-6` Z-spin difference is also well below the
rough `3.15e-5` operator-norm bound obtained from `||Z||=22.87` and the
projector difference. The historical `1e-6` Z-spin gate remains recorded as
failed, but it should not become a general stationarity criterion.

## Accepted M5b Facts

- Ordinary energy convergence occurs at sweep 60.
- Energy, minimum-overlap, energy-change, projector, and S2 gates remain
  satisfied after their first crossings at sweeps 36, 56, 60, 62, and 76.
- The Z-spin gate is satisfied only on sweeps 83--88.
- The minimum residual is `6.3260414e-7` at sweep 134; sweep 140 gives
  `6.3275583e-7`.
- The largest returned-energy mismatch is `2.20e-13 Ha`; the largest
  orthonormality/idempotency error is `4.56e-15`.
- Energy fluctuations after sweep 80 are at roundoff scale.

## Next Boundary

Do not add another cache milestone. Before promotion to `main`, perform one
docs-only compaction pass: preserve a compact final design/evidence account,
the reusable Be runner, and the reduced 140-sweep data, while removing
intermediate pass reviews and append-only status material made redundant by
the final account. Git history remains the detailed process archive. Do not
change production source or tests during compaction.

After promotion, convergence work should occur on a separate branch from the
new `main`. Do not merge the old March worktree wholesale because it predates
the sliced correction and cache program. Reconcile and port only the reviewed
determinant-variational/history pieces needed for matched tests. The first
cheap diagnostic should compare the sweep-140 state under a tighter local SCF
tolerance and one full-space Fock/polishing step before choosing a production
algorithm change.

-- hfdmrg-paper-manager@faraday30.ps.uci.edu
