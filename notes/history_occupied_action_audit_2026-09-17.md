# Historical history-audit occupied-action correction

Baseline: `b6af4140e8aece644bd045613fc44b738438f588` (`v0.1.0`).
This change is unstaged and uncommitted.

## Consumer contract

The full-basis RHF and UHF physical-audit routines are called after bootstrap,
for each history proposal, and after each cleanup segment. Their energy,
orthogonality, and spin values participate in proposal acceptance; their
occupied--virtual residuals drive the outer history stopping test and are
published in events and results. They do not build the history basis, reduced
Hamiltonian, reduced Fock matrices, Pulay error matrices, preconditioner, or
candidate orbitals. Those reduced-model operations remain unchanged.

The history entry points do not accept the separate `frozen_occupied` control.
Consequently there is no frozen history-audit variant to alter. The existing
frozen-occupied implementation and its 69-test regression group remain
unchanged and passed in the complete suite.

## Exact algebra

For density-density RHF, the physical audit now forms only the occupied action

`F C = H C + 2 diag(V d) C - (V .* (C C')) C`, with
`d_i = sum_a C[i,a]^2`.

The exchange action is evaluated in 64-row blocks. Each block forms
`C_block*C'`, multiplies it elementwise by the corresponding rows of `V`, and
contracts directly with the requested occupied columns. No full `N x N`
density or Fock temporary is formed. The energy is the identical contraction
`dot(C, F*C + H*C)`. UHF shares the Hartree vector, evaluates independent
alpha and beta exchange actions, and retains the original one-half spin
energy factors, residual definitions, `K` normalization, and `S^2` formula.

The sliced history route uses the shared `_sliced_hartree_data` preparation and
`_sliced_occupied_action!` recurrence directly. It does not construct
`_SlicedFockAuditContext`, because history energy and residual do not consume
that context's row bounds, diagonal/neighbor blocks, or preconditioner data.
The target-residual/frozen audit still constructs the full context and retains
its setup reasons, scale, neighbor blocks, counters, and preparation semantics.

For a multi-slice layout, history scratch consists of packed per-slice Hartree
blocks, per-slice occupied densities, and one exchange block for the current
slice pair (with its same-sized density block). Cumulative exchange allocation
is still quadratic in total length,
and if the layout has one full-size slice then an exchange block is itself
`N x N`. The density-density path owns at most a `64 x N` exchange panel (whose
dimensions coincide with `N x N` only when `N <= 64`) rather than separate
full density and Fock matrices. The corresponding sliced claim applies only
when the maximum slice is smaller than the full space; it is not a blanket
statement for an arbitrary sliced layout.
Reduced history-space density and Fock matrices remain intentionally present:
they are small and are required for proposal generation.

The dense formulas remain in `test/history_occupied_audit.jl` as independent
small-system oracles, not as a reachable solver fallback.

## Numerical and decision evidence

The focused tests cover stationary and nonstationary RHF, unequal-spin UHF,
an empty beta sector, occupied-gauge rotations, density and sliced backends,
occupied Fock actions, energy, residual normalizations, Gram errors, and spin.
They also construct real reduced-history RHF and UHF proposals and verify that
the dense and occupied-action audits make identical acceptance decisions.
The complete pre-extraction suite retained the end-to-end history
accepted/rejected trajectories and result shapes. Final focused tests verify
that the extracted helper preserves the context's RHF/UHF actions, Hartree
pack, row scale, diagonal/neighbor data, and setup status. Observed oracle
differences were at Float64 roundoff; no bounded fixture changed a decision
branch.

Validation on Julia 1.12.6 with four BLAS threads:

- focused ordinary after the final sliced extraction: 50/50;
- focused checked bounds after the final sliced extraction: 50/50;
- complete checked-bounds package regression on the exact final tree:
  2297/2297;
- whitespace: clean.

## Resources

The reproducible runner and exact records are in
`validation/history_occupied_audit/`. On the warmed `N=600`, `(n_alpha,
n_beta)=(6,5)` synthetic density fixture, RHF allocation changed from
14,533,264 to 492,832 bytes and UHF from 29,031,008 to 933,904 bytes. The
matched elapsed samples were 1.207 versus 1.286 ms for RHF and 2.797 versus
1.977 ms for UHF. These are bounded call-level measurements, not production
throughput claims. Together with the `N=100` deterministic history/oracle
tests, this small-system evidence is the resource and behavioral acceptance
basis. No repository fixture explicitly identified as q=5 was available; the
existing small fixtures exercise the same history audit and acceptance seam.

On an 80-dimensional fixture with 40 slices of size two, using six alpha and
five beta orbitals, constructing `_SlicedFockAuditContext` in history cost
803,600 bytes and 0.579 ms for RHF, and 1,604,224 bytes and 1.117 ms for UHF.
The action-only preparation reduced those measurements to 397,328 bytes and
0.369 ms for RHF, and 775,248 bytes and 0.670 ms for UHF. The dense sliced
oracles used 211,488 bytes/0.128 ms and 483,952 bytes/0.195 ms respectively.
The action-only correction therefore removes the confirmed redundant setup,
but on this tiny fixture it remains slower and allocates more cumulatively than
forming dense matrices. Its benefit is ownership bounded by slice geometry,
not the roughly 30-fold allocation reduction measured for the density backend.

Before the user clarified that larger-case RSS was unnecessary, one bounded
fresh-process synthetic preflight used `N=2400`, for which the five permanent
input objects occupied 138,451,480 bytes. The dense UHF process peaked at
904,298,496 bytes RSS and the occupied-action process at 554,188,800 bytes,
a difference of 350,109,696 bytes. Fresh-call timings (0.532 and 0.817 s) and
allocation counters include compilation and therefore are not warmed-kernel
comparisons. Peak RSS includes Julia, compilation, input construction, and the
audit; it is not a direct solver-workspace measurement. This sample is retained
as optional context only and is not an acceptance gate.

History queues and reduced-model storage are unchanged. No q9 input or run was
found in the HFDMRG checkout, and no q9 or full production calculation was
launched. The reported 9.1 GiB and 45.2/1.4 MiB preliminary figures therefore
remain external, unverified motivation. No measured or extrapolated q9 saving
is claimed. No application-level dense Fock check
was found or changed; if a calling application constructs such a check, that
allocation can still determine its end-to-end peak independently of this
library correction.
