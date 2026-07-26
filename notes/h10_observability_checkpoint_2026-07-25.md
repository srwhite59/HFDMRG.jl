# H10 observability and cache-coordinate acceptance

Date: 2026-07-25

Status: accepted bounded observability and first-sweep H10 window preflight
after the post-reorthogonalization cache-coordinate repair. This is an
implementation and preflight acceptance record, not a completed H10 solver
result.

The accepted implementation adds a private run-local diagnostics collector for
environment spectra, local iteration and damping history, and eigensolver
timing. The scratch-oriented H10 preflight runner adds a meter for cache
routing, retained storage, and backend phase timing. The public solver API and
`SweepInfo` are unchanged.

## Accepted lineage

- Cache-coordinate repair:
  `d6abbec15391204aec5be32b20368802018af00b`.
- Repair documentation head:
  `a39652d641c27d70a5d39b0c9513dad68103781b`.
- Original branch-only observability checkpoint:
  `34883f5f93f40c1ffb8f7fbbe0a60b2521fba745`.
- Observability replayed on the repaired parent:
  `ff45284ae3caba7be24dff8b10cccd2c4fd12bd5`.
- Accepted preflight runner SHA-256:
  `8ec8e7a836f828a2dec19921cd3d807398907591a0b0fd97038e82c94dd19655`.

The repair preserves the existing final-basis calculation while applying its
post-SVD reorthogonalization to both the old-block and center recurrence maps.
The resulting maps now describe the stored `phi_new` for common-H and split-H
one-body caches and for every interaction-backend absorption. The realistic
H10 preflight therefore corrected a coordinate inconsistency affecting the
one-body, density-density, and target-residual recurrences.

## Diagnostics-off default-route parity

The deterministic eight-route harness has SHA-256
`a7ea53707df19efa84394c8ca026c9ee39394d8033d330a2140d9ab28408d1a3`.
It covers density-density RHF, common-H UHF, and split-H UHF; cached-sliced
RHF, common-H UHF, and split-H UHF; target-residual UHF; and
projection-sliced UHF. No route passes an observer or `_diagnostics`.

Fresh repair-only results from `a39652d` were serialized with SHA-256
`9b640fa86f7cd98b796b3207d95c92ab7f2cf88f5bf7591ffa0b2c57af99040e`.
Running the same harness at `ff45284` passed `isequal` against that artifact,
proving exact returned alpha orbitals, beta orbitals, and energy on every
route. The matching per-route result hashes are:

| Route | Result SHA-256 |
|---|---|
| density-density RHF | `3dfc7e10544564caf7d50c933c9655ed90753e7b02263c11e853af3bcdfb20f7` |
| density-density common-H UHF | `808fa73a7995315b20ec82f6094d06d2dbccb6e5882119caa5b0b88129ed7ec2` |
| density-density split-H UHF | `25ae5deefea09d159384c032b6cd58fcef9d6b5fa5362ee9eff92a039df00ad1` |
| cached-sliced RHF | `0f6f96b5683032d309cd7c5c95c3d2cf65ead59318ed1d278a5797b0a7445497` |
| cached-sliced common-H UHF | `eb7765fec0cf6d6a4b2a88e96ff5f7f12ff8a4f52ab4947e4399c0d9ed83ebde` |
| cached-sliced split-H UHF | `58ca4e985538cf934f5db71aa0c910ef974a09a199b9935b7396501e2d914d68` |
| target-residual UHF | `b10fc8d3b69ec30b07719f8daf70fb03c9d97bbc3aa6321a1c7be33593bc1c3d` |
| projection-sliced UHF | `5ee2e436e803eadc1290da30f487625416d6769913b01fba2cb925b02e1edc18` |

## H10 selected-window oracle

The natural and odd-even runs used the unchanged one-body maximum, per-spin
Fock Frobenius, and window-energy gates of `2e-11`. All selected captures at
global windows 1, 200, 399, and 796 passed:

| Order | Captures | One-body maximum | Alpha Fock Frobenius | Beta Fock Frobenius | Energy maximum |
|---|---:|---:|---:|---:|---:|
| natural | 16 | `2.4868995751603507e-14` | `2.071465337656388e-12` | `2.0695867964933997e-12` | `5.364597654988756e-13 Ha` |
| odd-even | 14 | `3.694822225952521e-13` | `2.244820382612786e-12` | `2.244855284901849e-12` | `1.6981971384666394e-12 Ha` |

Supplementary maximum basis-Gram errors were
`1.3323158351614485e-13` natural and `3.071984385641888e-14` odd-even.
Maximum Fock-antisymmetry ratios were `3.448805104399223e-18` and
`1.3602658727951171e-18`, respectively, against
`1.1368683772161603e-13`. Odd-even ordering covariance errors were
`1.8474111129762605e-13 Ha` in energy and
`5.1437201588831186e-14` in Fock norm; natural-order covariance was exactly
zero.

Both orders recorded exactly two cache initializations, 398 initialization
absorptions, 796 sweep absorptions, and 796 cache-local sweep windows.
Fallback and projection-window counts were zero in initialization and sweep.

## Package and backend validation

- The complete integration package suite passed all displayed testsets,
  277/277.
- The repeated post-reorthogonalization recurrence test passed 88/88.
- The hidden run-diagnostics test passed 30/30.
- The maintained fixed/ragged cached-versus-projection verifier passed.
  Its maximum window Fock discrepancy was
  `5.10702591327572e-15`, and its maximum sweep-energy discrepancy was
  `1.9539925233402755e-14 Ha`.
- `git diff --check` passed.

## Work, allocation, and retained storage

These Julia 1.12.6 measurements used one Julia thread and one BLAS thread on
`rh310l.ps.uci.edu`.

| Order | Total wall | Solver excluding observer | Observer | Allocated bytes | Retained cache-array bytes | Cache `summarysize` | Observer-boundary RSS |
|---|---:|---:|---:|---:|---:|---:|---:|
| natural | `6.598680667 s` | `5.351367 s` | `1.247313667 s` | `9404278416` | `986637584` | `1015057296` | `7308492800` |
| odd-even | `7.339424667 s` | `6.190468583 s` | `1.148956084 s` | `10675441920` | `1020306016` | `1048725728` | `8512241664` |

The runner's final whole-process maximum-RSS messages were `7311048704` bytes
natural and `8515026944` bytes odd-even, after the remaining output work.

Phase counters and Julia-level timings were:

| Order | Initialization total | Cache init | Initialization absorption | Sweep absorption | Window construction | Fock contraction | Spin eigensolves |
|---|---:|---:|---:|---:|---:|---:|---:|
| natural | `3.470904375 s` | 2 / `0.000749208 s` | 398 / `3.249601664 s` | 796 / `1.249429691 s` | 796 / `0.367349162 s` | 3538 / `0.056557405 s` | 5484 / `0.094856162 s` |
| odd-even | `3.984075959 s` | 2 / `0.000819208 s` | 398 / `3.740189335 s` | 796 / `1.315576282 s` | 796 / `0.507435967 s` | 3098 / `0.064721070 s` | 4604 / `0.081232076 s` |

The independent selected-window oracle work was
`2.120323125 s` and `2033623472` allocated bytes for natural order, and
`1.85678525 s` and `1790165680` bytes for odd-even order. It is separate from
the solver/observer split above.

## Rank and local-iteration summary

The complete fixed left and right physical edges retained rank 4.
For natural order, initialized right-environment ranks ranged from 1 to 10
(258 of 399 were rank 10), left-to-right grown-left ranks ranged from 1 to 5
(351 of 398 were rank 5), and right-to-left grown-right ranks ranged from 2
to 5 (338 of 397 were rank 5). For odd-even order, the corresponding ranges
were 1 to 10 (327 of 399 rank 10), 1 to 8 (275 of 398 rank 5 and 62 rank 8),
and 1 to 5 (367 of 397 rank 5).

| Order | Microiterations | 2 / 3 / 4 iterations | Locally converged windows | Cap-exhausted windows | Windows with rises | Rise updates |
|---|---:|---:|---:|---:|---:|---:|
| natural | 2742 | 184 / 74 / 538 | 340 | 456 | 43 | 43 |
| odd-even | 2302 | 296 / 290 / 210 | 597 | 199 | 234 | 234 |

Natural-order damping transitions were `1.0 -> 1.0` in 728 windows,
`1.0 -> 0.5` in 43, and `0.5 -> 0.5` in 25. Odd-even transitions were
`1.0 -> 1.0` in 457 windows, `1.0 -> 0.5` in 225,
`0.5 -> 0.5` in 105, and `0.5 -> 0.25` in 9.

## Scope

These measurements are forced first-sweep observability and selected-window
preflights. Their returned orbitals and energies are not H10 endpoints. They
do not constitute a conventional-UHF comparison, a complete scientific
HFDMRG run, a stationarity or stability result, or a physical-energy
benchmark. Full matched-solver execution remains subject to a separately
frozen execution manifest.
