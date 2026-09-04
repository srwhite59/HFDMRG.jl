# Packet 10A5R2A same-machine performance reconciliation

Packet 10A5R2A confirms the compact-UHF performance stop with a fresh frozen
PPP denominator and localizes it to the remaining recipient entry-build and
publication representation. Packet 10A5R2 evidence is preserved unchanged;
the recipient H1000 nonlinear solve was not rerun and production source was
not modified.

## Contemporaneous H1000 decision

Frozen PPP commit `0bcd4f11831a98ae4bf4fb3d09eac5968df9d96c` was cloned
read-only into machine-local `/tmp` storage and run under Julia 1.12.6 with
four BLAS threads, direct block-native atomic-Neel initialization,
`state_maxdim=16`, `state_cutoff=1e-12`, and `energy_tolerance=1e-11`.
The internally timed nonlinear boundary is the accepted
`_independent_spin_stored_nonlinear_sweep!` call and includes terminal
coefficient reconstruction while excluding fixture construction,
initialization, replay, producer audit, and record writing.

Fresh PPP converged in 10 half-sweeps to
`-2282.257073261816 Ha` in `95.0519 s`. The recorded whole-process wall time,
including small compilation warmup and construction, was `128.9791 s`.
The machine load averages moved from `3.34/3.53/3.12` to
`5.72/4.22/3.43`; peak RSS was 1,582,235,648 bytes.

The unchanged recipient result was `289.4157 s` for 12 half-sweeps and
represented energy `-2282.2570732618165 Ha`. The contemporaneous ratio is
therefore `3.04482`, above the mandatory `1.5` stop. The earlier historical
PPP denominator of `123.7554 s` was conservative rather than unfair. The
H100 PPP timing used at R1B was not a useful warmed scaling denominator: its
30.18-second nonlinear time is inconsistent with the fresh warm H1000
per-center rate and was dominated by one-time or cold-run overhead.

## Schedule and proposal normalization

Recipient convergence visits 11,989 centers: 1,000 on the first half-sweep
and 999 on each of eleven later half-sweeps. PPP visits 10,000 centers in ten
terminal-complete half-sweeps. The 1.1989 schedule factor explains only part
of the total ratio. Recipient time is `24.1401 ms/visit`; PPP time is
`9.50519 ms/visit`, leaving a `2.53968` normalized per-center regression.

Every recipient center uses one Fock, one eigensolution per spin, one
post-factorization true-energy evaluation, and one publication. PPP uses the
same Fock/eigensolution/publication counts but three true-energy certificates
per center. The recipient remains slower despite doing fewer energy
certificates, so proposal count and convergence semantics are not the cause.

## Single-half-sweep localization

One qualification-only first-half H1000 probe was run per engine. The
recipient took `21.7758 s` with sampling and allocated 6,035,817,984 bytes;
PPP took `6.41474 s` and allocated 1,605,834,480 bytes. This is a `3.39` time
ratio and `3.76` allocation ratio for exactly 1,000 centers in each engine.

Recipient phase timing is decisive:

- factorization, entry build, and moving-root work: `16.5849 s` (76.16%);
- incoming plus candidate center preparation: `4.19347 s` (19.26%);
- incoming Fock and energy: `0.27366 s`;
- proposal copies and two eigensolutions: `0.30845 s`;
- candidate true energy: `0.27736 s`;
- publication, root copy, preflight, and terminal invariant combined:
  approximately `0.1374 s`.

PPP's tagged half-sweep reports `1.02221 s` for 4,000 preparations,
`0.69498 s` for 11,988 state-block factorizations, `1.05148 s` for 3,996
periodic core complements, `0.22781 s` for 2,000 local eigensolutions, and
`0.11304 s` for 4,000 UHF Fock calls. Its full stored local-update scope is
`5.95543 s`.

Sampling inside the recipient's dominant phase identifies
`_uhf_physical_block`, `seed_uhf_interval_block`, `_seed_cross_hartree`,
`build_uhf_outer_block`, `_fill_merged!`, packed-pair transforms, generic
matrix products, and cross-Hartree SVD compression. It does not identify the
accepted direct two-sided center Fock as the bottleneck.

## Root cause and correction seam

R1B successfully replaced the center contraction, but its adapted entry
lifecycle retained `UHFOuterBlock` as two `RHFOuterBlock` owners. Every
accepted move therefore creates a physical block, reconstructs merged dense
per-spin packed-pair fields, transforms and allocates new durable pair fields,
and recompresses transformed cross-Hartree factors.

Frozen PPP instead routes `_prepare_independent_move!` through
`_build_independent_spin_entry!` and
`_independent_interaction_projection!`. Its reusable
`_IndependentSpinEntryBuildWorkspace` maintains separated exchange
`internal`/`connector` records, cross-Hartree capacities, completed-core
capacity, and exact-sized durable publication. The R1B provenance description
claimed this entry-build family at an algorithmic level, but the recipient
adaptation retained the slower RHF-block publication arithmetic rather than
the donor's separated-record projection.

This is not chi growth: recipient maximum spin rank/pair rank remains 5/15
from H100 to H1000, recipient workspace grows only from 50.91 MB to 54.44 MB,
and PPP workspace from 41.84 MB to 43.63 MB. It is also not compilation,
system load, proposal count, or a timing-boundary mismatch; the paired
single-half probes were compiled first and executed over identical visit
counts.

The smallest correction seam is the entry owner used by
`_uhf_build_candidate!`: replace `_uhf_physical_block` plus
`build_uhf_outer_block` with a contracted port of PPP's reusable
`_IndependentSpinEntryBuildWorkspace`, separated exchange/cross-Hartree
projection, and exact-sized durable publication. Retain the accepted direct
two-sided center, public facade, recipient schedule, monotone transaction, and
historical backends. No correction is implemented in R2A.

## Scope and readiness

R2A ran one fresh full PPP H1000 solve and one single-half H1000 probe per
engine. It ran no recipient full trajectory, RHF case, producer audit,
H100/H200/H500 scaling point, package suite, GaussletBases case, or larger
benchmark. The compact public transfer remains **not ready** pending a
separate authorized entry-build correction packet and bounded requalification.
