# Packet 10A5R1A fast RHF review gate

Packet 10A5R1A passes the bounded fast-RHF replacement gate. The compact RHF
facade now reaches one separated-record center implementation. It stores
block-local interaction records, physical-center data, propagated channels,
and reusable cross-exchange maps; it never constructs a total-center
`pair_field[K,K]`. The old `PreparedRHFCenter` remains private only because
the unmigrated compact UHF implementation owns two instances. That shared
slow owner is explicitly reserved for Packet 10A5R1B removal.

The center arithmetic is adapted from the frozen PPP periodic-interaction
center algorithms, but it runs inside HFDMRG's simpler two-group lifecycle.
Whole-trajectory call-count or half-sweep identity with PPP is intentionally
not required. The accepted invariant is no redundant work per visited center,
scientific equivalence in the same finite representation, and the performance
bound.

## Lifecycle reconciliation

The H100 fixture has 1,803 sites and phase one, hence 101 traversal groups: an
initial 3-site group followed by 100 full 18-site groups. The first physical
half-sweep visits 100 adjacent centers. A terminal publication turns the root
around and the next half-sweep begins one center inward, so each later
half-sweep visits 99 centers. Nine converged half-sweeps therefore contain
`100 + 8*99 = 892` visits. This differs from PPP's terminal-complete schedule
without duplicating work inside a recipient center.

The audit did find one real duplication in the initial recipient transaction:
an unchanged baseline was factorized and measured before every full proposal.
The accepted-full path now follows the donor ordering and constructs that
baseline lazily only if fallback or rejection is needed. Each of the 892
accepted-full H100 visits performs exactly:

- two center preparations (incoming and post-factorization proposal);
- one materialized Fock construction;
- one eigensolution;
- one post-factorization true-energy evaluation; and
- one publication.

The nine additional preparations and Fock constructions are the one
convergence-invariant evaluation at the end of each half-sweep. PPP's explicit
accepted-update counters are likewise one Fock, eigensolution, true-energy
trial, and publication per visit. Its lower-level timing counts include a
separate scheduled preparation that HFDMRG folds into incoming preparation;
those implementation-layer counts are not compared as public work units.

## H100 confirmation

Julia 1.12.6 with four BLAS threads produced:

- represented energy `-158.82775662630124 Ha`;
- convergence after nine half-sweeps at the `1e-11 Ha` target;
- final same-parity changes `-1.31e-12` and `8.53e-14 Ha`;
- 892 accepted-full updates, zero fallback, and zero rejection;
- maximum state/pair ranks 9/45;
- primary solve time `5.845225791 s`, or 0.39--0.42 times the frozen PPP
  `13.913151458--15.136047125 s` records;
- primary solver allocation `213,229,552` bytes, center owner `2,684,920`
  bytes, and peak RSS `1,055,604,736` bytes;
- one 50-dimer initialization and zero later dense-row reads.

Every local publication stayed within the floating-point monotonicity
envelope; the largest raw increase was `2.02e-12 Ha`, with zero excess beyond
that envelope. Update-disabled replay changed represented energy by
`-1.14e-13 Ha` in the reflected orientation and `2.84e-14 Ha` in the physical
orientation. The recipient result differs from the independently truncated
PPP trajectory by about `1.92e-11 Ha`; this is within the governing
finite-representation comparison policy and does not imply bitwise trajectory
identity or a global HF minimum.

## Validation and boundary

H2/H10 exact arithmetic and lifecycle parity, H20 finite and exact orientation
tests, terminal reconstruction, reflection, warmed allocation, and hostile
transaction paths pass. The 53-check compatibility suite confirms the exact
eight-name export surface, zero ambiguous public dispatch, compact UHF smoke,
and historical dense, sliced, cached, target-residual, and tuple-return routes.
The complete direct package suite passes under Julia 1.12.6, including all
historical, compact UHF, producer-observer, and allocation sections. The
compact engine is 5,599 physical production lines, below the 6,800-line bound.

No H1000 or 840-dimensional GaussletBases calculation was run. Compact UHF
remains the known slow implementation for Packet 10A5R1B; no R1A evidence is
an acceptance claim for its performance. All R1A changes remain unstaged and
uncommitted.
