# Packet 10A5R2B compact-UHF entry-lifecycle donor map

This packet is based on HFDMRG commit
`3e94fd707f5980f8c522d67351f174445b942155`.  Its authoritative algorithm
donor is frozen PPP commit
`0bcd4f11831a98ae4bf4fb3d09eac5968df9d96c`.  R2 and R2A evidence remains
diagnostic and is not packet input or acceptance evidence.

## Static source-range and symbol map

Every donor slice below is **adapted**, not byte-identical.  The source-range
hash is SHA-256 over the named frozen-PPP lines including their original line
endings.  Recipient hashes are recorded at the review gate.

| Donor range and source hash | Live donor contracts | Recipient disposition |
|---|---|---|
| `_IndependentSpinInteractionBlocks.jl:26-204`, `d0a84de3bcf728a8b72010f4b9bba9e7f75ecd76acad4712bf8707b6fb2db989` | `_IndependentSpinExchangeBlock`, `_IndependentSpinInteractionBlock`, initialization and validation | Adapt into compact-UHF-only separated spin and interaction records in `src/uhf_entries.jl`.  Use `RHFStateLink`, `UnitCellInteraction`, row/cell geometry, reflection flag, and existing provenance rather than PPP module-owned link/certificate types. |
| `_IndependentSpinCrossHartree.jl:4-66`, `970d824e012b57bc241bbbec2c8ad006572b64e158355246a53ca2d7873e31fd` | exact-sized factored cross-Hartree owner and validation | Already algorithmically present as `UHFCrossHartree`; retain the accepted R1B record and tighten ownership checks only where required. |
| `_IndependentSpinCrossHartree.jl:120-345`, `b7e45b5315550c16df84cb8628fdf821325e25fca00cfe52b44413a5023a11cf` | packed explicit densities; reusable seed/grow/transform factor capacity | Adapt into the one compact-UHF entry workspace.  Preserve the donor's `old_factor_rank + 3*width` temporary-column bound and exact-sized durable factor publication. |
| `_IndependentSpinCrossHartree.jl:419-487`, `75e4a86b87a9550683b857bc81fe1942cbd3125e5139d9f19d3f3b49ec10eafc` | reflection, spin swap, and direct cross-field evaluation | Direct cross-field evaluation is already present in the accepted R1B center.  Adapt only durable reflection/spin-swap ownership; do not replace the center algorithm. |
| `_IndependentSpinEntryBuildCapacity.jl:4-101`, `9b6c016f1c923e9d072481529f2a9a216ee93b22de9bd271d3d74c382cb72796` | `_IndependentSpinParentBuildCapacity`, `_IndependentSpinEntryBuildWorkspace` | Replace `UHFBlockBuildWorkspace` with one maximum-capacity entry owner containing separate alpha/beta parent exchange, H1, completed-core, cross-Hartree, pair-map, and interaction-entry scratch. |
| `_IndependentSpinCompletedCore.jl:5-200`, `08986a9e0a5cf8df731e90dbb2885b71e0515eb24f9adea36f93ddb6cfebe751` | `_IndependentSpinCompletedCoreResult`, capacity, `_independent_spin_completed_core_capacity!` | Adapt completed-projector exchange/Hartree/core-energy arithmetic to `RHFStateLink.completed_map` and recipient packed-pair conventions.  Publish distinct alpha/beta core matrices and one shared completed Hartree channel. |
| `_IndependentSpinCollectionLifecycle.jl:285-400`, `e7de8fdb25046833e2fa6b719b7155ce98695416b2a8586be19582f14b415848` | spin H1 projection; identity parent maps; `_bare_spin_parent!` | Adapt to `BandedOneBody`, `RHFStateLink.close_map`, current exact/finite selector, and recipient periodic MPO helpers.  Keep exchange `internal` and `connector` separate throughout. |
| `_IndependentSpinCollectionLifecycle.jl:455-493`, `5605e7f14d7f7f60f661e7e6d08b6b916821e66665b2d7f3a33375380b857ef5` | `_fill_parent_inherited!`, `_transferred_hartree!` | Adapt into reusable parent-core/channel construction with recipient phase and transfer accessors; no physical `RHFOuterBlock` seed is permitted. |
| `_IndependentSpinCollectionLifecycle.jl:561-677`, `a1848c33030ddedca68756e4bf122e3969ba6b00527cacea6e12a50c81b94600` | `_independent_interaction_projection!`, `_build_independent_spin_entry!` | This is the governing replacement arithmetic.  Adapt exact publication to recipient block/root ownership and combine it with recipient H1 edge projection without changing the arithmetic or introducing a second lifecycle. |
| `_IndependentSpinCollectionLifecycle.jl:1013-1078`, `c2768f6e86420168ad7880bc03ed497eb94b300494c8e117161d0aeaa007c121` | `_prepare_independent_move!`, `_replace_independent_entry!` | The current `_uhf_build_candidate!` already owns equivalent root copy/open/select/close/atomic publication semantics.  Replace only its physical-block plus `build_uhf_outer_block` call with the new reusable entry build; retain the accepted recipient schedule and monotone transaction. |

## Recipient symbol disposition before editing

- Retain the accepted R1B center owner in `src/uhf_centers.jl`, adapting its
  side reads from `RHFOuterBlock.pair_field/left_feature/right_feature` to
  separated exchange `internal/connector` and shared Hartree channels.
- Retain `UHFMovingRoot`, collection scheduling, `_uhf_build_candidate!`
  root semantics, terminal turnaround, nonlinear transaction, compact facade,
  and producer observer.
- Replace `UHFOuterBlock`'s two `RHFOuterBlock` owners with exact-sized
  compact-UHF H1, alpha/beta exchange, cross-Hartree, shared channel, and
  completed-core records.
- Replace `UHFBlockBuildWorkspace`, `_spin_publication`, and
  `build_uhf_outer_block` with the single reusable separated-entry workspace
  and donor projection/publication route.
- Delete compact-UHF reachability of `_uhf_physical_block`,
  `seed_uhf_interval_block`, `_fill_merged!`, dense per-spin `pair_field`, and
  RHF physical-block reconstruction.  RHF definitions remain untouched and
  reachable only from compact RHF.
- Retain one `Vector{UHFEntryBlock}` collection and one live `UHFMovingRoot`;
  no shadow collection, fallback builder, or comparison engine is allowed.

The source edits are confined to replacement `src/uhf_entries.jl`, the entry-build
call seam in `src/uhf_lifecycle.jl`, the separated side reads in
`src/uhf_centers.jl`, and focused compact-UHF tests.  The final net production
change must remain below 900 lines unless a measured review is requested well
before the 1,200-line hard boundary.

## Implemented replacement and recipient provenance

All transferred families are adapted; none is described as byte-identical.
The donor range hashes above identify the frozen source.  The consolidated
recipient owner is `src/uhf_entries.jl`, SHA-256
`616d4ffc5b68c0c6b84efe72f5898e36e75a8fe54e12eddd0eb1f13705ff77de`.
It implements separate exact-sized alpha/beta H1, exchange-internal,
exchange-connector, edge, link, and completed-core records; one shared
cross-Hartree factor and one shared Hartree channel; and one reusable
`_IndependentSpinEntryBuildWorkspace`.

The donor parent/cross/completed-core arithmetic was adapted to recipient
packed-pair scaling, `RHFStateLink`, `BandedOneBody`, `UnitCellInteraction`,
the existing two-group root schedule, and exact publication constructors.
The accepted R1B two-sided center was retained and changed only to consume
the separated entry fields.  The compact RHF helper signatures were
generalized at two private exchange kernels so they can consume the new spin
record; compact RHF arithmetic and ownership did not change.  The obsolete
446-line `src/uhf_blocks.jl` owner is deleted.  No compact-UHF source references
`UHFOuterBlock`, `RHFOuterBlock`, `UHFBlockBuildWorkspace`,
`build_uhf_outer_block`, `_uhf_physical_block`, `_fill_merged!`, or
`pair_field`.

The final compact production inventory is 5,482 lines, a net increase of 391
lines over the accepted 5,091-line R1B engine and below the 6,800-line bound.
The separately retained historical dense/sliced/cached/target-residual/history
production inventory is 5,692 lines.  There is one compact UHF lifecycle, one
collection, one moving root, and one entry workspace.

## Numerical, ownership, and performance result

Focused H2/H10/H20 exact and finite tests pass in both orientations, including
dense-oracle Fock/energy parity, reflection, spin swap, terminal behavior,
failure atomicity, and zero warmed nonfactor center preparation/Fock/action
allocation.  H100 converges in 12 half-sweeps and 1,189 visits to
`-164.3389859225942 Ha` in 8.89381 seconds.  Each visit has exactly one Fock,
one eigensolution per spin, one post-factorization true-energy evaluation, and
one publication; later dense-row reads are zero.  Reusable lifecycle ownership
is 28,059,400 bytes.

The paired H1000 single-half probe takes 7.967677 seconds and allocates
224,706,464 bytes, versus fresh frozen PPP's 6.414744 seconds and
1,605,834,480 bytes.  The time ratio is 1.24209 and the allocation ratio is
0.13993, passing the 1.25 gate and authorizing one full qualification.

That single H1000 solve converges in 12 half-sweeps and 11,989 visits to
`-2282.2570732618124 Ha`, 7.28e-12 Ha from the accepted PPP value.  It takes
133.340498 seconds versus fresh PPP's 95.051920 seconds, a raw ratio of
1.40282 in the required classification band.  PPP used 10,000 visits.  After
normalizing for the recipient's accepted 1.1989 schedule ratio, the recipient
cost is 11.1219 ms/visit versus PPP's 9.50519 ms/visit, ratio 1.17009.  Thus
the remaining whole-trajectory excess is schedule/convergence semantics, not
duplicated entry arithmetic or a missing hot path.

The final state has alpha/beta ranks 5/5 and pair ranks 15/15, one
initialization, zero later dense-row reads, one collection/root, and one
Fock/true-energy/publication plus one eigensolution per spin per visited
center.  The physical/reflected fixed replay changes energy by 3.18e-12 Ha;
the physical nonlinear replay by -4.55e-12 Ha; and the spin-swapped replay by
-1.36e-12 Ha.  Reusable lifecycle ownership is 30,358,712 bytes, collection
storage is 33,411,776 bytes, nonlinear workspace summary size is 30,951,760
bytes, and peak RSS is 1,260,371,968 bytes.

## Validation and readiness

The complete direct, checked-bounds, and fresh isolated package suites pass
under Julia 1.12.6.  `solve_hfdmrg` retains 19 methods, recursive ambiguity
detection returns zero, and the public export surface remains exactly eight
names apart from the module binding.  Whitespace checks pass.  The package
Manifest remains deleted as it was before this packet; no PID marker or
repository scratch remains.

R2B is ready for public-transfer review.  The raw full-trajectory 1.25 target
is not met, but the explicitly required 1.25--1.5 reconciliation identifies
the excess as the accepted recipient schedule: the paired half-sweep and
normalized per-center gates both pass, allocation is substantially lower than
PPP, and no redundant per-center work remains.  No producer audit, producer
residual, RHF H1000, GaussletBases 840-dimensional calculation, or larger
benchmark was run, and no global-HF-minimum or producer-stationarity claim is
made.
