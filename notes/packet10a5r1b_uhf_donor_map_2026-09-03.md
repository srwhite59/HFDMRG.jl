# Packet 10A5R1B compact UHF donor map

Frozen donor: PPP `0bcd4f11831a98ae4bf4fb3d09eac5968df9d96c`.
Recipient base: HFDMRG `83cc6f112b0f477736f87e702d98402727ff3cfc`.
This map was recorded before any R1B production-source edit.

## Live donor owners

- `_IndependentSpinTwoSidedCenters.jl` (745 lines, SHA-256
  `5fe6aae1595ee0b6b780a5101e9a31603007494e62b5d1fa26a632f307cb5a99`):
  `_IndependentSpinSideView`, `_PreparedIndependentOneSpinTwoSided`,
  `_IndependentSpinTwoSidedWorkspace`,
  `_PreparedIndependentSpinTwoSided`, `_prepare_one_spin_two_sided!`,
  `_self_spin_fock!`, `_internal_cross_field!`,
  `_independent_cross_fields!`, `_prepare_independent_spin_two_sided!`,
  `_independent_spin_two_sided_fock_core!`, the checked Fock wrapper, and
  the matrix action. This is the authoritative center-arithmetic slice.
- `_IndependentSpinCrossHartree.jl` (505 lines, SHA-256
  `99a1949824e92f1df17ba43866f31731a5ef6ad3a34184c300b50e7ab9dfe3d2`)
  and `_IndependentSpinInteractionBlocks.jl` (310 lines, SHA-256
  `1aaf458e547535ba1d1fa32ebd5cb9327d77b52dc6d1d6b4dd2122eb4aeb4716`):
  compressed block-local opposite-spin Hartree factors, independent exchange
  records, shared Hartree connectors, reflection, and storage ownership.
- `_IndependentSpinEntryBuildCapacity.jl` (103 lines, SHA-256
  `a9c4cc3527d8fef15b1f669ce881be820407bfeea4ad69795a22681254fd5dc8`)
  and `_IndependentSpinCompletedCore.jl` (287 lines, SHA-256
  `73307166e1475e5dcfc145a54acfbb1d7e1caa760ee854970008e5568d49e062`):
  one reusable entry-build capacity and separated completed-core ownership.
- `_IndependentSpinLifecycleFoundation.jl` (542 lines, SHA-256
  `8620066248c0b736fb2ce4c0f11fa003d1a67c019088323a3310f9371cfd452f`),
  `_SpinResolvedStateLinks.jl` (721 lines, SHA-256
  `454129caa36834981c4c90b45c1951920110a3f8ffe6c36560c483ec7c6801b4`),
  and `_LocalProductInitialization.jl` (507 lines, SHA-256
  `be4ddc0e0ebc576a484bb2bdb7c63166cb24b5229feba87cf3e33b4a739579c1`):
  independent moving roots, close/open operations, reversible spin links,
  block-local product initialization, reflection, and spin swap.
- `_IndependentSpinCollectionLifecycle.jl` (1,486 lines, SHA-256
  `6f2aaebf06ed07a705f6b0706fc1d29fd97b630087df5921d4ba0a5f465f3bf6`):
  `_build_independent_spin_entry!`,
  `_initialize_independent_spin_composite`,
  `_prepare_independent_move!`, `_replace_independent_entry!`,
  `_IndependentSpinLifecycleCenterCapacity`,
  `_prepare_independent_center!`, terminal reflection/reconstruction, and
  `_independent_spin_fixed_state_half_sweep!`. Only their invariants and the
  function-level arithmetic already represented by recipient owners are
  required; the whole PPP collection/schedule module is not transferred.
- `_IndependentSpinMonotoneTransaction.jl` (396 lines, SHA-256
  `e4a65c79d9ff1d850e16c50a6a2988d42c2aecfc738dc3a6e36a48e2338cdb6c`):
  `_independent_transaction_energy!`,
  `_factor_independent_transaction_candidate!`,
  `_publish_independent_transaction!`, common chord/geodesic control, and
  `_independent_spin_monotone_transaction!`.
- `_IndependentSpinStoredNonlinear.jl` (292 lines, SHA-256
  `b9f3eb17f9e43dc5355666a4fae22806609ee422dc67a04c2065f9e8a59e4835`):
  reusable nonlinear workspace, per-half traversal, terminal turn, evidence,
  and convergence accounting.

All transferred recipient code will be marked adapted unless an individual
slice remains byte-identical. Module names, owner types, field access, the
two-group schedule, public facade, and result records necessarily differ.

## Recipient mapping

- `UHFOuterBlock` and `UHFCrossHartree` in `src/uhf_blocks.jl` already own the
  donor-equivalent separated alpha/beta block spaces, shared Hartree channel,
  independent exchange records, and compressed block-local cross-Hartree
  factors. `UHFBlockBuildWorkspace` is the recipient reusable entry capacity.
- `UHFMovingRoot`, `UHFFixedLifecycle`, and `UHFLifecycleWorkspace` in
  `src/uhf_lifecycle.jl` already own one exact-sized collection, one moving
  independent-spin root, close/open motion, block-native atomic-Neel
  initialization, terminal turnaround, reflection, spin swap, and terminal
  reconstruction. They replace the corresponding PPP collection wrappers;
  no second collection or reflected operator copy is needed.
- `_monotone_uhf_center_transaction!` and `_uhf_solver_loop!` in
  `src/nonlinear_uhf.jl` map to the donor monotone transaction and stored
  nonlinear sweep. R1B retains the accepted one-Fock/common-geodesic contract
  and makes the unchanged baseline lazy on the accepted-full path, as in the
  R1A reconciliation and donor transaction ordering.
- `FastRHFPreparedCenter`, direct periodic-center panels and exchange maps,
  block feature records, channel transfer, pair packing, state-link
  factorization, and the accepted recipient traversal schedule are reusable
  R1A/shared machinery. UHF adds only independent-spin ownership and the
  donor direct opposite-spin Hartree contractions.

## Obsolete recipient structures

`src/uhf_centers.jl` currently owns `PreparedUHFCenter`, cap-sized
`cross_alpha`/`cross_beta` packed-pair panels, `_append_internal_cross!`,
`_append_physical_cross!`, `_append_component_cross!`, and
`_prepare_center_cross!`. The nonlinear transaction prepares this expanded
representation for incoming, unchanged baseline, and full proposal. The
entire expanded-panel center is replaced; none of those structures or helpers
may remain reachable.

Once UHF no longer embeds `PreparedRHFCenter`, the superseded dense compact
center in `src/interaction_centers.jl` has no production owner. Its tests will
be moved to the fast center contract and the file removed from the compact
include graph. Block-local `RHFOuterBlock.pair_field` records remain valid
independent self-interaction records; no total-center pair field remains.

## Bounded transfer and ownership plan

1. Add one maximum-capacity fast UHF center owner containing two independent
   one-spin prepared records, square spin Fock/field workspaces, packed
   left/right density vectors, compressed-factor vectors, and one shared MPO
   channel vector. Port the donor prepare, self-spin, internal-factor, channel,
   physical-center, Fock, true-energy, and action arithmetic.
2. Wire `UHFLifecycleWorkspace` to that owner and prove H2/H10 exact arithmetic
   and both-orientation fixed lifecycle behavior before changing nonlinear
   acceptance.
3. Adapt the donor accepted-full transaction ordering: incoming Fock and two
   eigensolutions, one post-factorization true-energy trial, one publication;
   build baseline and common-parameter geodesic candidates only on fallback.
4. Delete the expanded UHF center and, after its last owner disappears, the
   obsolete dense compact center. Keep one compact UHF route and one compact
   RHF route.
5. Gate H20 numerical/reflection/spin-swap/replay/ownership/allocation before
   one matched H100 run. Package-wide suites follow only if H100 is at most
   1.25 times PPP.

At `state_maxdim=16`, each spin center dimension is at most 68. The replacement
owns square `O(68^2)` spin fields and block-side packed vectors, not
`K=2346` by thousands-of-columns cross panels. Center memory should remain a
few megabytes and independent of chain length; collection memory remains
linear in the number of published blocks. Warm nonfactor prepare/Fock/action
must allocate zero bytes.

The current compact slice is 5,599 physical production lines. Removing the
339-line expanded UHF center and 575-line obsolete dense compact center before
adding the consolidated donor arithmetic projects a final 5,400--6,100 lines,
below the 6,800-line gate. No source duplication or third backend is planned.
