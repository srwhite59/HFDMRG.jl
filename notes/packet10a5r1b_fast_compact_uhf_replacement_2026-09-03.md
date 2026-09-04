# Packet 10A5R1B fast compact UHF replacement

Packet 10A5R1B reaches its unstaged review gate. The compact UHF route now
uses the frozen PPP independent-spin direct two-sided center algorithm and
reusable ownership contracts, consolidated into HFDMRG's existing typed
facade and two-group schedule. The historical dense, sliced, cached,
target-residual, and history implementations are unchanged.

The replacement owns independent alpha/beta H1, exchange, state-link, and
packed self-pair spaces plus shared MPO Hartree channels and compressed
block-local opposite-spin factors. It constructs square per-spin center
fields directly. It does not own a cap-sized expanded cross panel, dense
alpha-by-beta pair map, persistent common-union field, global coefficient
panel, second collection, or reflected operator copy. The old compact
`PreparedUHFCenter` path and the obsolete dense total-center RHF pair-field
owner have been removed rather than retained as comparison engines.

## Numerical and performance result

The frozen H20 hydrogen fixture passes practical finite and exact controls in
both orientations. Finite physical/reflected represented energies agree to
`7.11e-15` Ha, exact orientations to `1.78e-13` Ha, and the finite spin-swap
projectors align within `4.60e-10`. Every case converges, update-disabled
replay energy is at roundoff, and no later dense rows are read.

On the matched H100 atomic-Neel fixture, Julia 1.12.6 and four BLAS threads,
the recipient converges in 12 half-sweeps to `-164.33898592259422` Ha in
`22.5111` seconds. Frozen PPP converges in 10 half-sweeps to
`-164.33898592259447` Ha in `30.1823` seconds. The ratio is `0.74584`, below
the `1.25` gate, and the energy difference is `2.56e-13` Ha. The recipient's
1,189 visits each perform one Fock, one eigensolution per spin, one
post-factorization true-energy evaluation, and one publication. There are no
fallbacks or rejections. Center storage is 5,455,352 bytes; total reusable
workspace summary size is 50,909,088 bytes, versus PPP's 41,844,243 bytes.
Peak recipient RSS is 1.380 GB versus 1.724 GB in the matched PPP run.

The compact production slice is 5,091 physical lines. Warm center
preparation, Fock construction, and action allocate zero bytes in the focused
gate. Initialization is block-native and occurs once; subsequent dense-row
reads are zero.

## Validation

- H2/H10 exact center arithmetic and lifecycle gates pass.
- H20 finite/exact orientation, reflection, spin-swap, replay, ownership, and
  allocation gates pass on synthetic and frozen hydrogen fixtures.
- The matched H100 recipient/PPP timing and ownership gate passes.
- The complete direct suite passes under Julia 1.12.6.
- The complete checked-bounds suite passes under Julia 1.12.6.
- Fresh isolated `Pkg.test()` passes under Julia 1.12.6.
- `solve_hfdmrg` retains 19 methods with no ambiguity, and the public surface
  retains exactly eight exported names apart from Julia's module binding.
- Whitespace checks pass; no package Manifest, PID marker, or repository
  scratch remains.

No H1000 or 840-dimensional GaussletBases calculation was run. No claim is
made about producer-Hamiltonian stationarity or a global HF minimum. The
function-level donor mapping and honest adapted-slice hashes are recorded in
`validation/packet10a5r1b/provenance.toml`.
