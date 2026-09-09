# Packet 10A7C compact-RHF residual-transfer correction

## Scope and recovered base

This packet was prepared in the recovered linked worktree
`/Users/srw/dmrgtmp/hfdmrg-packet10a7a-recovered-20260908`, branch
`packet10a7a-recovered-c5c21cf`, from
`c5c21cf5aca72d9650387141b989a9bb8feaa489`. The accepted, unstaged Packet
10A7B occupation-default changes and six recovered Packet 10A6D documentation
deltas remain separate. The older untracked diagnostics absent at recovery
were not reconstructed.

10A7C changes compact RHF only. Compact UHF, all historical backends,
selection, optimizer, stopping policy, schedules, public APIs, and scientific
controls are unchanged.

## Ownership finding and correction

`RHFOuterBlock` stores two interaction faces, but a live left environment is
read through `right_feature`/`right_core`, and a live right environment is read
through `left_feature`/`left_core`. Fast RHF center preparation is the only
production interaction-field consumer. Root reopening, terminal
reconstruction, and `producer_energy` use state links/root coordinates rather
than interaction faces. Reflection swaps the faces.

The former shared builder always propagated both faces. During rightward
publication, its unused left/outward face crossed the entire accumulated left
block; during leftward publication, its unused right/outward face crossed the
entire accumulated right block. Initialization had the same growing-length
work. The needed inward face crosses only the newly absorbed physical group.

The builder now accepts a private `required_face` ownership mode. Every RHF
lifecycle initialization and publication site requests its inward face, so
each accepted block performs exactly one one-cell channel transfer. The unused
face is zeroed. The existing private `:both` mode remains solely as the dense
block arithmetic oracle and compatibility constructor; no compact solver
route selects it. Pair fields, H1, completed-state energy, links, edges, and
durable block shape retain their prior contracts. There is no alternate engine
or global reconstruction.

Three counters in `BlockBuildWorkspace` expose channel-transfer calls,
residual matrix steps, and maximum requested power. They cover block
initialization/publication only and do not count the fixed one- and two-cell
center contractions.

## Frozen PPP donor provenance

Authoritative donor commit:
`0bcd4f11831a98ae4bf4fb3d09eac5968df9d96c`.

- `src/PeriodicInteractionStateBlocks.jl`, `_absorb!`, lines 408--440:
  side-specific ownership and one propagation of the inherited Hartree
  channel per absorbed physical cell. Source-range SHA-256:
  `25c63f82cccaa8fa64532562a8f8a0a406559aaf4edf04dad56d4abd1570511c`.
- `src/PeriodicInteractionCenters.jl`, `_transfer_vector!`, lines 452--468:
  one-step tail multiplication and orientation-dependent residual matrix
  contraction. Source-range SHA-256:
  `e3053ac58c85e9d3c188a07296e7fda2bbd3981786a3e4a8586aebc112056d2f`.

The recipient is adapted, not byte-identical. Adaptations retain HFDMRG's
bidirectional `RHFOuterBlock`, packed pair field, result/facade contracts, and
builder workspace; select the required inward face at the six existing
lifecycle call sites; and retain `:both` only for block-oracle comparison.
No PPP module or external repository was modified.

## Permanent numerical gates

`test/compact_rhf_residual_transfer.jl` uses a 57-site synthetic fixture with
residual rank 2, a nonsymmetric residual transfer, phase 16, a partial prefix,
terminal geometry, an explicit four-atom partition, and H1 bandwidth 3. It
checks:

- one-face versus full block H1, packed interaction, inward feature/core;
- dense physical center energy and Fock, plus matrix action;
- reflection round trip and link-based reopening through repeated half-sweeps;
- exact physical projector and energy parity in both orientations;
- producer-energy parity and state nonmutation;
- invalid face, hostile capacity, inactive face zeroing, and transfer counts.

Focused ordinary and checked-bounds runs each passed 622 assertions across the
RHF foundation, lifecycle, terminal-fragment, accepted 10A7B, and new 10A7C
tests. Focused ordinary and checked-bounds nonlinear RHF runs each passed 286
assertions, including H20 finite/exact orientations and transaction paths.

## Synthetic scaling evidence

Julia 1.12.6, four BLAS threads, checked bounds, finite retained capacity 2,
12 warmed fixed-state half-sweeps:

| groups / sites | visits | new residual steps | old-formula steps | time/visit | bytes/visit |
| --- | ---: | ---: | ---: | ---: | ---: |
| 10 / 165 | 96 | 96 | 432 | 0.694 ms | 36,203 |
| 20 / 345 | 216 | 216 | 2,052 | 0.707 ms | 36,670 |

The maximum transfer power was 1 in both cases. Allocation includes state-link
factorization and exact-sized durable publication, not merely channel
propagation. The nearly unchanged per-visit cost and exact work counts establish
constant transfer-composition work per published center; the old formula grows
quadratically over a sweep.

These are deliberately synthetic fixed-rank measurements, not PPP physical
reproductions or physical-chain performance claims. No H100, H1000, huge
benchmark, full isolated suite, or stopping-policy experiment was run. The
focused exact-count and two-size evidence made a larger physical run
unnecessary at this gate.

## Packet-owned paths

- `src/interaction_blocks.jl`
- `src/fixed_state_lifecycle.jl` (10A7C call-site hunks only; also contains
  accepted 10A7B hunks)
- `test/runtests.jl` (one 10A7C include beside the accepted 10A7B include)
- `test/compact_rhf_residual_transfer.jl`
- `notes/packet10a7c_compact_rhf_residual_transfer_2026-09-08.md`
- `validation/packet10a7c/evidence.toml`
- `validation/packet10a7c/run_residual_transfer_scaling.jl`
- `validation/packet10a7c/SHA256SUMS`

All work remains unstaged and uncommitted for review.
