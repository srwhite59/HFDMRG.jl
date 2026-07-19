# Paper-Manager Review: Cached Sliced Milestone 0

Date: 2026-07-19
Reviewed commit: `edf398758f9df537dfe2d9e9297cd701be28c2e8`
Decision: **conditionally accepted for Milestones 1 and 2**

## Bottom Line

No blocking algebraic or scientific error was found. The proposed design is
appropriately narrow and addresses the measured defects rather than hiding
them with solver-policy changes.

The accepted core choices are:

1. retain `P`, `W`, `Wswap`, and `Vijkl` rather than introduce another cache
   family;
2. derive the absorption map from final reorthogonalized `phi_new`;
3. transform aggregate `W/Wswap`, add only newly absorbed complete-slice
   sources, and reconstruct the small `Vijkl` from `W/P`;
4. preserve direct rebuilding and projection as independent oracles and
   partial-slice fallbacks;
5. use an explicit slice-aligned partition while preserving legacy integer
   `blocksize` behavior; and
6. implement all nine ordered left/center/right Coulomb and exchange sectors
   without assuming additional raw-integral symmetry.

The leading-term reduction is a credible design estimate, not yet a timing
claim. It is acceptable that an absorption is linear in total slice count and
a complete chain remains quadratic. Window-local Fock work must become
independent of total slice count at fixed retained/center ranks; window
assembly remains a separately measured global-range cost.

## Required Clarifications Before Source Work

Record these in the maintained status memo. They do not require another design
review unless they force a wider API or algorithm change.

### 1. Define Pair Ordering

Define the matrix view explicitly:

```text
M(Vst)[(a,b),(c,d)] = Vst[a,b,c,d]
```

using the same column-major pair ordering as `_pair_map_mat!`. State the exact
transpose used for `Wswap`. This removes implementation ambiguity without
adding a public abstraction.

### 2. Freeze Partition Semantics

For the new opt-in route:

- `block_partition !== nothing` overrides the numeric `blocksize` value;
- `nblockcenter >= 1` is required for this structured route, while legacy
  behavior remains unchanged;
- the partition must cover exactly `1:N` with nonempty contiguous ranges;
- a supplied `SliceLayout` must end at `N`; and
- for a sliced backend, its dimensions/offsets must exactly match the
  partition layout. A mismatch must throw rather than silently enter global
  fallback.

The provisional spelling `block_partition=layout` is accepted. Do not add a
new partition type or backend lifecycle method.

### 3. Add Absolute Milestone-2 Gates

The scaling-ratio gate is necessary but insufficient. Also require, for the
52-by-9 rank-four engineering fixture:

```text
complete initial cache chain                <= 16 s
forward-plus-reverse absorption cache work  <= 32 s
cumulative cache-chain allocation           <= 512 MiB
```

Exclude the persistent input `V6` from the allocation figure and state exactly
what retained-state allocation is included. These are provisional engineering
gates. Missing one requires a measured explanation and renewed performance
review, not a solver-policy change.

### 4. Strengthen The Ordered-Sector Oracle

Before Milestone 3 acceptance, isolate all nine ordered sectors:

```text
LL LC LR
CL CC CR
RL RC RR
```

Use distinct unsymmetrized `Vst` and `Vts` blocks plus symmetric,
non-idempotent RHF and UHF densities. In particular, test `LC/CL`, `LR/RL`,
and `CR/RC` independently so an orientation swap cannot cancel in a full
random comparison.

### 5. Strengthen Fallback Oracles

Milestone 2 must test:

- left and right absorptions that split one physical slice, with nonzero
  old-center cross products;
- a deliberately off-span final old-row block that fails the reconstruction
  check and invokes direct rebuilding; and
- an observable count proving that the accepted aligned route invokes no
  fallback.

Milestone 3 must delegate the entire split-window RHF/UHF interaction to the
projection oracle rather than retain duplicate global exchange or Coulomb
implementations.

## Durability And Bloat

The 617-line Milestone-0 memo is acceptable as an initial mathematical design
record because it defines tensor orientation, alternatives, cost, gates, and
deletion targets. Do not grow it append-only. Future milestone updates should
replace current status sections or rely on Git history.

Before production source work, commit the currently untracked execution
contract and this acceptance review:

```text
notes/cached_sliced_performance_program_2026-07-19.md
notes/cached_sliced_performance_paper_review_2026-07-19.md
```

in a narrow documentation-only commit. Leave unrelated `JuliaStyle.md`
untracked. Do not amend `edf3987`.

The net-growth and deletion plan is accepted. Failure by Milestone 3 to remove
the duplicated RHF/UHF global sliced paths is a review stop, not merely a
missed cleanup preference.

## Authority To Continue

After recording the five clarifications and committing the execution contract,
proceed through Milestone 1 and Milestone 2 with separate narrow commits. Do
not stop after Milestone 1 if its gates pass. Stop after Milestone 2 for
paper-manager review of numerical parity, scaling, time, allocation, fallback
counts, and source growth before beginning the window-local Fock rewrite.

No production implementation or Be performance claim is accepted by this
review; only the design and the authority to execute Milestones 1 and 2 are
accepted.

-- hfdmrg-paper-manager@faraday30.ps.uci.edu
