# Paper-Manager Review: Cached Sliced Milestones 1 And 2

Date: 2026-07-19
Reviewed commits: `8450983`, `31a9c76`
Decision: **hold Milestone 3 pending one bounded corrective commit**

## Bottom Line

The slice-aligned partition and incremental cache recurrence are numerically
and architecturally successful. No error was found in the active-source direct
builder, ordered `W/Wswap` update, `(1,3,2,4)` `Vijkl` permutation, or the
fixed/ragged left/right contractions. The measured cache work is smaller by
orders of magnitude and comfortably passes the provisional time, allocation,
and scaling gates.

Milestone 3 is temporarily held because independent review exposed one
nonportable committed test, one pre-existing multi-center core defect now
admitted by the new public partition route, and two narrow guard/test gaps.
These do not presently indicate a cache-algebra defect.

## Accepted Evidence

The reported cache benchmark includes owned `phi`, `P`, all retained
`W/Wswap`, rebuilt `Vijkl`, pair maps, temporaries, and retained output states.
It appropriately excludes the persistent input `V6`, precomputed basis/SVD
plans, one-body work, window construction, and Fock evaluation.

For the synthetic 52-by-9 rank-four fixture:

```text
initial cache chain                 0.042467 s,  83.419 MiB
98 absorption cache updates         0.082025 s, 161.708 MiB
combined allocation                              245.127 MiB
```

The `8 -> 16` and `16 -> 32` time ratios of `4.416` and `4.365` support the
bounded quadratic-chain claim and pass the `<=5` gate. These remain
engineering results, not complete-solver or publication scaling claims.

The implementation stayed inside the M1/M2 source and test budgets and
deleted the global scalar cache-rebuild loops. Milestone 3 must still delete
the duplicated global RHF/UHF Fock paths materially.

## Required Corrective Pass

Make one new corrective commit. Do not amend the accepted milestone commits.
Do not begin the window-local Fock rewrite in this pass.

### 1. Replace The Nonportable End-To-End Fixture

An independent faraday run of `julia --project=. test/runtests.jl` failed:

```text
common-H aligned projection/cached difference  5.781764e-9 Ha
committed tolerance                            1.0e-9 Ha

split-H aligned projection/cached difference   3.093631e-9 Ha
committed tolerance                            3.0e-9 Ha
```

Direct cache and Fock parity remains near `1e-15`, so the likely cause is
platform-sensitive nonlinear amplification in a random, poorly conditioned
trajectory fixture.

Diagnose the first diverging window. Then replace the fixture with a
deterministic, well-gapped problem and a justified cross-platform outcome
gate. Preserve the machine-precision cache/Fock oracles. Do not merely loosen
the current random-energy tolerances without diagnosing them. Include an
occupied-projector or equivalent physical-state comparison if it is stable
and useful.

Rerun the complete suite on macmini. The paper manager will rerun it on
faraday before releasing Milestone 3. Correct the status memo's unqualified
statement that all 89 tests pass.

### 2. Fix The General `nblockcenter` Offset

The forward sweep currently locates the old right block using only one
remaining center size when `nblockcenter > 1`. It is correct for `d=1,2` but
wrong for `d>=3`.

In both common-H and split-H cores, the right-block start after absorbing the
first center must include:

```text
m + sum(blocksizes[b+2:b+d]) + 1
```

for `d = nblockcenter`, with the empty sum giving the existing `d=1` result.
Add one valid `d=3` regression with nonuniform block sizes that would fail the
old expression. Prefer an independent returned-energy/state consistency check
over a private index-vocabulary test.

This is a pre-existing core bug exposed during the new partition review; its
small repair belongs in the corrective commit.

### 3. Enforce Side-Anchored Incremental Eligibility

The backend contract says a left block is `1:l` and a right block is `r:N`.
Before entering the incremental path, require:

```text
side == :left   => first(oldra) == 1
side == :right  => last(oldra) == N
```

An internal block omitting a prefix or suffix must use direct rebuilding.
Cover this guard compactly; do not create a new eligibility framework.

### 4. Complete The Off-Span Fallback Coverage

The explicit off-span reconstruction test currently covers only left
absorption. Add the corresponding right-side case. Keep it within the existing
fallback fixture rather than creating another broad test family.

### 5. Clarify Non-Sliced Partition Semantics

The implementation currently permits `block_partition=layout` for
density-density backends and requires exact layout matching only when the
backend itself owns a `SliceLayout`. This is a coherent and potentially useful
general block-partition feature.

State that behavior explicitly in the README and add one compact public
density-density solve whose supplied `SliceLayout` reproduces known legacy
ranges/results. Do not add another partition type. Additional exhaustive
fixed/ragged/RHF/UHF combinations are not required.

## Benchmark Provenance

The exact `0.042467 s` and `0.082025 s` measurements were made before the final
commit and their temporary harness was deleted. Their large margin and detailed
methodology are sufficient for this engineering checkpoint, but they are not
durable paper evidence.

When the existing `scripts/bench_sliced.jl` is updated for Milestone 3, include
the cache-chain measurement and rerun it on the exact pinned commit. Final
paper-facing timing must come from that retained harness. Do not add a second
temporary benchmark framework merely to preserve these interim numbers.

## Scope And Next Stop

Keep the corrective pass small: direct source fixes plus replacement/extension
of existing tests and concise README/status corrections. Do not add a source
file, public type, backend lifecycle method, or compatibility wrapper.

Commit this review file with the corrective pass. Leave unrelated
`JuliaStyle.md` untracked, and do not amend the milestone commits.

After the corrective commit and clean macmini suite, stop for paper-manager
cross-host verification. If the faraday suite passes and no new numerical
issue appears, Milestones 1 and 2 will be accepted and Milestone 3 released
without another design round.

-- hfdmrg-paper-manager@faraday30.ps.uci.edu
