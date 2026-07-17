# HFDMRG Target-Space Interaction Residual Kickoff

## Status

This is a design-first assignment. It is not yet authority to edit source,
tests, scripts, `README.md`, or `Project.toml`.

The first deliverable is a concise design memo. Implementation begins only
after user review.

## Motivation

The established HFDMRG density-density route accepts a diagonal gausslet
interaction `V[i,j]`. The July 17 CR2 experiment demonstrated that it is useful
to add a small exact four-index residual acting in a fixed target orbital
space:

```text
large working basis:       N = 7069
target occupied union:     m = 46
unique target pairs:       1081
```

The CR2-local solver projects each moving-window density into the 46-orbital
space, contracts the residual direct and exchange terms there, and lifts the
Fock correction back into the moving window. The eight-sweep calculation was
stable and improved the screened Cr2 HF energy by about `0.439 mHa`, but it had
to copy the HFDMRG sweep because HFDMRG has no interaction hook for this mixed
representation.

The desired durable capability is:

```text
density-density V[i,j]
+ a limited target-space four-index residual
```

It must not become a global `N^4` interaction interface.

## Required Reading

Repository:

- `AGENTS.md`
- `README.md`
- `JuliaStyle.md`
- `src/HFDMRG.jl`
- `src/backend_api.jl`
- `src/backends/density_density.jl`
- relevant entrypoint and sweep code in `src/core.jl`

Archive/history:

- `~/Dropbox/chatarchive/handoff/software_packets/hfdmrg.md`
- `~/Dropbox/chatarchive/codex_sessions/capsules/2026/2026-03-07__look-at-hf-dmrg-jl__019cca40.md`
- `~/Dropbox/chatarchive/reports/software_reviews/qiu_white_atomic_hf_exactification_review_for_cartesian_chromium_2026-04-14.md`

CR2 executable oracle:

- `~/Dropbox/codexhome/work/cr2/reports/cr2_q7_occupied_space_four_index_residual_2026-07-17.md`
- `~/Dropbox/codexhome/work/cr2/runs/2026-07-17_cr2_occupied_residual/run_screened_occupied_residual_hfdmrg.jl`
- especially `pair_to_full`, `jk`, `residual_state`,
  `add_residual_fock!`, and the construction of `Qwindow`

Do not modify the CR2 files.

## Existing Architectural Seam

The HFDMRG sweep core is already interaction-agnostic. Backends own:

- block initialization;
- block absorption;
- moving-window interaction state;
- RHF and UHF Fock additions.

The preferred direction is therefore a narrow backend or additive backend
component. In a moving window, the demonstrated operation is:

```text
D_target = Qwindow' * D_window * Qwindow
(J_target, K_target) = contract(residual, D_target)
DeltaF_window = Qwindow * (J_target - K_target) * Qwindow'
```

The design should determine how `Qwindow` is maintained through block
projection without copying the sweep and without forming a global four-index
tensor.

## First Assignment: Read-Only Design Audit

Write:

`notes/target_space_interaction_residual_design_2026-07-17.md`

Do not edit any other file.

The memo must answer:

1. What is the smallest concrete HFDMRG data contract?
   - global-to-target map `Q`;
   - target-space interaction residual;
   - pair ordering, symmetry, normalization, and scalar type;
   - dense pair matrix versus full target tensor versus factorized channels.

2. Should the first implementation be a specific
   density-density-plus-target-residual backend, or a generic additive backend?
   Prefer the specific form unless a second live consumer justifies generic
   composition.

3. How should target vectors be represented in left block, center, right
   block, and moving-window coordinates?

4. Give explicit RHF and UHF direct, exchange, Fock, and energy formulas with
   the repository's index convention. State all factors of two.

5. How will the implementation guarantee that a zero residual reproduces the
   existing density-density route exactly?

6. What are the expected costs in `N`, target rank `m`, window dimension, and
   residual factor rank? Estimate the Cr2 memory footprint and per-window
   contraction cost for `N=7069`, `m=46`, and block size `128`.

7. Is factorization needed for the first correct implementation? Inspect the
   saved residual spectrum or benchmarkable CR2 representation if practical;
   do not assume compression without evidence.

8. What is the smallest validation ladder?
   - tiny random symmetric-orbital-integral oracle;
   - RHF and UHF Fock/energy parity against explicit transformed four-index
     integrals;
   - zero-residual regression;
   - parity with the saved CR2 fixed-state correction;
   - one representative HFDMRG sweep or accepted restart;
   - timing and allocation/memory comparison.

9. Which existing code becomes unnecessary after adoption? The CR2 copied
   sweep should be retired only after production parity is demonstrated.

10. List exact proposed files, public names, and an added-line budget. Avoid a
    broad framework, new artifact schema, or unrelated cleanup.

## Guardrails

- No source or test edits in the first pass.
- No copied or forked sweep engine.
- No global `N x N x N x N` storage.
- No assumption that the residual is positive semidefinite; it is a signed
  difference between exact and approximate interactions.
- No truncation of factors without measured energy and Fock error.
- Preserve RHF, common-H UHF, split-H UHF, and sliced-backend behavior.
- Do not change GaussletBases or CR2.
- Preserve the pre-existing untracked `JuliaStyle.md`.
- Performance is part of acceptance.

## First Handback

Report:

- the recommended minimal representation and backend boundary;
- the exact Fock/energy convention;
- whether factorization appears worthwhile;
- expected Cr2 cost;
- proposed implementation slices and line budget;
- ambiguities requiring user judgment;
- files changed, which must be only the requested design memo;
- final git status.

Sign:

`-- hfdmrg-manager@<host>`
