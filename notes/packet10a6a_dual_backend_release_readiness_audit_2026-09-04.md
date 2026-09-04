# Packet 10A6A dual-backend integration and release-readiness audit

Packet 10A6A is a source-free audit of frozen candidate commit
`698d18e98967b7e1de18a57c2f65397de2b111f8` on
`packet10a3-compact-uhf-dual-backend`.  No production or test file was changed
and no numerical suite, H1000 case, GaussletBases calculation, external test,
merge, rebase, tag, or push was run.

## Ancestry and integration state

Local `main` is `f344a698c72be2e25a2d6e3c658cdcf6ee1c870d` and is an ancestor of
the candidate.  The candidate is exactly 11 commits ahead and zero commits
behind.  A direct fast-forward is mechanically possible; no solver-history
reconciliation or merge commit is needed.  No remote or branch upstream is
configured.

The candidate worktree is intentionally not clean: the already preserved
tracked `Manifest.toml` deletion is unstaged, and the older untracked
`validation/packet10a5/h1000_rhf.toml` and
`validation/packet10a5/run_h1000_qualification.jl` remain untouched.  Their
SHA-256 hashes are respectively
`7bde224cd25b091493dbfb32e174b66e9cf361d8e31bf178a217faf0330b004a`
and
`7b5febd90ef947ba9787358d16e68570cecdf8d6991258b911daef0c4e4f8b96`.
The index is clean.

The stale primary checkout remains at `18afb6f476c1817727e59f636e4f33fb912d10f1`
on `packet10a1-compact-rhf-foundation`, with four modified tracked files and
17 untracked files.  Its eight pre-Packet-10 historical untracked files were
inspected read-only and remain distinct from the nine later Packet 10A3
transfer files:

- `JuliaStyle.md`, SHA-256 `a298dc2da12b03e1bcb19b6ae398343bc598784a59085cbebd37db4d1614c7c4`;
- `notes/angular_be_g0_continuation_result_2026-08-05.md`, `b75e14018907567bc7247c8aca97aca2753e0c59013c100f0f9b948862b7e685`;
- `notes/axis_pinned_be_history_uhf_validation_design_2026-08-05.md`, `5161d264619af97fae84815141eb5fcb045b2dea71f5a02acfb1528e2cd9e992`;
- `notes/frozen_local_solve_design_2026-08-04.md`, `a3aae93b1f64bc40ed3242225b49c71b0aac23a8f4bc6816faee232b8d993a05`;
- `notes/history_accelerated_density_uhf_basin_validation_plan_2026-08-04.md`, `76811e2d99ea8fd9dd84cbc36d7cd56f2e65b665019ccfde73fadb4fdf329e54`;
- `notes/history_accelerated_density_uhf_wl_qualification_plan_2026-08-05.md`, `98cb9d98586a58cf8606fb43760f7cf075dfb80a6defc3546a08de3c96263ef0`;
- `notes/history_accelerated_sliced_design_2026-08-05.md`, `cc3e463df2e2445e300bf41e7e7797bd3b1e9dde5386b9f1c8ad8d483a2c1bf2`;
- `notes/uhf_guidance_donor_hill_reconciliation_2026-08-05.md`, `80a7796dd2d0f46e64e9aa96b97a80818ea010167e808212cc77afd6c7005a10`.

None is a transfer-scaffolding cleanup target.  The stale primary must not be
used as the fast-forward integration checkout.

## Loaded production inventory

`src/HFDMRG.jl` contains 27 includes and the tree contains exactly 28 source
files including that module owner.  Every source file is loaded exactly once.
There is no orphan `.jl` production owner.  “Physical” below means nonblank,
non-comment-only lines; total lines include comments and blanks.

| Category | Files | Total | Physical | Ownership |
|---|---:|---:|---:|---|
| Shared module and timing | 2 | 291 | 254 | `HFDMRG.jl`, `timing.jl` |
| Compact shared records | 4 | 251 | 234 | unit-cell interaction, pairs, state links, banded one-body |
| Compact RHF | 4 | 2,498 | 2,381 | RHF entries/centers/lifecycle/nonlinear transaction |
| Compact independent-spin UHF | 4 | 2,214 | 2,124 | separated UHF entries, two-sided centers, lifecycle, nonlinear transaction |
| Compact producer observer | 1 | 519 | 496 | streamed observational producer energy |
| Historical dense core | 3 | 1,581 | 1,461 | backend API, dense sweep core, density-density backend |
| Historical sliced/target/cached | 6 | 2,490 | 2,356 | layout, projection, cached, cached Coulomb, target residual, audit |
| Historical frozen/history optional routes | 4 | 1,621 | 1,578 | frozen occupied and private RHF/UHF/sliced history drivers |
| **All loaded production** | **28** | **11,465** | **10,884** | |

Compact production totals 5,482/5,235 total/physical lines; intentionally
retained historical chemistry production totals 5,692/5,395; shared module
and timing totals 291/254.  The similar compact and historical sizes reflect
two different result and ownership contracts.  They are not duplicate
implementations of the same facade.

The only currently visible superseded transfer material is the two-file
untracked Packet 10A5 RHF diagnostic above.  It is not loaded production and
is superseded by the committed R2 record/runner.  It may be archived or
removed only after an explicit ownership decision.  Historical dense,
sliced, cached, target-residual, frozen/history files and their tests are live
chemistry capability, not cleanup candidates.

## Hot-path isolation proof

The loaded tree has no `src/uhf_blocks.jl`, `UHFOuterBlock`,
`PreparedUHFCenter`, `_prepare_center_cross!`, `_append_component_cross!`,
`_uhf_physical_block`, `build_uhf_outer_block`, or `UHFBlockBuildWorkspace`.
The only PPP references in production are donor-provenance comments in the
fast RHF and UHF center owners.  `Project.toml` has no PPP, trimmed-HFDMRG, or
GaussletBases dependency.

The typed compact method resolves to `src/nonlinear_rhf.jl`; matrix RHF and
UHF methods resolve to `src/core.jl`.  The compact lifecycle records own one
collection and one live root.  Compact results retain those records and do
not own global occupied-coefficient panels.  `producer_energy` has methods
only for `HFDMRGResult`/`UHFDMRGResult` plus `ProducerHamiltonian`; historical
tuples have no producer fallback.  Historical resource/geodesic “fallback”
labels and compact monotone geodesic fallback are internal same-engine
algorithms, not cross-backend calls.

## Capability matrix

| Capability | Dispatch/result contract | Preservation evidence |
|---|---|---|
| Compact unit-cell RHF | `BandedOneBody, UnitCellInteraction`, `spin=:rhf`; `HFDMRGResult` | R1A and R2 H100/H1000 records; dimer initialization |
| Compact unit-cell UHF | same typed method, `spin=:uhf`; `UHFDMRGResult` | R1B/R2B exact, finite, H100, and accepted H1000 evidence; atomic-Neel initialization |
| Producer energy | compact result plus `ProducerHamiltonian`; `ProducerEnergy` | observational streamed evaluator; state nonmutation; no residual/stationarity claim |
| Historical dense RHF | `H,V,C`; global `(Cup,Cdn,E)` tuple | dense tests and three-matrix GaussletBases handshake |
| Historical dense UHF | `H,V,Cup,Cdn` or `Hup,Hdn,V,Cup,Cdn`; global tuple | common/split-one-body tests and four-matrix handshake |
| Historical sliced RHF/UHF | fixed/ragged layouts and explicit backend owners; global tuple | projection and convenience-route tests |
| Cached and target residual | cached, cached-Coulomb, and density-plus-target-residual methods; global tuple | direct/cached parity, route-count, zero-residual, and split-H tests |
| Frozen/history | density frozen-occupied and private density/sliced RHF/UHF history drivers | retained source plus ordinary/checked/fresh suite coverage; deliberately unexported |

The accepted GaussletBases matrix-dispatch evidence remains applicable because
R2B changed only typed compact files.  The three-matrix RHF smoke returned
`13.02625601311897 Ha`; the one-electron four-matrix UHF smoke returned
`4.916191887700468 Ha`.  The old Pass 583 value is a separate three-electron
record with incomplete retained controls.  The governing preservation
certificate is the clean-versus-candidate two-alpha/one-beta comparison at
`39.38381696539275 Ha`, with zero energy/projector differences and identical
coefficient and projector hashes.

## Package, documentation, and release surface

- Package name/UUID/version are `HFDMRG`,
  `47a4102e-5223-4fb2-9ef6-33cce25f1843`, and `0.1.0`.
- Runtime dependencies contain only the `LinearAlgebra` stdlib.  Test extras
  are `Random`, `Serialization`, and `Test`.  There is no package `[compat]`
  table, including no declared Julia 1.12 compatibility.
- A bounded Julia 1.12.6 load resolved `pathof(HFDMRG)` to this exact
  candidate.  `solve_hfdmrg` has 19 methods, recursive ambiguity detection
  returns zero, and the eight exports are `BandedOneBody`, `HFDMRGResult`,
  `ProducerEnergy`, `ProducerHamiltonian`, `UHFDMRGResult`,
  `UnitCellInteraction`, `producer_energy`, and `solve_hfdmrg`.
- Documentation comprises 19 files with a Documenter 1 environment, strict
  doctests/checkdocs, and manual, algorithm, example, reference, and developer
  sections.  It has no deployment configuration.
- Six examples cover historical dense RHF/UHF, observer/checkpoint,
  target-residual, cached sliced RHF, and compact producer energy.  There is
  no standalone compact-UHF public example.
- Ten test files contain 4,234 lines and cover compact and historical routes.
  The final R2B direct, checked-bounds, ordinary package, and fresh isolated
  suites are already accepted; none was rerun in this audit.
- README installation, local-development, route, return-contract, numerical,
  example, and test guidance is accurate.  It explicitly says license,
  citation metadata, CI, hosted docs, and a public acceleration interface are
  separate decisions.
- `LICENSE`, `CITATION.cff`, `.github`/CI, package compat bounds, and hosted
  documentation are absent.  The tracked Manifest records Julia 1.12.2 and a
  path dependency on the package itself; its preserved unstaged deletion is
  the appropriate library-package policy but has not yet been landed.

## Release classification

Essential before fast-forward/publication:

1. Select and add a license; its terms require owner direction and cannot be
   inferred by an implementation packet.
2. Add bounded package compatibility, at minimum Julia 1.12 and the relevant
   stdlib compatibility, and retain version `0.1.0` unless release policy says
   otherwise.
3. Land the existing package `Manifest.toml` deletion and verify that clean
   package tests do not recreate repository state.
4. Add a Julia 1.12 CI workflow for an isolated package test and a strict
   Documenter build.  CI must use the repository commit, not a sibling dirty
   checkout.
5. Add one small runnable compact-UHF example and link it from README/docs so
   both typed public modes have an executable user path.
6. Revalidate the GaussletBases three-/four-matrix dispatch from clean,
   explicitly pinned source snapshots as described below, then run final
   isolated package/docs gates and require an empty integration worktree.

Optional later cleanup includes citation metadata, hosted documentation,
release tagging, benchmark automation, pruning stale administrative worktree
registrations, and consolidating old validation prose.  The two untracked
Packet 10A5 files may be archived or removed after explicit review.  None of
these optional items justifies changing solver arithmetic.

The following must not be removed: historical dense and sliced backends,
cached and cached-Coulomb owners, target-residual support, frozen-occupied
logic, private density/sliced RHF/UHF history drivers, GaussletBases-compatible
matrix signatures, tuple-return compatibility, or their regression tests and
adapters.

## Proposed Packet 10A6B sequence

1. Resolve release metadata only: owner-selected license, `[compat]`, and the
   package Manifest deletion.  Gate with TOML parsing, UUID/version checks,
   clean dependency resolution in a temporary depot, and no repository
   Manifest regeneration.
2. Add the minimal compact-UHF example and narrow README/docs link.  Gate with
   the example, documentation doctests, eight-export/19-method/zero-ambiguity
   checks, and no source diff.
3. Add CI for Julia 1.12 isolated `Pkg.test()` and strict Documenter.  Exercise
   the exact local commands once, with caches and manifests outside the tree.
4. For GaussletBases, materialize both exact Git commits into fresh `/tmp`
   snapshots (or clean detached worktrees), assert both are clean, and run
   Julia with the clean GaussletBases snapshot as the explicit active project.
   Load HFDMRG by absolute `Base.include` into a unique module, pass that module
   through the adapter's `hfmod` argument, and assert `pathof`/`functionloc`
   resolves inside the frozen HFDMRG snapshot before invoking either matrix
   handshake.  Never use `Pkg.develop` against the stale primary checkout or
   rely on ambient `LOAD_PATH`.  Compare the accepted tuple shapes, energies,
   projectors, and hashes; keep external trees read-only.
5. Run one clean direct, checked-bounds, fresh isolated package suite and docs
   build; verify checksums, whitespace, clean status, Manifest absence, and no
   source/test changes beyond the explicitly reviewed example/metadata/CI
   scope.
6. Fast-forward `main` with `--ff-only` from a new clean integration worktree,
   never from the stale primary.  Recheck the commit hash and zero-behind
   ancestry immediately before integration.  Tagging, remotes, and pushing
   remain separately authorized actions.

Public `main` should ultimately be a direct fast-forward of this lineage after
the bounded 10A6B release-hardening work.  It should not fast-forward at the
current audit gate because the license, compat, CI, Manifest policy, clean
GaussletBases integration gate, and dirty-worktree resolution are unfinished.
No numerical or architectural reconciliation is indicated.
