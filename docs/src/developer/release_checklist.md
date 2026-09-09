# Final release checklist

Use this checklist from a clean checkout of the exact proposed release commit.
Items marked “locally verified” describe accepted candidate evidence; they do
not replace a future clean-hosted run.

## Locally verified in the candidate lineage

- [x] MIT `LICENSE` with `Copyright (c) 2026 Steven R. White` and retained PPP
  donor provenance.
- [x] Julia 1.12 and applicable stdlib compatibility bounds; package UUID and
  version unchanged.
- [x] No package `Manifest.toml`; temporary validation environments leave no
  repository scratch.
- [x] Public historical RHF/UHF and compact RHF/UHF result contracts are
  documented and dispatch-isolated.
- [x] Nineteen `solve_hfdmrg` methods, zero recursive ambiguities, and eight
  exported names.
- [x] Isolated checked-bounds package suite and strict Documenter build passed
  locally under Julia 1.12.6.
- [x] Clean pinned GaussletBases three-/four-matrix handshake passed with
  matching occupations, tuple shapes, energies, and hashes.
- [x] H2/H10/H20/H100 recipient gates and final H1000 RHF/UHF extraction
  records are retained with controls and per-run provenance.
- [x] Distributed H10/H20 physical compact RHF/UHF fixtures pass source-array
  parity, hostile-loader tests, and clean-checkout reproductions without PPP,
  GaussletBases, or coefficient checkpoints at runtime.

## Required at the public-release boundary

- [ ] Freeze and report the exact release commit; require `main` to be an
  ancestor and use only an explicitly authorized fast-forward.
- [ ] Run the configured test and strict-docs CI jobs from a clean GitHub
  checkout. Record workflow links and the exact commit; do not infer hosted CI
  success from local validation.
- [ ] Confirm the checkout has no generated Manifest, build output, raw log,
  PID, checkpoint, or unrelated untracked file.
- [ ] Recheck README installation commands, public example links, license link,
  19-method/zero-ambiguity dispatch, and eight-name export surface.
- [ ] Describe H10/H20 as the distributable physical reproduction boundary.
  If H100/H1000 paper-scale timings are advertised as independently
  reproducible, distribute their inputs separately; until then, cite the
  checked-in records and state that private PPP fixtures/checkpoints are not
  distributed.
- [ ] Keep represented energy distinct from producer energy and avoid claims
  of producer-Hamiltonian stationarity or a global HF minimum.

## Separate owner or publication decisions

- Citation metadata and its authorship remain an owner decision.
- Tagging, pushing, hosted-documentation deployment, registry submission, and
  publication each require separate authorization.
- Paper figures and manuscript changes are owned by the separate paper-writing
  workflow, not this repository release checklist.
