# Final publication preparation audit

Date: 2026-09-13

This source-free audit covers the exact 149-commit ancestry and tracked tree at
accepted executable baseline
`5fec81573739b719b9823814a6c414ee2cbd87d2`. The private remote `main` matched
that commit before this unstaged preparation. Private hosted run `34415477615`
passed the isolated checked-bounds package job and strict documentation job in
10 minutes 29 seconds. No numerical calculation was rerun.

## Current tree and history boundary

The prior credential and large-object audit remains applicable. The current
tree and reachable history contain no private keys, tokens, credential files,
embedded authenticated URLs, raw logs, PID files, checkpoints, caches, or
unexpected large objects. The largest blobs are the authorized H20 and H10
physical fixtures at 1,142,960 and 323,976 bytes; the expected packed transfer
is about 1.31 MB.

Ordinary author emails and machine-local paths occur in Git metadata,
source-free development notes, and validation provenance. They disclose the
author's institutional email, local username, scratch conventions, and former
checkout names, but no credentials or private file contents. They are useful
for reproducing the provenance chain but are not required for ordinary package
use. The user has selected publication of the complete development history and
administrative evidence. This preparation retains it unchanged and does not
rewrite or delete history.

The transfer/process records under `notes/packet10a*` and machine-specific
runners under `validation/packet10a*` are intentionally retained with the full
history. Some require a frozen private PPP or GaussletBases checkout and cannot
execute from an ordinary public checkout; their value is provenance, not a
promise of standalone execution. The compact TOML results, hashes, public
H10/H20 fixture provenance, donor maps, and final algorithm/capability records
preserve the evidentiary chain.

Two deleted reference modules were the only unresolved ownership question in
the exact history:

- `reference/HF_dmrg_legacy.jl`, blob
  `e62cf1122b78c7240e96c333ef9fb524b2e5c928` (23,131 bytes);
- `reference/Thouless.jl`, blob
  `b0ae9a0ffa5da7cbeefaa1ea3ac4c14eaf85aba0` (26,100 bytes).

Both entered in initial scaffold commit
`9c33eb28aa08484850c3a3a3f61838889af83a9b`, have no embedded copyright,
license, or source attribution, and were deleted in commits `88c0d582...` and
`bfebe9a...`. The user has now explicitly confirmed that both are his or are
cleared for distribution. This resolves the rights blocker for retaining the
full history. The tiny deleted `test/HFnn.jl` and `test/Util.jl` compatibility
stubs and duplicated initial `backends/backend_api.jl` do not present a
separate substantive ownership concern.

Recommendation: retain the complete current source, tests, public fixtures,
license, provenance records, capability-preserving historical backends, and
the complete existing Git history. The confirmed deleted reference modules no
longer block that choice.

## GitHub Actions disclosure boundary

`.github/workflows/ci.yml` triggers on every push and pull request. It checks
out repository content, installs Julia 1.12, runs the package test suite in a
temporary environment with bounds checks, and builds strict documentation. It
has no deployment, release, artifact-upload, external path, or repository-secret
step, and does not use `pull_request_target`.

Workflow permissions are not declared explicitly and therefore inherit the
repository or organization default. Before public visibility, explicitly
limiting the token to `contents: read` is recommended as this complete,
top-level addition only:

```yaml
permissions:
  contents: read
```

The two third-party actions use mutable major-version references
(`actions/checkout@v4` and `julia-actions/setup-julia@v2`); pinning reviewed
commit SHAs is optional supply-chain hardening. Neither issue invalidates the
accepted private run. This preparation does not change the workflow.

## Citation and release material

`CITATION.cff` is deliberately software-only. Its sole author is Steven R.
White. The user explicitly confirmed that he is the sole software citation
author. The record uses repository URL
`https://github.com/srwhite59/HFDMRG.jl`, MIT license, and package version
`0.1.0`. It contains no DOI, ORCID, manuscript title, preferred paper citation,
release date, or inferred contributor.

`RELEASE_NOTES.md` describes proposed v0.1.0, the two dispatch-isolated backend
families, result contracts, public H10/H20 reproduction boundary, and accepted
private hosted baseline. It distinguishes recipient validation from frozen PPP
and paper-calculator evidence and does not attribute paper headline timings to
the recipient package. Its installation instructions distinguish adding the
library to a Julia environment from cloning the tagged repository required to
run shipped examples in place.

## Authorized next actions and optional follow-up

The user has authorized public visibility and a v0.1.0 GitHub release. The
coordinator will perform those actions only after this metadata-only release
commit passes hosted CI. The accepted run `34415477615` applies to the parent
executable baseline; this audit does not claim that final-commit CI has already
passed.

Optional later decisions are:

1. Decide whether and when paper title/DOI/authors should be added as a
   preferred citation; none is assumed here.
2. Decide whether explicit read-only Actions permission and action-SHA pinning
   are required before visibility changes.
3. Optionally pursue Zenodo archival DOI, hosted documentation, and Julia
   General registry submission after the first public release.

The live historical dense, sliced, cached, target-residual, frozen/history, and
GaussletBases-compatible routes are chemistry capabilities and are not cleanup
targets.
