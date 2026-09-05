# Packet 10A6B non-license release preparation

Packet 10A6B is based on accepted dual-backend candidate commit
`698d18e98967b7e1de18a57c2f65397de2b111f8` plus the separately committed
source-free Packet 10A6A audit at
`f6b2a24d3bbc96713146c2228dd54666ab20a6bc`. This packet changes no production
solver or test source and remains unstaged and uncommitted for review.

## Release metadata and public guidance

`Project.toml` retains package UUID
`47a4102e-5223-4fb2-9ef6-33cce25f1843`, version `0.1.0`, and its sole runtime
dependency `LinearAlgebra`. It now declares Julia 1.12 compatibility and the
applicable Julia 1.12/1.11 stdlib versions: `LinearAlgebra = "1.12"` and the
test-only `Random`, `Serialization`, and `Test` stdlibs at `"1.11"`. The docs
environment also requires Julia 1.12. The tracked root `Manifest.toml` deletion
is explicitly part of this library-package preparation; validation used
temporary projects and a machine-local first depot, and did not regenerate a
root Manifest.

`examples/compact_uhf.jl` is a small deterministic typed unit-cell UHF example
using the block-native atomic-Neel start. It converged in five half-sweeps to
represented energy `-0.45785804275905684`, with maximum alpha/beta ranks
`(1, 1)`, in 5.67 seconds including compilation under Julia 1.12.6. README and
manual links now
distinguish this compact `UHFDMRGResult` contract from the historical
global-coefficient tuple example.

`.github/workflows/ci.yml` adds Julia 1.12 jobs for one fresh temporary-project
package suite with `--check-bounds=yes` passed to the test subprocess, and one
strict local Documenter build. It contains no deployment, remote, tag, or
publication step. `docs/make.jl` ignores only the private `HFDMRG.TimeG`
submodule during `checkdocs=:exports`; all top-level public exports remain
strictly checked.

## Validation

All commands used Julia 1.12.6 (`julialauncher +release`) and a writable
machine-local first depot, with the normal depot available read-only after it.
The CI-equivalent commands were:

```sh
julia --startup-file=no --project=. examples/compact_uhf.jl
julia --startup-file=no -e 'using Pkg; Pkg.activate(; temp=true); Pkg.develop(path=pwd()); Pkg.test("HFDMRG"; julia_args=["--check-bounds=yes"])'
julia --startup-file=no --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=docs docs/make.jl
```

The isolated checked-bounds suite passed in approximately 263.8 seconds. The
strict Documenter build passed doctests, cross-references, exported-docstring
coverage, and rendering in 9.69 seconds after dependencies were resolved.
Generated `docs/Manifest.toml` and `docs/build/` were removed after validation.
A bounded surface check loaded HFDMRG from this exact worktree and confirmed
19 `solve_hfdmrg` methods, zero recursive ambiguities, and eight exports.

The pinned GaussletBases gate materialized exact Git archives in `/tmp`, never
loaded a dirty checkout, and used:

- GaussletBases `891bf3e84653197b90b930d7d5787be17d8fb998`;
- corrected historical HFDMRG `19e11fc2b12a142138dc3324f585d8aa92d46098`;
- candidate HFDMRG `f6b2a24d3bbc96713146c2228dd54666ab20a6bc`.

The accepted G10/count-6, 840-dimensional input hashes were reproduced. The
three-matrix RHF case used matching one-alpha/one-beta occupations and its
recorded blocksize 16, returning `13.02625601311897 Ha`. The four-matrix UHF
preservation certificate used matching two-alpha/one-beta occupations and its
recorded blocksize 64, returning `39.38381696539275 Ha`. Clean and candidate
energies differed by exactly zero; occupied-coefficient and projector hashes
were identical in both cases. Method locations resolved inside the respective
frozen `src/core.jl` snapshots. Full controls, tuple shapes, inputs, hashes, and
paths are in `validation/packet10a6b/pinned_gausslet_handshake.toml`.

The first strict docs attempt exposed that Documenter 1.19 recursively checked
the private timing submodule; the final configuration fixes that documentation
boundary without changing exports or source. An initial integration harness
attempt used the UHF blocksize for RHF and was rejected before evidence was
accepted; the retained record asserts each route's exact accepted controls and
energy. These were validation-harness corrections, not solver changes.

## Attribution and remaining decisions

HFDMRG names Steven R. White as package author, and repository history contains
only Steven/Steve White identities with two historical email spellings.
GaussletBases commit `891bf3e...` carries an MIT license naming Steven R. White.
The frozen PPP donor commit `0bcd4f11831a98ae4bf4fb3d09eac5968df9d96c`
has no `LICENSE`, `COPYING`, or `NOTICE` file. Function-level donor provenance
and hashes are retained in the Packet 10A5R1/R2B records, but provenance is not
a license grant. The concrete public-release blocker is therefore an explicit
owner decision establishing HFDMRG license terms and confirming that those
terms cover the adapted PPP donor slices. This packet does not infer or add
terms. Citation metadata and hosted documentation remain optional later
release choices; authorship metadata was not changed.

Subject to that license decision, no numerical, compatibility, CI, docs,
Manifest-policy, or integration blocker was found in the bounded preparation.
