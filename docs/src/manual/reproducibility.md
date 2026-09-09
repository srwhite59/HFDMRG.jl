# Paper-release reproducibility

HFDMRG's public workflows and its paper-facing performance evidence have
different reproducibility boundaries. Public examples and package tests are
self-contained in this repository. The large-chain extraction records also
use a frozen external fixture and must not be presented as standalone public
examples.

## Public package checks

Record the exact HFDMRG commit before running a release check:

```sh
git rev-parse HEAD
julia --version
```

From a clean checkout, resolve the library environment and run the same
isolated checked-bounds suite configured in CI:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no -e 'using Pkg; Pkg.activate(; temp=true); Pkg.develop(path=pwd()); Pkg.test("HFDMRG"; julia_args=["--check-bounds=yes"])'
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The public historical examples are `examples/quickstart.jl` for RHF and
`examples/basic_uhf.jl` for UHF. They return global
`(psiup, psidn, energy)` tuples. The public compact examples are
`examples/compact_uhf.jl` and `examples/compact_producer_energy.jl`; typed
unit-cell solves return `HFDMRGResult` or `UHFDMRGResult` and retain compact
state rather than global coefficient matrices. The remaining examples cover
observation/restart, target-residual, and cached sliced routes.

The bounded physical H10/H20 fixture is also public and self-contained:

```sh
julia --project=. examples/physical_hydrogen_chain.jl h10 rhf
julia --project=. examples/physical_hydrogen_chain.jl h10 uhf
julia --project=. examples/physical_hydrogen_chain.jl h20 rhf
julia --project=. examples/physical_hydrogen_chain.jl h20 uhf
```

The loader validates its schema, dimensions, element type, byte order, array
offsets, and per-array and whole-file SHA-256 checksums before constructing an
HFDMRG input. Public execution imports neither PPP nor GaussletBases and reads
no coefficient checkpoint.

These commands verify supported package behavior. They are not substitutes
for recording the machine, Julia version, BLAS thread count, exact controls,
and timing boundary of a performance run.

## Accepted extraction evidence

The checked-in records distinguish recipient execution from the benchmark
oracle:

| Scale | Recipient evidence | Meaning |
|---|---|---|
| H2/H10 | compact exact-arithmetic and lifecycle tests; Packet 10A5R1A/R1B notes | Small exact donor-arithmetic and ownership parity. |
| H20 | `validation/packet10a5r1b/h20_uhf.toml` plus compact RHF/UHF tests | Finite/exact orientation, reflection, spin swap, replay, and truncation-envelope checks. |
| Physical H10/H20 | `validation/packet10a7f/evidence.toml` and four per-run TOML records | Clean-checkout RHF/UHF reproduction using only the distributed primitive arrays and declared package dependencies. |
| H100 | `validation/packet10a5r1a/h100_rhf.toml` and `validation/packet10a5r1b/h100_uhf.toml` | Recipient science, per-center call counts, ownership, allocation, and matched performance. |
| H1000 RHF | `validation/packet10a5r2/h1000_rhf.toml` | Recipient dimer-start qualification: `-2227.128386923008 Ha`, nine half-sweeps, and 88.9003 seconds. |
| H1000 UHF | `validation/packet10a5r2b/h1000_uhf_qualification.toml` and `performance_reconciliation.toml` | Recipient atomic-Neel qualification: `-2282.2570732618124 Ha`, twelve half-sweeps, and 133.3405 seconds. |

Frozen PPP commit
`0bcd4f11831a98ae4bf4fb3d09eac5968df9d96c` is an external donor and
benchmark oracle, not recipient code. The RHF record compares to PPP's frozen
72.9910-second record. The final UHF comparison is contemporaneous on the same
machine: PPP took 95.0519 seconds for 10,000 visits, while HFDMRG took 133.3405
seconds for 11,989 visits. The accepted single-half ratio is 1.2421 and the
normalized per-visit ratio is 1.1701. The older
`h1000_uhf_performance_stop.toml` describes the superseded recipient entry
lifecycle and is diagnostic history, not a current performance result. No
large timing recorded for a superseded recipient path is a PPP timing.

The H1000 controls were Julia 1.12.6, four BLAS threads,
`state_maxdim=16`, `state_cutoff=1e-12`, and energy tolerance `1e-11`, with a
dimer start for RHF and atomic-Neel start for UHF. The TOML record for each run
is authoritative for initialization, half-sweep count, ranks, allocations,
workspace ownership, RSS, replay, and call counts. Timing comparisons are
machine- and load-specific and should not be reconstructed from values with a
different timing boundary.

The public H10/H20 reproductions use Julia 1.12.6 with four BLAS threads,
`state_maxdim=16`, `state_cutoff=1e-12`, requested energy tolerance `1e-11`,
physical orientation, and explicit atom intervals. RHF uses the dimer-product
start; UHF uses the atomic-Neel start. Compact results report requested and
effective stopping tolerances separately, plus the reason for any effective
Float64 floor. Projector replay changes are interpreted at their appropriate
square-root scale rather than as an energy tolerance.

## Energy boundary

Compact `represented_energy` is the energy of the compressed unit-cell
representation used by the solve. `producer_energy(result, producer)` is a
separate, streamed observation of a caller-supplied uncompressed
density-density Hamiltonian. It does not mutate compact state, certify
stationarity for the producer Hamiltonian, evaluate a producer residual, or
establish a global Hartree--Fock minimum. Historical tuple results retain their
electronic-energy contract and do not dispatch through the producer observer.
For the distributed fixture, `represented_energy` is electronic; the reported
represented total adds the recorded nuclear constant. `ProducerEnergy.total`
already includes that constant. Their difference measures the compact
interaction representation, not nonlinear convergence.

The accepted H1000 RHF record includes one nonmutating producer-energy audit.
The final H1000 UHF R2B record explicitly performed no producer audit; its
accepted qualification is represented-energy, replay, ownership, and
performance evidence only.

## Distributed physical H10/H20 fixture

`examples/data/hydrogen_chains_v1/manifest.toml` describes explicit contiguous
little-endian Float64 arrays. H10 and H20 were generated independently from
GaussletBases commit
`c680bbbcce5e7ff9ef6405c60fec0d3c2fabb2f0` using the pinned construction and
adapter provenance at frozen PPP commit
`0bcd4f11831a98ae4bf4fb3d09eac5968df9d96c`. The inputs use hydrogen STO-6G,
3.6-bohr spacing, 0.2-bohr grid spacing, 18 sites per atom, two padding cells,
open boundaries, atom-centered phase zero, one-body tolerance `1e-8`, and
interaction tolerance `1e-5`. Exact atom intervals, units, source hashes, and
reference origins are recorded in the manifest and validation evidence.

The binary payload stores the one-body bands, compact interaction entries,
the full ordered producer interaction, the producer one-body matrix, and the
nuclear constant. The full ordered interaction is deliberate: independently
evaluated opposite triangles are physically symmetric but not bitwise equal,
so triangle packing would modify the pinned producer input.

| Case | Represented total (Ha) | Producer total (Ha) | Matched frozen reference (Ha) | Half-sweeps |
|---|---:|---:|---:|---:|
| H10 RHF | -4.2015409365127425 | -4.201540681021279 | -4.2015409365128 | 9 |
| H10 UHF | -4.7512946690001066 | -4.751294570003836 | -4.751294669000901 | 12 |
| H20 RHF | -8.458420930767339 | -8.458419731998816 | -8.458420930773865 | 8 |
| H20 UHF | -9.559260520715455 | -9.55925914860868 | -9.559260520716624 | 12 |

All four solves converged under the locked controls and passed stationary or
update-disabled replay checks. The matched references certify the represented
Hamiltonian and physical basin. Producer totals are separate observational
values; they are not construction inputs and do not claim producer-Hamiltonian
stationarity or a global Hartree--Fock minimum.

## External and machine-local boundary

The large-chain runners under `validation/packet10a5r1*` and
`validation/packet10a5r2*` accept a frozen PPP root and import its
`validation/packet8f7a/run_block_native_product.jl` fixture. They also require
the matching GaussletBases environment. The frozen PPP checkout, raw logs,
PID files, serialized terminal states, and machine-local checkpoints are not
part of the HFDMRG distribution and must not be implied to be available.

The distributed H10/H20 fixture now provides an independent public compact
RHF/UHF physical reproduction. The checked-in H100/H1000 TOML records and
hashes remain auditable, but those paper-scale inputs and private checkpoints
are not distributed and their timings are not independently reproducible from
a clean HFDMRG checkout. Do not copy the private PPP module tree wholesale or
represent machine-local artifacts as public inputs.

The Packet 10A6B GaussletBases handshake record separately pins exact clean
source snapshots and matching RHF/UHF occupations. It certifies the historical
matrix/tuple compatibility boundary; it is not a compact-engine benchmark.

Local Packet 10A6B checks passed under Julia 1.12.6. A future GitHub release
must still observe the configured CI jobs from a clean checkout at the exact
candidate commit. Local success does not claim that an unrun GitHub workflow
has passed.
