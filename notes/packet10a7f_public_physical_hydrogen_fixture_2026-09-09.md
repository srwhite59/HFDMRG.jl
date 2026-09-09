# Packet 10A7F: public physical hydrogen fixtures

## Outcome

Packet 10A7F adds independently generated H10 and H20 physical inputs for the
compact typed RHF/UHF engine. A clean temporary candidate package reproduced
all four matched represented-Hamiltonian reference basins using only files in
the HFDMRG distribution and declared Julia dependencies. Runtime execution
did not import PPP or GaussletBases, use a path callback, read a coefficient
checkpoint, or use a reference result as a construction input.

The generated coefficients, producer arrays, and retained reference results
are distributed under the HFDMRG MIT terms by explicit owner authorization.
External repositories were read-only and were not relicensed or modified.

## Pinned construction provenance

- Candidate base: `68364562d78c80da31ae351417bddc88138f155f`.
- Frozen algorithm/adapter and reference oracle: PPP
  `0bcd4f11831a98ae4bf4fb3d09eac5968df9d96c`.
- Physical producer: GaussletBases
  `c680bbbcce5e7ff9ef6405c60fec0d3c2fabb2f0`.
- Geometry: open H chain, STO-6G, 3.6 bohr separation, 0.2 bohr grid,
  18 sites per atom, two-bohr end padding, atom-centered phase zero.
- Construction tolerances: one-body `1e-8`, interaction `1e-5`.
- Energies are hartree; lengths are bohr. The manifest records explicit atom
  intervals and all source/content hashes.

H10 and H20 were constructed independently. The validation-only generator
checks the exact donor source hashes before generation. Source-parity
validation reconstructed each pinned input and compared every one-body band,
compact interaction entry, ordered producer-interaction entry, producer
one-body entry, nuclear constant, and atom interval against the loaded public
fixture.

## Format and size

The public format is `hfdmrg-physical-hydrogen-chain-v1` with contiguous
little-endian Float64 primitive arrays described by TOML. The loader verifies
schema, basename-only paths, byte order, element type, dimensions, contiguous
offsets, finite values, exact file length, atom partition, per-array hashes,
and whole-file hashes.

- H10: 183 sites, 323,976 bytes,
  SHA-256 `deec5f641062e4d96f0be6643d0f20187c9f2053c8ed5ab079bf788711999df1`.
- H20: 363 sites, 1,142,960 bytes,
  SHA-256 `28db3add9400afe663b5f113b6e8cfc005b19f4bfded624a091eec9a14ffa957`.
- Total binary payload: 1,466,936 bytes (approximately 1.40 MiB).

The full ordered producer interaction is retained. Independently evaluated
opposite triangles are physically symmetric but not bitwise identical;
triangle packing would therefore change the pinned producer input.

## Clean-checkout reproductions

All runs used Julia 1.12.6, four BLAS threads, `state_maxdim=16`,
`state_cutoff=1e-12`, requested energy tolerance `1e-11`, physical
orientation, explicit atom partitions, and the locked dimer RHF or atomic-Neel
UHF starts. Requested and effective energy tolerances were both `1e-11`, with
reason `requested_energy_tolerance`.

| Case | Represented total (Ha) | Matched reference delta (Ha) | Producer total (Ha) | Half-sweeps | Maximum rank |
|---|---:|---:|---:|---:|---:|
| H10 RHF | -4.2015409365127425 | 5.77e-14 | -4.201540681021279 | 9 | 5 |
| H10 UHF | -4.7512946690001066 | 7.94e-13 | -4.751294570003836 | 12 | 5/5 |
| H20 RHF | -8.458420930767339 | 6.53e-12 | -8.458419731998816 | 8 | 8 |
| H20 UHF | -9.559260520715455 | 1.17e-12 | -9.55925914860868 | 12 | 5/5 |

All four solves converged and passed stationary or update-disabled replay
checks. The largest absolute replay-energy change was `7.11e-15 Ha`; projector
changes were between `1.14e-11` and `2.24e-10`, within the recorded
finite-representation envelopes. Per-run TOML records contain ranks,
discarded weights, acceptance counts, replay diagnostics, allocations, and
timings.

The represented total is compact electronic energy plus the recorded nuclear
offset. Producer total is a separate streamed observation of the uncompressed
producer input and already includes the offset. The difference is a compact
interaction-representation effect, not a stopping error. No
producer-Hamiltonian stationarity or global-minimum claim is made.

## Validation

- Pinned source parity: passed for H10 and H20.
- Focused loader/schema tests: 27 assertions passed ordinarily and 27 passed
  with `--check-bounds=yes`, including checksum, trailing-byte, offset,
  element-count, partition-gap, and path-traversal rejection.
- Clean temporary package execution: four physical solves and four producer
  audits passed without external runtime imports or checkpoints.
- Reproduction record validation: matched references and case-specific
  conservative interaction-compression envelopes passed.
- The fresh isolated checked-bounds package suite passes 2,222 assertions,
  including the 27 fixture-loader assertions and both declared fixture
  stdlib dependencies (`SHA` and `TOML`). Packet 10A7F changes no solver
  source, optimizer, schedule, tolerance, API method, or historical backend.

## Public boundary and remaining release work

The H10/H20 fixture is a public physical reproduction, not a performance
benchmark. H100/H1000 inputs, raw logs, terminal checkpoints, and private PPP
paths remain machine-local and are not reproducible from the public package.
Future GitHub CI results, an exact release commit, citation metadata, tagging,
pushing, registry publication, and hosted documentation remain separate
release actions or owner decisions.
