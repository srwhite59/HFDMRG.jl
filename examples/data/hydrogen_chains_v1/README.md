# Pinned physical hydrogen-chain inputs

This directory contains independently generated H10 and H20 open-chain inputs
for the compact HFDMRG RHF/UHF examples. The data are distributed under the
repository's MIT license, Copyright (c) 2026 Steven R. White.

`manifest.toml` is the authoritative text schema. Each `.f64le` file is a
contiguous sequence of the primitive `Float64` arrays listed in its case
record, in Julia column-major order and little-endian byte order. Array and
whole-file SHA-256 hashes, byte offsets, element counts, and dimensions are
checked before an HFDMRG owner is constructed. No Julia object serialization,
coefficient checkpoint, PPP source, or GaussletBases source is embedded.

The physical source is GaussletBases commit
`c680bbbcce5e7ff9ef6405c60fec0d3c2fabb2f0`, through the adapter and compressor
frozen at PPP commit `0bcd4f11831a98ae4bf4fb3d09eac5968df9d96c`.
GaussletBases is MIT licensed with the same copyright notice. PPP remains an
external provenance and benchmark oracle; public execution uses only these
data, HFDMRG, and declared Julia standard-library dependencies.

The ordered producer Vee matrices are stored in full. The pinned producer's
opposite-triangle evaluations can differ at floating-point roundoff, so
triangle packing would not preserve bitwise producer-input parity.

Run from the repository root:

```sh
julia --project=. examples/physical_hydrogen_chain.jl h10 rhf
julia --project=. examples/physical_hydrogen_chain.jl h10 uhf
julia --project=. examples/physical_hydrogen_chain.jl h20 rhf
julia --project=. examples/physical_hydrogen_chain.jl h20 uhf
```

The represented energy is electronic and belongs to the compressed compact
Hamiltonian. The producer total uses the stored ordered physical H1/Vee and
includes the nuclear-repulsion constant. Neither energy is a claim of a global
Hartree--Fock minimum, and the producer observation is not a stationarity or
residual calculation.
