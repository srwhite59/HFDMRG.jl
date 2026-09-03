# Packet 10A4: compact producer audit and public surface

Date: 2026-09-03

Status: unstaged, uncommitted review candidate based directly on Packet 10A3
commit `0c2d35dc4e054898a750ddbbcebd8921d0bbfd3e`.

## Outcome

Packet 10A4 transfers the accepted block-streamed `producer_energy` evaluator
for compact RHF and independent-spin UHF results. It preserves all historical
dense, sliced, cached, target-residual, and history routes unchanged. Concrete
argument types continue to isolate compact and historical solver dispatch, and
the producer operation has no method for historical tuple results.

The public surface is exactly `solve_hfdmrg`, `producer_energy`,
`BandedOneBody`, `UnitCellInteraction`, `ProducerHamiltonian`,
`HFDMRGResult`, `UHFDMRGResult`, and `ProducerEnergy`. Compact lifecycle state
records remain available by qualification but are not exported.

## Numerical convention and ownership

For alpha and beta projectors, the producer total is the explicit constant
plus separate one-body traces, one Hartree term from the total diagonal
density, and separate like-spin exchange terms. RHF uses equal spin
projectors. The returned decomposition and roundoff envelope are observations
of the final compact determinant; they do not replace `represented_energy` or
certify producer-Hamiltonian stationarity.

The evaluator opens exact state links outward from the live terminal root and
streams physical projector blocks. It does not reconstruct global occupied
coefficients, allocate a dense density/Fock/interaction owner, reverse the
collection, mutate lifecycle state, or call a historical backend. Time is
`O(N^2)`; storage is linear in `N` for spin densities plus cap-sized scratch.

## Compatibility

All historical RHF/UHF matrix and sliced signatures retain their global
`(psiup, psidn, energy)` tuple contract. The GaussletBases handshake, controls,
and adapter are unchanged. Typed compact RHF/UHF retains `HFDMRGResult` and
`UHFDMRGResult`, respectively, and gains only the optional producer audit.
There is no cross-backend fallback in either direction.

## Validation

Julia 1.12.6 bounded validation passed:

- exact public export set, 19 `solve_hfdmrg` methods, and zero ambiguities;
- checked-bounds dual-backend dispatch boundaries: 53/53 assertions,
  including all historical RHF/UHF route families and no producer fallback;
- compact producer example, with represented and producer energies agreeing
  to roundoff on its exact small fixture;
- checked-bounds compact foundation: 41/41 assertions;
- checked-bounds focused producer suite: 107/107 assertions;
- H2/H10/H20 RHF/UHF same-state dense-oracle decomposition;
- physical/reflected state, nonorthogonal spin, spin swap, terminal fragment,
  zero-rank, stale/hostile workspace, nonfinite input, mutation invariance,
  warmed allocation, and BLAS-preservation gates;
- representative warmed public producer evaluation: 240 bytes;
- donor source and focused-test SHA-256 equality;
- no historical backend, PPP, or GaussletBases reference in the producer
  owner; and
- whitespace check.

The documentation build was attempted but not run because Documenter is not
installed in the existing `docs` environment. No dependency installation or
repository-local manifest was created for this optional check.

No complete long numerical suite, H1000 calculation, large benchmark,
producer residual, optional operator representation, history change, external
repository mutation, merge, rebase, tag, or push was performed.

## Files

Production changes are `src/producer_energy.jl` and narrow include/export
wiring in `src/HFDMRG.jl`. Test changes are
`test/compact_producer_energy.jl`, its include in `test/runtests.jl`, the exact
public-surface assertion in `test/compact_rhf_compatibility.jl`, and test-only
`Serialization` in `Project.toml`.

README, numerical-convention/reference/manual pages, and
`examples/compact_producer_energy.jl` document and exercise the bounded public
operation while retaining historical route documentation.
