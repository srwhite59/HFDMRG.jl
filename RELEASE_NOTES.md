# HFDMRG.jl v0.1.1 release notes

HFDMRG.jl v0.1.1 fixes compact and historical solver correctness and retains
the existing Julia 1.12 interfaces.

## Solver corrections

- Compact RHF and UHF include one-body hopping between retained environments
  when the band spans the explicit center. `BandedOneBody` validates band
  dimensions, finiteness, and trailing padding for both matrices and views,
  and copies the coefficients.
- Historical dense and sliced sweeps report the energy of the orbitals
  returned to the caller and passed to the observer. Damping remains part of
  local SCF. Sweep convergence compares measured energies, and `maxiter`
  must be a positive integer. Historical hopping ranges use exact structural
  zeros, preserving small nonzero couplings under rescaling.
- A binding compact state-rank cap selects the lowest-loss allowed
  whole-cluster choice. When the cutoff is reachable, state selection keeps
  the smallest sufficient active rank.
- Finite compact UHF applies `state_cutoff` to the discarded probability
  `1 - Pσ` independently in each spin. The existing per-spin result
  fields now report that probability instead of `1 - Pσ²`. Finite RHF
  retains the closed-shell loss `1 - P²`, and exact selection retains its
  prior behavior. The same numeric cutoff can give different finite-UHF
  ranks, energies, and certificates than v0.1.0.

## Install and run

With Julia 1.12, install the published tag into an active project:

```sh
julia -e 'using Pkg; Pkg.add(url="https://github.com/srwhite59/HFDMRG.jl", rev="v0.1.1")'
```

For the shipped examples, clone the tagged source and run from its root:

```sh
git clone --branch v0.1.1 --depth 1 https://github.com/srwhite59/HFDMRG.jl.git
cd HFDMRG.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. examples/physical_hydrogen_chain.jl h10 uhf
```

## Evidence and reproducibility

The distributed H10/H20 fixtures support compact RHF and UHF reproductions
from the public checkout. The physical example distinguishes represented
electronic energy from the producer audit and its nuclear offset. Paper-scale
PPP timings, archived calculations, and machine-local checkpoints have
separate provenance. The [reproducibility guide](docs/src/manual/reproducibility.md)
describes the public input and evidence boundary. Prior release notes are
available at [v0.1.0](https://github.com/srwhite59/HFDMRG.jl/releases/tag/v0.1.0).
