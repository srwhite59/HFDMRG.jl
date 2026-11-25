# HFDMRG.jl v0.1.2 release notes

HFDMRG.jl v0.1.2 fixes three solver edge cases and clarifies two compact
diagnostics. It supports Julia 1.12.

## Solver corrections

- Historical automatic partitioning fits the requested sweep window by
  reducing the edge spans and interior chunk size when needed. Infeasible
  windows and invalid controls produce a preflight error. Existing feasible
  partitions keep their geometry.
- Compact RHF and UHF accept interactions with local terms and zero
  intercell channels.
- Compact RHF supports two partial groups with an explicit physical atom
  partition, including the two-site case. The workspace uses the real
  physical span without padding.

## Compact diagnostics

- The producer audit checks both triangles outside the declared
  `h1_bandwidth` once per call, with constant extra storage. For a lazy or
  structurally banded input with a guaranteed band, explicitly set
  `trust_h1_bandwidth=true` to skip indexing outside the band. The error
  envelope covers arithmetic on included terms.
- Replay reports `replay_root_coordinate_change` in RHF and
  `replay_alpha_root_coordinate_change` and
  `replay_beta_root_coordinate_change` in UHF. These are entrywise
  differences in the root's local bases. A basis rotation can change the
  diagnostic even for the same physical projector. The original
  `replay_*projector_change` properties and example output keys remain
  available. The replay convergence gate still uses represented energy.

## Install and run

With Julia 1.12, install the published tag into an active project:

```sh
julia -e 'using Pkg; Pkg.add(url="https://github.com/srwhite59/HFDMRG.jl", rev="v0.1.2")'
```

For the shipped examples, clone the tagged source and run from its root:

```sh
git clone --branch v0.1.2 --depth 1 https://github.com/srwhite59/HFDMRG.jl.git
cd HFDMRG.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. examples/physical_hydrogen_chain.jl h10 rhf
julia --project=. examples/physical_hydrogen_chain.jl h10 uhf
```

## Evidence and reproducibility

The distributed H10/H20 fixtures support compact RHF and UHF reproductions.
The physical example distinguishes represented electronic energy from the
producer total and nuclear offset. Paper-scale PPP runs retain separate
provenance. See the [reproducibility guide](docs/src/manual/reproducibility.md)
and the prior [v0.1.1 release](https://github.com/srwhite59/HFDMRG.jl/releases/tag/v0.1.1).
