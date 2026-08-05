# Developer Notes

HFDMRG has one sweep engine and several interaction representations. Start with
[Backend architecture](@ref) before changing a backend or adding a numerical
route. That guide defines the five-operation seam, state lifetimes, trace
convention, ownership rules, scaling expectations, and validation ladder.

## Repository map

- `src/HFDMRG.jl`: module and include order.
- `src/core.jl`: block geometry, sweep engine, local SCF, convergence, and
  observer calls.
- `src/backend_api.jl`: backend extension seam.
- `src/backends/`: density, target-residual, sliced projection, and cached
  sliced implementations.
- `src/slice_layout.jl`: fixed/ragged physical slice layout.
- `src/history_accelerated_*.jl`: private provisional history research paths.
- `examples/`: executable user workflows.
- `test/`: compact numerical and end-to-end regression tests.
- `scripts/`: engineering benchmarks and cached/projection verification.
- `reference/HF_dmrg_legacy.jl`: historical reference only; it is not an active
  regression oracle.

## Local commands

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=. scripts/verify_cached_vs_projection.jl
julia --project=docs docs/make.jl
```

The Makefile provides corresponding `make test`, `make verify`, and benchmark
conveniences. Heavy fixtures and logs belong in machine-local scratch rather
than the repository.

## Documentation policy

Public documentation should describe stable supported contracts. Private
research features may have algorithm pages when their numerical constraints
matter to maintainers, but those pages must be labeled provisional and must not
advertise a public call. Avoid undifferentiated API dumps; add a type or method
to the reference only when users need it to select or observe a supported
route.

This local Documenter site is not deployed. License, citation metadata, CI, and
a hosted documentation workflow remain separate release decisions.
