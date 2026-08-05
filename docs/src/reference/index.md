# Reference

The reference is deliberately curated. It documents supported solve routes,
solver controls, observation, physical slice layouts, and the backend
constructors needed to select a supported interaction representation. It does
not expose every internal helper or state type.

- [Solvers and controls](@ref) contains the `solve_hfdmrg` docstrings, route
  table, keyword contract, return semantics, and `SweepInfo`.
- [Layouts and backends](@ref) contains `SliceLayout`, sliced representation
  helpers, and explicit backend constructors.

`solve_hfdmrg` is the only exported name. Supported types such as
`HFDMRG.SliceLayout` are intentionally qualified.

## Module

```@docs
HFDMRG
```
