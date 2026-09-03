# Reference

The reference is deliberately curated. It documents supported solve routes,
solver controls, observation, physical slice layouts, and the backend
constructors needed to select a supported interaction representation. It does
not expose every internal helper or state type.

- [Solvers and controls](@ref) contains the solve-route table, the distinct
  historical and compact return contracts, producer audit, keyword contract,
  and `SweepInfo`.
- [Layouts and backends](@ref) contains `SliceLayout`, sliced representation
  helpers, and explicit backend constructors.

The exported API comprises `solve_hfdmrg`, `producer_energy`,
`BandedOneBody`, `UnitCellInteraction`, `ProducerHamiltonian`,
`HFDMRGResult`, `UHFDMRGResult`, and `ProducerEnergy`. Historical types such
as `HFDMRG.SliceLayout` are intentionally qualified.

## Module

```@docs
HFDMRG
```
