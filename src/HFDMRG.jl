"""
HFDMRG implements a DMRG-style Hartree-Fock sweep.

The compact RHF route owns one exact-sized block collection and one live root,
uses block-native dimer initialization, visits terminal-complete two-block
centers, and publishes at most one monotone one-Fock update per center.  State
selection is controlled explicitly by `state_maxdim` and `state_cutoff`.

Public entrypoints:
- solve_hfdmrg(one_body::BandedOneBody, operator::UnitCellInteraction;
  kwargs...): compact unit-cell restricted HF.
- solve_hfdmrg(H, V, psiup0, psidn0; kwargs...): unrestricted HF using the
  retained historical density-density engine.
- solve_hfdmrg(Hup, Hdn, V, psiup0, psidn0; kwargs...): unrestricted HF with

The superseded dense, sliced, cached, and target-residual RHF signatures fail
with deterministic migration errors.  They never enter the global-coefficient
engine.  Their UHF overloads and underlying backend utilities remain available
during the staged transfer.
"""
module HFDMRG

include("timing.jl")
using .TimeG: @timeg

# Compact unit-cell RHF numerical engine.
const UNIT_CELL_WIDTH = 18
include("unit_cell_interaction.jl")
include("pairs.jl")
include("state.jl")
include("one_body.jl")
include("interaction_blocks.jl")
include("interaction_centers.jl")
include("fixed_state_lifecycle.jl")
include("nonlinear_rhf.jl")

include("backend_api.jl")
include("slice_layout.jl")
include("core.jl")
include("backends/density_density.jl")
include("history_accelerated_rhf.jl")
include("history_accelerated_uhf.jl")
include("frozen_occupied.jl")
include("backends/density_density_target_residual.jl")
include("backends/sliced_basis.jl")
include("backends/sliced_basis_cached.jl")
include("backends/sliced_basis_cached_coulomb.jl")
include("sliced_fock_audit.jl")
include("history_accelerated_sliced.jl")

export BandedOneBody, HFDMRGResult, RHFCompactState, UnitCellInteraction,
    solve_hfdmrg

end # module
