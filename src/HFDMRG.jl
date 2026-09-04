"""
HFDMRG implements a DMRG-style Hartree-Fock sweep.

HFDMRG retains two scientifically distinct implementations behind the public
`solve_hfdmrg` name.  Concrete argument types select the implementation; one
implementation never falls back to the other.

The compact unit-cell route owns one exact-sized block collection and one live
root, uses block-native dimer or atomic-Neel initialization, visits
terminal-complete two-block centers, and publishes at most one monotone
one-Fock update per center. Restricted and independent-spin unrestricted
states share the unit-cell operator but not orbital or exchange-pair spaces.
State selection is controlled explicitly by `state_maxdim` and
`state_cutoff`, and compact results retain block/root state without global
occupied-coefficient panels.

The historical dense/sliced route accepts global dense or sliced Hamiltonian,
interaction, and occupied-coefficient inputs.  It preserves global
coefficient output and the private history driver needed by atoms, diatomics,
and GaussletBases.

Public entrypoints:
- solve_hfdmrg(one_body::BandedOneBody, operator::UnitCellInteraction;
  kwargs...): compact unit-cell restricted or independent-spin unrestricted
  HF, selected by the `spin` keyword.
- solve_hfdmrg(H, V, psiup0; kwargs...): historical dense restricted HF.
- solve_hfdmrg(H, backend, psiup0; kwargs...): historical sliced, cached, or
  target-residual restricted HF.
- solve_hfdmrg(H, V, psiup0, psidn0; kwargs...): unrestricted HF using the
  retained historical density-density engine.
- solve_hfdmrg(Hup, Hdn, V, psiup0, psidn0; kwargs...): unrestricted HF with
  spin-dependent one-body matrices.
- producer_energy(result, producer): observational uncompressed-Hamiltonian
  energy audit for compact RHF or UHF results only.

Historical calls return `(psiup, psidn, energy)` with global occupied
coefficients. Compact calls return `HFDMRGResult` or `UHFDMRGResult` with a
compact state handle. `producer_energy` does not change either solver
implementation or reconstruct global occupied coefficients.
"""
module HFDMRG

include("timing.jl")
using .TimeG: @timeg

# Compact unit-cell RHF/UHF numerical engine.
const UNIT_CELL_WIDTH = 18
include("unit_cell_interaction.jl")
include("pairs.jl")
include("state.jl")
include("one_body.jl")
include("interaction_blocks.jl")
include("fast_rhf_centers.jl")
include("fixed_state_lifecycle.jl")
include("nonlinear_rhf.jl")
include("uhf_entries.jl")
include("uhf_centers.jl")
include("uhf_lifecycle.jl")
include("nonlinear_uhf.jl")
include("producer_energy.jl")

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

export BandedOneBody, HFDMRGResult, ProducerEnergy, ProducerHamiltonian,
    UHFDMRGResult, UnitCellInteraction, producer_energy, solve_hfdmrg

end # module
