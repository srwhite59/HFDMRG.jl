"""
HFDMRG implements a DMRG-style Hartree-Fock sweep.

Algorithm summary:
- Build a local L-C-R superblock basis from left/right environment blocks and
  the bare center sites.
- Perform a few SCF micro-iterations: build rho from occupied orbitals, ask the
  backend to add interaction mean-field terms to the Fock matrix, diagonalize,
  occupy, and mix densities.
- Update block bases via SVD truncation of occupied orbitals (DMRG analogy),
  and sweep left-to-right and right-to-left until energy converges.

Public entrypoints:
- solve_hfdmrg(H, V, psiup0; kwargs...): restricted HF using the density-density
  backend.
- solve_hfdmrg(H, V, psiup0, psidn0; kwargs...): unrestricted HF using the
  density-density backend.
- solve_hfdmrg(Hup, Hdn, V, psiup0, psidn0; kwargs...): unrestricted HF with
  spin-dependent one-body Hamiltonians and the same density-density backend.

Inputs:
- H, V: N x N real/symmetric matrices (one-body Hamiltonian and interaction).
- Hup, Hdn: N x N up/down one-body Hamiltonians for spin-dependent UHF.
- psiup0/psidn0: N x Nspin matrices with orthonormal columns (initial orbitals).
"""
module HFDMRG

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
include("history_accelerated_sliced.jl")

export solve_hfdmrg

end # module
