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

Inputs:
- H, V: N x N real/symmetric matrices (one-body Hamiltonian and interaction).
- psiup0/psidn0: N x Nspin matrices with orthonormal columns (initial orbitals).
"""
module HFDMRG

include("backend_api.jl")
include("core.jl")
include("backends/density_density.jl")

export solve_hfdmrg

end # module
