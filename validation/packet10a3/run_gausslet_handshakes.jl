using LinearAlgebra
using GaussletBases

const HFDMRG_SOURCE = normpath(joinpath(@__DIR__, "..", "..", "src",
    "HFDMRG.jl"))
isfile(HFDMRG_SOURCE) || error("explicit HFDMRG source is missing")
Base.include(Main, HFDMRG_SOURCE)
const RecipientHFDMRG = Main.HFDMRG

function run_handshake(radial_ops, nelec)
    adapter = build_atomic_injected_angular_hfdmrg_hf_adapter(radial_ops;
        nelec)
    result = run_atomic_injected_angular_hfdmrg_hf(adapter;
        hfmod=RecipientHFDMRG, nblockcenter=1, maxiter=1,
        blocksize=min(size(adapter.hamiltonian, 1), 16),
        cutoff=0.0, scf_cutoff=Inf, verbose=false)
    result.psiup isa Matrix{Float64} ||
        error("GaussletBases adapter did not receive global alpha coefficients")
    result.psidn isa Matrix{Float64} ||
        error("GaussletBases adapter did not receive global beta coefficients")
    isfinite(result.energy) || error("handshake returned nonfinite energy")
    (adapter, result)
end

function main()
    length(ARGS) == 1 ||
        error("usage: run_gausslet_handshakes.jl GAUSSLET_ROOT")
    gausslet_root = realpath(ARGS[1])
    project_root = realpath(dirname(Base.active_project()))
    project_root == gausslet_root ||
        error("GaussletBases must be the explicit active project")

    rb = build_basis(RadialBasisSpec(:G10;
        count=6,
        mapping=AsinhMapping(c=0.15, s=0.15),
        reference_spacing=1.0,
        tails=3,
        odd_even_kmax=2,
        xgaussians=[XGaussian(alpha=0.2)]))
    grid = radial_quadrature(rb; accuracy=:medium, quadrature_rmax=12.0)
    radial_ops = atomic_operators(rb, grid; Z=2.0, lmax=2)

    restricted_adapter, restricted = run_handshake(radial_ops, 2)
    restricted_adapter.solver_mode == :restricted_closed_shell ||
        error("restricted fixture did not select the three-matrix route")
    restricted.psiup == restricted.psidn ||
        error("restricted handshake returned unequal spins")

    unrestricted_adapter, unrestricted = run_handshake(radial_ops, 1)
    unrestricted_adapter.solver_mode == :unrestricted ||
        error("unrestricted fixture did not select the four-matrix route")
    size(unrestricted.psiup, 2) == 1 ||
        error("unrestricted alpha occupation changed")
    size(unrestricted.psidn, 2) == 0 ||
        error("unrestricted beta occupation changed")

    println("gausslet_root=", gausslet_root)
    println("hfdmrg_source=", HFDMRG_SOURCE)
    println("basis_dim=", size(restricted_adapter.hamiltonian, 1))
    println("rhf_tuple_contract=true")
    println("rhf_energy=", restricted.energy)
    println("uhf_tuple_contract=true")
    println("uhf_energy=", unrestricted.energy)
end

main()
