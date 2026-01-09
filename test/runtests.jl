using Test
using LinearAlgebra
using Random
using HFDMRG

old_load_path = copy(LOAD_PATH)
try
    testdir = @__DIR__
    if !(testdir in LOAD_PATH)
        pushfirst!(LOAD_PATH, testdir)
    end

    include(joinpath(testdir, "..", "reference", "HF_dmrg_legacy.jl"))
    Legacy = HF_dmrg

    @testset "Public API" begin
        @test isdefined(HFDMRG, :solve_hfdmrg)
        @test :solve_hfdmrg in names(HFDMRG, all = false)
    end

    @testset "Sliced backend stub" begin
        layout = HFDMRG.SliceLayout([2, 3])
        backend = HFDMRG.SlicedBasisBackend(layout, nothing)
        phi = zeros(2, 1)
        err = ErrorException("SlicedBasisBackend not implemented yet")
        @test_throws err HFDMRG.vee_init_block(:left, 1:2, 3:5, phi, backend)
    end

    function orthonormal_cols(rng, n, m)
        Q = Matrix(qr(randn(rng, n, m)).Q)
        Q[:, 1:m]
    end

    @testset "HF-DMRG regression" begin
        rng = MersenneTwister(1234)
        N = 16
        A = randn(rng, N, N)
        H = (A + A') / 2
        B = randn(rng, N, N)
        V = (B + B') / 2

        maxiter = 1
        blocksize = 2
        cutoff = 1e-8

        Nup = 1
        Ndn = 1
        psiup0 = orthonormal_cols(rng, N, Nup)
        psidn0 = orthonormal_cols(rng, N, Ndn)
        _, _, e_new = solve_hfdmrg(H, V, psiup0, psidn0;
            maxiter = maxiter, blocksize = blocksize, cutoff = cutoff, verbose = false)
        _, _, e_legacy = Legacy.dohfdmrg(H, V, psiup0, psidn0, false, 1,
            blocksize, maxiter, cutoff; verbose = false)
        @test isapprox(e_new, e_legacy; atol = 1e-6, rtol = 0)

        Nup_r = 1
        psiup_r = orthonormal_cols(rng, N, Nup_r)
        _, _, e_new_r = solve_hfdmrg(H, V, psiup_r;
            maxiter = maxiter, blocksize = blocksize, cutoff = cutoff, verbose = false)
        _, _, e_legacy_r = Legacy.dohfdmrg(H, V, psiup_r, psiup_r, true, 1,
            blocksize, maxiter, cutoff; verbose = false)
        @test isapprox(e_new_r, e_legacy_r; atol = 1e-6, rtol = 0)
    end
finally
    empty!(LOAD_PATH)
    append!(LOAD_PATH, old_load_path)
end
