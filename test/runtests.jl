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

    function slice_local(p, nj)
        n = (p - 1) ÷ nj + 1
        a = p - (n - 1) * nj
        n, a
    end

    function v6_lookup(V6, p, q, r, s, nj)
        np, a = slice_local(p, nj)
        nq, c = slice_local(q, nj)
        nr, b = slice_local(r, nj)
        nslice, d = slice_local(s, nj)
        if np == nr && nq == nslice
            return V6[a, b, c, d, np, nq]
        end
        0.0
    end

    @testset "Public API" begin
        @test isdefined(HFDMRG, :solve_hfdmrg)
        @test :solve_hfdmrg in names(HFDMRG, all = false)
    end

    @testset "Sliced backend stub" begin
        nj = 2
        ns = 3
        layout = HFDMRG.SliceLayout(fill(nj, ns))
        V6 = zeros(nj, nj, nj, nj, ns, ns)
        backend = HFDMRG.SlicedBasisBackend(layout, V6)
        phi = zeros(nj, 1)
        err = ErrorException("SlicedBasisBackend not implemented yet")
        @test_throws err HFDMRG.vee_init_block(:left, 1:nj, nj + 1:nj * ns, phi, backend)
    end

    @testset "Sliced RHF Fock" begin
        nj = 2
        ns = 3
        N = nj * ns
        V6 = randn(nj, nj, nj, nj, ns, ns)
        rho = randn(N, N)
        rho = (rho + rho') / 2
        F1 = zeros(N, N)
        HFDMRG.sliced_add_fock_r!(F1, rho, V6, nj, ns)

        F2 = zeros(N, N)
        for p = 1:N, q = 1:N, r = 1:N, s = 1:N
            F2[p, q] += rho[r, s] * (2.0 * v6_lookup(V6, p, q, r, s, nj) -
                                     v6_lookup(V6, p, r, q, s, nj))
        end
        @test maximum(abs.(F1 .- F2)) < 1e-10
    end

    @testset "Sliced UHF Fock" begin
        nj = 2
        ns = 3
        N = nj * ns
        V6 = randn(nj, nj, nj, nj, ns, ns)
        rhoup = randn(N, N)
        rhoup = (rhoup + rhoup') / 2
        rhodn = randn(N, N)
        rhodn = (rhodn + rhodn') / 2
        Fup1 = zeros(N, N)
        Fdn1 = zeros(N, N)
        HFDMRG.sliced_add_fock_uhf!(Fup1, Fdn1, rhoup, rhodn, V6, nj, ns)

        Fup2 = zeros(N, N)
        Fdn2 = zeros(N, N)
        for p = 1:N, q = 1:N, r = 1:N, s = 1:N
            vdir = v6_lookup(V6, p, q, r, s, nj)
            vex = v6_lookup(V6, p, r, q, s, nj)
            Fup2[p, q] += (rhoup[r, s] + rhodn[r, s]) * vdir - rhoup[r, s] * vex
            Fdn2[p, q] += (rhoup[r, s] + rhodn[r, s]) * vdir - rhodn[r, s] * vex
        end
        @test maximum(abs.(Fup1 .- Fup2)) < 1e-10
        @test maximum(abs.(Fdn1 .- Fdn2)) < 1e-10
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
