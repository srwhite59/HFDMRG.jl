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

    function vblocks_lookup(layout, Vblocks, p, q, r, s)
        np, a = HFDMRG.slice_local(layout, p)
        nq, c = HFDMRG.slice_local(layout, q)
        nr, b = HFDMRG.slice_local(layout, r)
        nslice, d = HFDMRG.slice_local(layout, s)
        if np == nr && nq == nslice
            return Vblocks[np][nq][a, b, c, d]
        end
        0.0
    end

    function orthonormal_cols(rng, n, m)
        Q = Matrix(qr(randn(rng, n, m)).Q)
        Q[:, 1:m]
    end

    @testset "Public API" begin
        @test isdefined(HFDMRG, :solve_hfdmrg)
        @test :solve_hfdmrg in names(HFDMRG, all = false)
    end

    @testset "Split one-body UHF API" begin
        rng = MersenneTwister(301)
        N = 12
        A = randn(rng, N, N)
        H = (A + A') / 2
        B = randn(rng, N, N)
        V = 0.01 * ((B + B') / 2)
        Nup = 1
        Ndn = 1
        psiup0 = orthonormal_cols(rng, N, Nup)
        psidn0 = orthonormal_cols(rng, N, Ndn)

        _, _, e_old = solve_hfdmrg(H, V, psiup0, psidn0;
            maxiter = 2, blocksize = 2, cutoff = 1e-8, verbose = false)
        _, _, e_same = solve_hfdmrg(H, H, V, psiup0, psidn0;
            maxiter = 2, blocksize = 2, cutoff = 1e-8, verbose = false)
        @test e_same == e_old
        _, _, e_scf_cutoff = solve_hfdmrg(H, V, psiup0, psidn0;
            maxiter = 2, blocksize = 2, cutoff = 1e-8, scf_cutoff = 1e-8,
            verbose = false)
        @test e_scf_cutoff == e_old

        diag_up = [3.0, -2.0, 1.0, 4.0, 0.5, 2.0, 6.0, 5.0]
        diag_dn = [1.0, 4.0, 3.0, 0.25, -1.5, 2.0, 5.0, 6.0]
        Hup = Matrix(Diagonal(diag_up))
        Hdn = Matrix(Diagonal(diag_dn))
        V0 = zeros(length(diag_up), length(diag_up))
        psiup_exact = zeros(length(diag_up), 1)
        psiup_exact[argmin(diag_up), 1] = 1.0
        psidn_exact = zeros(length(diag_dn), 1)
        psidn_exact[argmin(diag_dn), 1] = 1.0
        _, _, e_split = solve_hfdmrg(Hup, Hdn, V0, psiup_exact, psidn_exact;
            maxiter = 2, blocksize = 1, cutoff = 1e-12, verbose = false)
        @test isapprox(e_split, minimum(diag_up) + minimum(diag_dn); atol = 1e-10,
            rtol = 0)
    end

    @testset "Sliced backend stub" begin
        nj = 2
        ns = 3
        N = nj * ns
        layout = HFDMRG.SliceLayout(fill(nj, ns))
        V6 = zeros(nj, nj, nj, nj, ns, ns)
        backend = HFDMRG.SlicedBasisBackend(layout, V6)
        phi = orthonormal_cols(MersenneTwister(1), nj, 1)
        state = HFDMRG.vee_init_block(:left, 1:nj, (nj + 1):N, phi, backend)
        @test state.ra == 1:nj
        @test size(state.phi) == size(phi)
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

    @testset "Sliced ragged RHF Fock" begin
        rng = MersenneTwister(202)
        dims = [1, 2, 3]
        layout = HFDMRG.SliceLayout(dims)
        ns = length(dims)
        N = layout.offs[end]
        Vblocks = [[randn(rng, dims[n], dims[n], dims[m], dims[m]) for m in 1:ns]
                   for n in 1:ns]
        rho = randn(rng, N, N)
        rho = (rho + rho') / 2
        vee = HFDMRG.SlicedVeeRagged(layout, Vblocks)

        F1 = zeros(N, N)
        HFDMRG.sliced_add_fock_r!(F1, rho, vee)
        F2 = zeros(N, N)
        for p = 1:N, q = 1:N, r = 1:N, s = 1:N
            F2[p, q] += rho[r, s] * (2.0 * vblocks_lookup(layout, Vblocks, p, q, r, s) -
                                     vblocks_lookup(layout, Vblocks, p, r, q, s))
        end
        @test maximum(abs.(F1 .- F2)) < 1e-10
    end

    @testset "Sliced ragged UHF Fock" begin
        rng = MersenneTwister(203)
        dims = [1, 2, 3]
        layout = HFDMRG.SliceLayout(dims)
        ns = length(dims)
        N = layout.offs[end]
        Vblocks = [[randn(rng, dims[n], dims[n], dims[m], dims[m]) for m in 1:ns]
                   for n in 1:ns]
        rhoup = randn(rng, N, N)
        rhoup = (rhoup + rhoup') / 2
        rhodn = randn(rng, N, N)
        rhodn = (rhodn + rhodn') / 2
        vee = HFDMRG.SlicedVeeRagged(layout, Vblocks)

        Fup1 = zeros(N, N)
        Fdn1 = zeros(N, N)
        HFDMRG.sliced_add_fock_uhf!(Fup1, Fdn1, rhoup, rhodn, vee)
        Fup2 = zeros(N, N)
        Fdn2 = zeros(N, N)
        for p = 1:N, q = 1:N, r = 1:N, s = 1:N
            vdir = vblocks_lookup(layout, Vblocks, p, q, r, s)
            vex = vblocks_lookup(layout, Vblocks, p, r, q, s)
            Fup2[p, q] += (rhoup[r, s] + rhodn[r, s]) * vdir - rhoup[r, s] * vex
            Fdn2[p, q] += (rhoup[r, s] + rhodn[r, s]) * vdir - rhodn[r, s] * vex
        end
        @test maximum(abs.(Fup1 .- Fup2)) < 1e-10
        @test maximum(abs.(Fdn1 .- Fdn2)) < 1e-10
    end

    @testset "Sliced window projection" begin
        rng = MersenneTwister(42)
        nj = 2
        ns = 4
        N = nj * ns
        layout = HFDMRG.SliceLayout(fill(nj, ns))
        V6 = randn(rng, nj, nj, nj, nj, ns, ns)
        backend = HFDMRG.SlicedBasisBackend(layout, V6)

        Lra = 1:2
        Cra = 3:6
        Rra = 7:8
        Lphi = orthonormal_cols(rng, length(Lra), 1)
        Rphi = orthonormal_cols(rng, length(Rra), 1)
        Lvee = HFDMRG.vee_init_block(:left, Lra, (Lra[end] + 1):N, Lphi, backend)
        Rvee = HFDMRG.vee_init_block(:right, Rra, 1:(Rra[1] - 1), Rphi, backend)
        win = HFDMRG.vee_window(Lvee, Rvee, Cra, backend)
        B = win.B
        superdim = size(B, 2)

        rho_sb = randn(rng, superdim, superdim)
        rho_sb = (rho_sb + rho_sb') / 2
        F_sb = zeros(superdim, superdim)
        HFDMRG.vee_add_fock_r!(F_sb, rho_sb, win)
        rho_full = B * rho_sb * B'
        G_full = zeros(N, N)
        HFDMRG.sliced_add_fock_r!(G_full, rho_full, V6, nj, ns)
        G_sb_ref = B' * G_full * B
        @test maximum(abs.(F_sb .- G_sb_ref)) < 1e-10

        rhoup_sb = randn(rng, superdim, superdim)
        rhoup_sb = (rhoup_sb + rhoup_sb') / 2
        rhodn_sb = randn(rng, superdim, superdim)
        rhodn_sb = (rhodn_sb + rhodn_sb') / 2
        Fup_sb = zeros(superdim, superdim)
        Fdn_sb = zeros(superdim, superdim)
        HFDMRG.vee_add_fock!(Fup_sb, Fdn_sb, rhoup_sb, rhodn_sb, win)
        rhoup_full = B * rhoup_sb * B'
        rhodn_full = B * rhodn_sb * B'
        Gup_full = zeros(N, N)
        Gdn_full = zeros(N, N)
        HFDMRG.sliced_add_fock_uhf!(Gup_full, Gdn_full, rhoup_full, rhodn_full, V6, nj, ns)
        Gup_sb_ref = B' * Gup_full * B
        Gdn_sb_ref = B' * Gdn_full * B
        @test maximum(abs.(Fup_sb .- Gup_sb_ref)) < 1e-10
        @test maximum(abs.(Fdn_sb .- Gdn_sb_ref)) < 1e-10
    end

    @testset "Cached sliced window parity" begin
        rng = MersenneTwister(31)
        nj = 2
        ns = 4
        N = nj * ns
        layout = HFDMRG.SliceLayout(fill(nj, ns))
        V6 = randn(rng, nj, nj, nj, nj, ns, ns)
        backend_proj = HFDMRG.SlicedBasisBackend(layout, V6)
        backend_cached = HFDMRG.SlicedBasisBackendCached(layout, V6)

        Lra = 1:2
        Cra = 3:6
        Rra = 7:8
        Lphi = orthonormal_cols(rng, length(Lra), 1)
        Rphi = orthonormal_cols(rng, length(Rra), 1)
        Lvee_proj = HFDMRG.vee_init_block(:left, Lra, (Lra[end] + 1):N, Lphi, backend_proj)
        Rvee_proj = HFDMRG.vee_init_block(:right, Rra, 1:(Rra[1] - 1), Rphi, backend_proj)
        win_proj = HFDMRG.vee_window(Lvee_proj, Rvee_proj, Cra, backend_proj)

        Lvee_cached = HFDMRG.vee_init_block(:left, Lra, (Lra[end] + 1):N, Lphi, backend_cached)
        Rvee_cached = HFDMRG.vee_init_block(:right, Rra, 1:(Rra[1] - 1), Rphi, backend_cached)
        win_cached = HFDMRG.vee_window(Lvee_cached, Rvee_cached, Cra, backend_cached)

        superdim = size(Lphi, 2) + length(Cra) + size(Rphi, 2)
        rho = randn(rng, superdim, superdim)
        F_proj = zeros(superdim, superdim)
        F_cached = zeros(superdim, superdim)
        HFDMRG.vee_add_fock_r!(F_proj, rho, win_proj)
        HFDMRG.vee_add_fock_r!(F_cached, rho, win_cached)
        @test maximum(abs.(F_proj .- F_cached)) < 1e-10

        rhoup = randn(rng, superdim, superdim)
        rhodn = randn(rng, superdim, superdim)
        Fup_proj = zeros(superdim, superdim)
        Fdn_proj = zeros(superdim, superdim)
        Fup_cached = zeros(superdim, superdim)
        Fdn_cached = zeros(superdim, superdim)
        HFDMRG.vee_add_fock!(Fup_proj, Fdn_proj, rhoup, rhodn, win_proj)
        HFDMRG.vee_add_fock!(Fup_cached, Fdn_cached, rhoup, rhodn, win_cached)
        @test maximum(abs.(Fup_proj .- Fup_cached)) < 1e-10
        @test maximum(abs.(Fdn_proj .- Fdn_cached)) < 1e-10

        rho2 = randn(rng, superdim, superdim)
        F1 = zeros(superdim, superdim)
        F2 = zeros(superdim, superdim)
        HFDMRG.vee_add_fock_r!(F1, rho2, win_cached)
        HFDMRG.vee_add_fock_r!(F2, rho2, win_cached)
        @test maximum(abs.(F1 .- F2)) < 1e-12

        rhoup2 = randn(rng, superdim, superdim)
        rhodn2 = randn(rng, superdim, superdim)
        Fup1 = zeros(superdim, superdim)
        Fdn1 = zeros(superdim, superdim)
        Fup2 = zeros(superdim, superdim)
        Fdn2 = zeros(superdim, superdim)
        HFDMRG.vee_add_fock!(Fup1, Fdn1, rhoup2, rhodn2, win_cached)
        HFDMRG.vee_add_fock!(Fup2, Fdn2, rhoup2, rhodn2, win_cached)
        @test maximum(abs.(Fup1 .- Fup2)) < 1e-12
        @test maximum(abs.(Fdn1 .- Fdn2)) < 1e-12
    end

    @testset "Sliced end-to-end sweep" begin
        rng = MersenneTwister(7)
        nj = 2
        ns = 4
        N = nj * ns
        layout = HFDMRG.SliceLayout(fill(nj, ns))
        H = randn(rng, N, N)
        H = (H + H') / 2
        V6 = 0.01 * randn(rng, nj, nj, nj, nj, ns, ns)
        backend = HFDMRG.SlicedBasisBackend(layout, V6)
        Nup = 2
        psiup0 = orthonormal_cols(rng, N, Nup)
        rho0 = psiup0 * psiup0'
        F0 = copy(H)
        HFDMRG.sliced_add_fock_r!(F0, rho0, V6, nj, ns)
        energy0 = tr(rho0 * (F0 + H))
        _, _, energy = solve_hfdmrg(H, backend, psiup0;
            maxiter = 2, blocksize = 2, cutoff = 1e-8, verbose = false)
        @test isfinite(energy)
        @test energy <= energy0 + 1e-6 * max(1.0, abs(energy0))
    end

    @testset "Sliced ragged end-to-end sweep" begin
        rng = MersenneTwister(91)
        dims = [2, 3, 2, 3]
        layout = HFDMRG.SliceLayout(dims)
        ns = length(dims)
        N = layout.offs[end]
        H = randn(rng, N, N)
        H = (H + H') / 2
        Vblocks = [[0.01 * randn(rng, dims[n], dims[n], dims[m], dims[m]) for m in 1:ns]
                   for n in 1:ns]
        vee = HFDMRG.SlicedVeeRagged(layout, Vblocks)
        backend = HFDMRG.SlicedBasisBackend(vee)
        Nup = 2
        psiup0 = orthonormal_cols(rng, N, Nup)
        rho0 = psiup0 * psiup0'
        F0 = copy(H)
        HFDMRG.sliced_add_fock_r!(F0, rho0, vee)
        energy0 = tr(rho0 * (F0 + H))
        _, _, energy = solve_hfdmrg(H, backend, psiup0;
            maxiter = 2, blocksize = 2, cutoff = 1e-8, verbose = false)
        @test isfinite(energy)
        @test energy <= energy0 + 1e-6 * max(1.0, abs(energy0))
    end

    @testset "Cached sliced ragged window parity" begin
        rng = MersenneTwister(92)
        dims = [1, 2, 3]
        layout = HFDMRG.SliceLayout(dims)
        ns = length(dims)
        N = layout.offs[end]
        Vblocks = [[randn(rng, dims[n], dims[n], dims[m], dims[m]) for m in 1:ns]
                   for n in 1:ns]
        backend_proj = HFDMRG.SlicedBasisBackend(layout, Vblocks)
        backend_cached = HFDMRG.SlicedBasisBackendCached(layout, Vblocks)

        Lra = 1:2
        Cra = 3:4
        Rra = 5:6
        Lphi = orthonormal_cols(rng, length(Lra), 1)
        Rphi = orthonormal_cols(rng, length(Rra), 1)
        Lvee_proj = HFDMRG.vee_init_block(:left, Lra, (Lra[end] + 1):N, Lphi, backend_proj)
        Rvee_proj = HFDMRG.vee_init_block(:right, Rra, 1:(Rra[1] - 1), Rphi, backend_proj)
        win_proj = HFDMRG.vee_window(Lvee_proj, Rvee_proj, Cra, backend_proj)

        Lvee_cached = HFDMRG.vee_init_block(:left, Lra, (Lra[end] + 1):N, Lphi, backend_cached)
        Rvee_cached = HFDMRG.vee_init_block(:right, Rra, 1:(Rra[1] - 1), Rphi, backend_cached)
        win_cached = HFDMRG.vee_window(Lvee_cached, Rvee_cached, Cra, backend_cached)

        superdim = size(Lphi, 2) + length(Cra) + size(Rphi, 2)
        rho = randn(rng, superdim, superdim)
        rho = (rho + rho') / 2
        F_proj = zeros(superdim, superdim)
        F_cached = zeros(superdim, superdim)
        HFDMRG.vee_add_fock_r!(F_proj, rho, win_proj)
        HFDMRG.vee_add_fock_r!(F_cached, rho, win_cached)
        @test maximum(abs.(F_proj .- F_cached)) < 1e-10

        rhoup = randn(rng, superdim, superdim)
        rhoup = (rhoup + rhoup') / 2
        rhodn = randn(rng, superdim, superdim)
        rhodn = (rhodn + rhodn') / 2
        Fup_proj = zeros(superdim, superdim)
        Fdn_proj = zeros(superdim, superdim)
        Fup_cached = zeros(superdim, superdim)
        Fdn_cached = zeros(superdim, superdim)
        HFDMRG.vee_add_fock!(Fup_proj, Fdn_proj, rhoup, rhodn, win_proj)
        HFDMRG.vee_add_fock!(Fup_cached, Fdn_cached, rhoup, rhodn, win_cached)
        @test maximum(abs.(Fup_proj .- Fup_cached)) < 1e-10
        @test maximum(abs.(Fdn_proj .- Fdn_cached)) < 1e-10
    end

    @testset "Cached ragged backend vs projection" begin
        rng = MersenneTwister(93)
        dims = [1, 2, 3, 2, 2]
        layout = HFDMRG.SliceLayout(dims)
        ns = length(dims)
        N = layout.offs[end]
        H = randn(rng, N, N)
        H = (H + H') / 2
        Vblocks = [[0.01 * randn(rng, dims[n], dims[n], dims[m], dims[m]) for m in 1:ns]
                   for n in 1:ns]
        backend_proj = HFDMRG.SlicedBasisBackend(layout, Vblocks)
        backend_cached = HFDMRG.SlicedBasisBackendCached(layout, Vblocks)
        Nup = 1
        Ndn = 1
        psiup0 = orthonormal_cols(rng, N, Nup)
        psidn0 = orthonormal_cols(rng, N, Ndn)
        _, _, e_proj = solve_hfdmrg(H, backend_proj, psiup0, psidn0;
            maxiter = 2, blocksize = 2, cutoff = 1e-8, verbose = false)
        _, _, e_cached = solve_hfdmrg(H, backend_cached, psiup0, psidn0;
            maxiter = 2, blocksize = 2, cutoff = 1e-8, verbose = false)
        @test isapprox(e_proj, e_cached; atol = 1e-9, rtol = 0)

        psiup0_r = orthonormal_cols(rng, N, Nup)
        _, _, e_proj_r = solve_hfdmrg(H, backend_proj, psiup0_r;
            maxiter = 2, blocksize = 2, cutoff = 1e-8, verbose = false)
        _, _, e_cached_r = solve_hfdmrg(H, backend_cached, psiup0_r;
            maxiter = 2, blocksize = 2, cutoff = 1e-8, verbose = false)
        @test isapprox(e_proj_r, e_cached_r; atol = 1e-9, rtol = 0)
    end

    @testset "Sliced convenience overloads" begin
        rng = MersenneTwister(11)
        nj = 2
        ns = 4
        N = nj * ns
        layout = HFDMRG.SliceLayout(fill(nj, ns))
        H = randn(rng, N, N)
        H = (H + H') / 2
        V6 = 0.01 * randn(rng, nj, nj, nj, nj, ns, ns)
        backend = HFDMRG.SlicedBasisBackend(layout, V6)
        Nup = 2
        psiup0 = orthonormal_cols(rng, N, Nup)
        _, _, e_backend = solve_hfdmrg(H, backend, psiup0; maxiter = 2, blocksize = 2, cutoff = 1e-8)
        _, _, e_layout = solve_hfdmrg(H, layout, V6, psiup0; maxiter = 2, blocksize = 2, cutoff = 1e-8)
        @test isapprox(e_backend, e_layout; atol = 1e-10, rtol = 0)

        psiup0_uhf = orthonormal_cols(rng, N, 1)
        psidn0 = orthonormal_cols(rng, N, 1)
        _, _, e_uhf_backend = solve_hfdmrg(H, backend, psiup0_uhf, psidn0;
            maxiter = 2, blocksize = 2, cutoff = 1e-8)
        _, _, e_uhf_split_backend = solve_hfdmrg(H, H, backend, psiup0_uhf, psidn0;
            maxiter = 2, blocksize = 2, cutoff = 1e-8)
        _, _, e_uhf_split_layout = solve_hfdmrg(H, H, layout, V6, psiup0_uhf, psidn0;
            maxiter = 2, blocksize = 2, cutoff = 1e-8)
        @test e_uhf_split_backend == e_uhf_backend
        @test e_uhf_split_layout == e_uhf_backend
    end

    @testset "Sliced end-to-end sanity" begin
        rng = MersenneTwister(19)
        nj = 2
        ns = 6
        N = nj * ns
        layout = HFDMRG.SliceLayout(fill(nj, ns))
        H = randn(rng, N, N)
        H = (H + H') / 2
        V6 = 0.01 * randn(rng, nj, nj, nj, nj, ns, ns)
        Nup = 3
        psiup0 = orthonormal_cols(rng, N, Nup)
        rho0 = psiup0 * psiup0'
        F0 = copy(H)
        HFDMRG.sliced_add_fock_r!(F0, rho0, V6, nj, ns)
        energy0 = tr(rho0 * (F0 + H))
        _, _, energy = solve_hfdmrg(H, layout, V6, psiup0;
            maxiter = 2, blocksize = 2, cutoff = 1e-8, verbose = false)
        @test isfinite(energy)
        @test energy <= energy0 + 1e-6 * max(1.0, abs(energy0))
    end

    @testset "Cached sliced backend vs projection" begin
        rng = MersenneTwister(23)
        nj = 2
        ns = 6
        N = nj * ns
        layout = HFDMRG.SliceLayout(fill(nj, ns))
        H = randn(rng, N, N)
        H = (H + H') / 2
        V6 = 0.01 * randn(rng, nj, nj, nj, nj, ns, ns)
        backend_proj = HFDMRG.SlicedBasisBackend(layout, V6)
        backend_cached = HFDMRG.SlicedBasisBackendCached(layout, V6)
        Nup = 1
        Ndn = 1
        psiup0 = orthonormal_cols(rng, N, Nup)
        psidn0 = orthonormal_cols(rng, N, Ndn)
        _, _, e_proj = solve_hfdmrg(H, backend_proj, psiup0, psidn0;
            maxiter = 2, blocksize = 2, cutoff = 1e-8, verbose = false)
        _, _, e_cached = solve_hfdmrg(H, backend_cached, psiup0, psidn0;
            maxiter = 2, blocksize = 2, cutoff = 1e-8, verbose = false)
        @test isapprox(e_proj, e_cached; atol = 1e-9, rtol = 0)
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
            maxiter, blocksize, cutoff, verbose = false)
        _, _, e_legacy = Legacy.dohfdmrg(H, V, psiup0, psidn0, false, 1,
            blocksize, maxiter, cutoff; verbose = false)
        @test isapprox(e_new, e_legacy; atol = 1e-6, rtol = 0)

        Nup_r = 1
        psiup_r = orthonormal_cols(rng, N, Nup_r)
        _, _, e_new_r = solve_hfdmrg(H, V, psiup_r;
            maxiter, blocksize, cutoff, verbose = false)
        _, _, e_legacy_r = Legacy.dohfdmrg(H, V, psiup_r, psiup_r, true, 1,
            blocksize, maxiter, cutoff; verbose = false)
        @test isapprox(e_new_r, e_legacy_r; atol = 1e-6, rtol = 0)
    end
finally
    empty!(LOAD_PATH)
    append!(LOAD_PATH, old_load_path)
end
