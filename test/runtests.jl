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

    mutable struct FockCountingBackend{B}
        backend::B
        calls::Int
        windows::Int
    end

    struct CountedFockWindow{W, B}
        window::W
        backend::B
    end

    HFDMRG.vee_init_block(side, ra, raV, phi, b::FockCountingBackend) =
        HFDMRG.vee_init_block(side, ra, raV, phi, b.backend)
    HFDMRG.vee_absorb_block(side, old, cra, Phi_old, Phi_C, phi, raV,
        b::FockCountingBackend) = HFDMRG.vee_absorb_block(side, old, cra,
        Phi_old, Phi_C, phi, raV, b.backend)

    function HFDMRG.vee_window(L, R, Cra, b::FockCountingBackend)
        b.windows += 1
        CountedFockWindow(HFDMRG.vee_window(L, R, Cra, b.backend), b)
    end

    function HFDMRG.vee_add_fock_r!(F, rho, win::CountedFockWindow)
        win.backend.calls += 1
        HFDMRG.vee_add_fock_r!(F, rho, win.window)
    end

    function HFDMRG.vee_add_fock!(Fup, Fdn, rhoup, rhodn, win::CountedFockWindow)
        win.backend.calls += 1
        HFDMRG.vee_add_fock!(Fup, Fdn, rhoup, rhodn, win.window)
    end

    @testset "Public API" begin
        @test isdefined(HFDMRG, :solve_hfdmrg)
        @test :solve_hfdmrg in names(HFDMRG, all = false)
    end

    @testset "Local Fock reuse" begin
        rng = MersenneTwister(719)
        layout = HFDMRG.SliceLayout(fill(2, 6))
        N = layout.offs[end]
        A = randn(rng, N, N); H = Matrix(Symmetric(A))
        B = randn(rng, N, N); V = 0.02 * Matrix(Symmetric(B))
        V6 = 0.002 * randn(rng, 2, 2, 2, 2, 6, 6)
        up = orthonormal_cols(rng, N, 1)
        dn = orthonormal_cols(rng, N, 1)
        Hdn = H + Diagonal(range(-0.05, 0.05; length = N))
        kw = (; maxiter = 1, blocksize = 2, cutoff = 0.0, verbose = false)
        density = HFDMRG.DensityDensityBackend(V)
        cached = HFDMRG.SlicedBasisBackendCached(layout, V6)
        routes = (
            (density, b -> HFDMRG.solve_hfdmrg_core(H, b, up, up;
                restricted = true, scf_cutoff = 0.0, kw...)),
            (density, b -> HFDMRG.solve_hfdmrg_core(H, b, up, dn;
                restricted = false, scf_cutoff = 0.0, kw...)),
            (density, b -> HFDMRG.solve_hfdmrg_core_split(H, Hdn, b, up, dn;
                scf_cutoff = 0.0, kw...)),
            (cached, b -> HFDMRG.solve_hfdmrg_core(H, b, up, dn;
                restricted = false, block_partition = layout,
                scf_cutoff = 0.0, kw...)),
        )
        for (backend, run) in routes
            counted = FockCountingBackend(backend, 0, 0)
            @test run(counted) == run(backend)
            @test counted.windows > 0
            @test counted.calls == 5 * counted.windows
        end

        counted = FockCountingBackend(density, 0, 0)
        HFDMRG.solve_hfdmrg_core(H, counted, up, up; restricted = true,
            scf_cutoff = Inf, kw...)
        @test counted.calls == 2 * counted.windows
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

    @testset "Per-sweep observer" begin
        rng = MersenneTwister(809)
        N = 12
        A = randn(rng, N, N)
        H = (A + A') / 2
        B = randn(rng, N, N)
        V = 0.01 * (B + B') / 2
        up0 = orthonormal_cols(rng, N, 1)
        dn0 = orthonormal_cols(rng, N, 1)
        kw = (; maxiter = 2, blocksize = 2, cutoff = 0.0,
            scf_cutoff = 1e-8, verbose = false)
        function full_energy(Hup, Hdn, V, info)
            rhoup = info.psiup * info.psiup'
            rhodn = info.psidn * info.psidn'
            occupation = diag(rhoup) + diag(rhodn)
            sum(Hup .* rhoup) + sum(Hdn .* rhodn) +
            0.5 * dot(occupation, V * occupation) -
            0.5 * sum(V .* (rhoup .* rhoup + rhodn .* rhodn))
        end

        baseline = solve_hfdmrg(H, V, up0, dn0; kw...)
        seen = Any[]
        observed = solve_hfdmrg(H, V, up0, dn0; kw...,
            observer = info -> (push!(seen, info); false))
        @test observed == baseline
        @test getproperty.(seen, :sweep) == [1, 2]
        @test all(info -> !info.converged, seen)
        @test size(seen[end].psiup) == size(up0)
        @test size(seen[end].psidn) == size(dn0)
        @test seen[end].energy == observed[3]
        @test seen[end].psiup === observed[1]
        @test seen[end].psidn === observed[2]
        @test all(info -> isapprox(info.energy, full_energy(H, H, V, info);
            atol = 1e-12, rtol = 0), seen)
        @test_throws ErrorException solve_hfdmrg(H, V, up0, dn0; kw...,
            observer = _ -> :invalid)
        @test_throws ErrorException solve_hfdmrg(H, V, up0, dn0; kw...,
            observer = _ -> error("observer failure"))

        stopped_seen = Any[]
        stopped = solve_hfdmrg(H, V, up0, dn0;
            maxiter = 5, blocksize = 2, cutoff = 0.0,
            scf_cutoff = 1e-8, verbose = false,
            observer = info -> (push!(stopped_seen, info); info.sweep == 2))
        @test getproperty.(stopped_seen, :sweep) == [1, 2]
        @test stopped == baseline

        split_seen = Any[]
        Hup = H + Diagonal(range(0, 0.1; length = N))
        Hdn = H - Diagonal(range(0, 0.1; length = N))
        split = solve_hfdmrg(Hup, Hdn, V, up0, dn0;
            maxiter = 1, blocksize = 2, cutoff = 0.0, scf_cutoff = 1e-8,
            verbose = false,
            observer = info -> (push!(split_seen, info); nothing))
        @test length(split_seen) == 1
        @test split_seen[1].energy == split[3]
        @test split_seen[1].psiup === split[1]
        @test split_seen[1].psidn === split[2]
        @test isapprox(split_seen[1].energy,
            full_energy(Hup, Hdn, V, split_seen[1]); atol = 1e-12, rtol = 0)

        Q = orthonormal_cols(rng, N, 2)
        zero_backend = HFDMRG.DensityDensityTargetResidualBackend(
            V, Q, zeros(3, 3))
        restricted_seen = Any[]
        restricted = solve_hfdmrg(H, zero_backend, up0;
            maxiter = 3, blocksize = 2, cutoff = Inf, verbose = false,
            observer = info -> (push!(restricted_seen, info); false))
        @test length(restricted_seen) == 1
        @test restricted_seen[1].converged
        @test restricted_seen[1].psiup === restricted_seen[1].psidn
        @test restricted_seen[1].psiup === restricted[1]
        @test restricted[1] === restricted[2]
    end

    @testset "Density-density target residual" begin
        rng = MersenneTwister(411)
        N, m = 10, 3
        P = m * (m + 1) ÷ 2
        Q = orthonormal_cols(rng, N, m)
        A = randn(rng, P, P)
        Rpair = 0.01 * (A + A') / 2
        backend = HFDMRG.DensityDensityTargetResidualBackend(zeros(N, N), Q, Rpair)
        pair(p, q) = max(p, q) * (max(p, q) - 1) ÷ 2 + min(p, q)
        function explicit_jk(D)
            J, K = zeros(m, m), zeros(m, m)
            for p = 1:m, q = 1:m, r = 1:m, s = 1:m
                J[p, q] += Rpair[pair(p, q), pair(r, s)] * D[r, s]
                K[p, q] += Rpair[pair(p, r), pair(q, s)] * D[r, s]
            end
            J, K
        end

        Qbad = copy(Q)
        Qbad[:, 1] .*= 1.01
        @test_throws ErrorException HFDMRG.DensityDensityTargetResidualBackend(
            zeros(N, N), Qbad, Rpair)
        Rbad = copy(Rpair)
        Rbad[1, 2] += 1e-4
        @test_throws ErrorException HFDMRG.DensityDensityTargetResidualBackend(
            zeros(N, N), Q, Rbad)
        @test_throws ErrorException HFDMRG.DensityDensityTargetResidualBackend(
            zeros(Int, N, N), Q, Rpair)

        Lphi = orthonormal_cols(rng, 2, 1)
        Rphi = orthonormal_cols(rng, 2, 1)
        L = HFDMRG.vee_init_block(:left, 1:2, 3:N, Lphi, backend)
        R = HFDMRG.vee_init_block(:right, 9:10, 1:8, Rphi, backend)
        OL = orthonormal_cols(rng, 3, 2)
        PhiL, PhiLC = OL[1:1, :], OL[2:3, :]
        Lphi_new = vcat(Lphi * PhiL, PhiLC)
        L = HFDMRG.vee_absorb_block(:left, L, 3:4, PhiL, PhiLC, Lphi_new,
            5:N, backend)
        OR = orthonormal_cols(rng, 3, 2)
        PhiRC, PhiR = OR[1:2, :], OR[3:3, :]
        Rphi_new = vcat(PhiRC, Rphi * PhiR)
        R = HFDMRG.vee_absorb_block(:right, R, 7:8, PhiR, PhiRC, Rphi_new,
            1:6, backend)
        win = HFDMRG.vee_window(L, R, 5:6, backend)
        Qwindow = vcat(Lphi_new' * Q[1:4, :], Q[5:6, :], Rphi_new' * Q[7:10, :])
        w = size(Qwindow, 1)

        rho = randn(rng, w, w)
        rho = (rho + rho') / 2
        F = zeros(w, w)
        HFDMRG.vee_add_fock_r!(F, rho, win)
        D = Qwindow' * rho * Qwindow
        J, K = explicit_jk(D)
        Fref = Qwindow * (2J - K) * Qwindow'
        @test isapprox(F, Fref; atol = 2e-12, rtol = 0)
        @test isapprox(tr(rho * F), 2 * sum(D .* J) - sum(D .* K);
            atol = 2e-12, rtol = 0)

        rhoup = randn(rng, w, w); rhoup = (rhoup + rhoup') / 2
        rhodn = randn(rng, w, w); rhodn = (rhodn + rhodn') / 2
        Fup, Fdn = zeros(w, w), zeros(w, w)
        HFDMRG.vee_add_fock!(Fup, Fdn, rhoup, rhodn, win)
        Da, Db = Qwindow' * rhoup * Qwindow, Qwindow' * rhodn * Qwindow
        Jtot = first(explicit_jk(Da + Db))
        Ka = last(explicit_jk(Da))
        Kb = last(explicit_jk(Db))
        @test isapprox(Fup, Qwindow * (Jtot - Ka) * Qwindow'; atol = 2e-12, rtol = 0)
        @test isapprox(Fdn, Qwindow * (Jtot - Kb) * Qwindow'; atol = 2e-12, rtol = 0)
        Eref = 0.5 * sum((Da + Db) .* Jtot) - 0.5 * sum(Da .* Ka) -
               0.5 * sum(Db .* Kb)
        @test isapprox(0.5 * sum(rhoup .* Fup) + 0.5 * sum(rhodn .* Fdn), Eref;
            atol = 2e-12, rtol = 0)

        H = randn(rng, N, N); H = (H + H') / 2
        V = 0.01 * randn(rng, N, N); V = (V + V') / 2
        zero_backend = HFDMRG.DensityDensityTargetResidualBackend(V, Q, zeros(P, P))
        up, dn = orthonormal_cols(rng, N, 1), orthonormal_cols(rng, N, 1)
        kw = (; maxiter = 1, blocksize = 1, cutoff = 1e-9, verbose = false)
        @test solve_hfdmrg(H, zero_backend, up; kw...) == solve_hfdmrg(H, V, up; kw...)
        @test solve_hfdmrg(H, zero_backend, up, dn; kw...) ==
              solve_hfdmrg(H, V, up, dn; kw...)
        Hup, Hdn = H + Diagonal(range(0, 0.1; length = N)), H - 0.1I
        @test solve_hfdmrg(Hup, Hdn, zero_backend, up, dn; kw...) ==
              solve_hfdmrg(Hup, Hdn, V, up, dn; kw...)
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
            F2[p, q] += rho[r, s] * (2.0 * v6_lookup(V6, p, r, q, s, nj) -
                                     v6_lookup(V6, p, r, s, q, nj))
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
            vdir = v6_lookup(V6, p, r, q, s, nj)
            vex = v6_lookup(V6, p, r, s, q, nj)
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
            F2[p, q] += rho[r, s] * (2.0 * vblocks_lookup(layout, Vblocks, p, r, q, s) -
                                     vblocks_lookup(layout, Vblocks, p, r, s, q))
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
            vdir = vblocks_lookup(layout, Vblocks, p, r, q, s)
            vex = vblocks_lookup(layout, Vblocks, p, r, s, q)
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

    @testset "Cached sliced ordered sectors" begin
        rng = MersenneTwister(31)
        d, ns, N = 2, 3, 6
        layout = HFDMRG.SliceLayout(fill(d, ns))
        raw = [[randn(rng, d, d, d, d) for _ = 1:ns] for _ = 1:ns]
        @test minimum([norm(raw[n][m] - permutedims(raw[m][n], (3, 4, 1, 2)))
            for n = 1:ns, m = 1:ns]) > 1e-6
        Q = orthonormal_cols(rng, N, N)
        rho = Matrix(Symmetric(Q * Diagonal(collect(range(0.1, 0.8; length = N))) * Q'))
        alpha = 0.37
        rhodn = alpha * rho
        @test all(issymmetric(D) && norm(D * D - D) > 1e-2 for D in (rho, rhodn))
        relerr(A, B) = norm(A - B, Inf) / max(1.0, norm(B, Inf))
        function make_window(backend, Lra, Cra, Rra, Lphi, Rphi)
            Lvee = HFDMRG.vee_init_block(:left, Lra, (Lra[end] + 1):N, Lphi, backend)
            Rvee = HFDMRG.vee_init_block(:right, Rra, 1:(Rra[1] - 1), Rphi, backend)
            HFDMRG.vee_window(Lvee, Rvee, Cra, backend)
        end
        function focks(win)
            F, Fup, Fdn = zeros(N, N), zeros(N, N), zeros(N, N)
            HFDMRG.vee_add_fock_r!(F, rho, win)
            HFDMRG.vee_add_fock!(Fup, Fdn, rho, rhodn, win)
            F, Fup, Fdn
        end
        eye = Matrix{Float64}(I, d, d)
        Lra, Cra, Rra = 1:2, 3:4, 5:6
        sector_routes = withenv("HFDMRG_BENCH_TIMING" => "1") do
            HFDMRG._bench_timing_reset!()
            for n = 1:ns, m = 1:ns
                Vblocks = [[zeros(d, d, d, d) for _ = 1:ns] for _ = 1:ns]
                Vblocks[n][m] .= raw[n][m]
                J, K = zeros(N, N), zeros(N, N)
                for p = 1:N, q = 1:N, r = 1:N, s = 1:N
                    v = vblocks_lookup(layout, Vblocks, p, q, r, s)
                    J[p, r] += rho[q, s] * v
                    K[p, s] += rho[r, q] * v
                end
                @test min(norm(J, Inf), norm(K, Inf)) > 1e-8
                backends = (HFDMRG.SlicedBasisBackend(layout, Vblocks),
                    HFDMRG.SlicedBasisBackendCached(layout, Vblocks))
                outputs = map(b -> focks(make_window(b, Lra, Cra, Rra, eye, eye)), backends)
                for (F, Fup, Fdn) in outputs
                    Kgot = (Fdn - Fup) / (1 - alpha)
                    Jgot = (Fup + Kgot) / (1 + alpha)
                    @test relerr(F, 2J - K) <= 1e-12
                    @test relerr(Jgot, J) <= 1e-12
                    @test relerr(Kgot, K) <= 1e-12
                end
                @test maximum(relerr(outputs[2][i], outputs[1][i]) for i = 1:3) <= 1e-12
            end
            HFDMRG._bench_timing_snapshot()
        end
        @test (sector_routes[:window_local_calls],
            sector_routes[:window_projection_calls]) == (9.0, 0.0)
        aligned = make_window(HFDMRG.SlicedBasisBackendCached(layout, raw),
            Lra, Cra, Rra, eye, eye)
        Falloc, Fupalloc, Fdnalloc = focks(aligned)
        fill!(Falloc, 0); fill!(Fupalloc, 0); fill!(Fdnalloc, 0)
        @test (@allocated HFDMRG.vee_add_fock_r!(Falloc, rho, aligned)) == 0
        @test (@allocated HFDMRG.vee_add_fock!(Fupalloc, Fdnalloc,
            rho, rhodn, aligned)) == 0
        split_proj = make_window(HFDMRG.SlicedBasisBackend(layout, raw),
            1:1, 2:5, 6:6, ones(1, 1), ones(1, 1))
        split_cached, split_routes = withenv("HFDMRG_BENCH_TIMING" => "1") do
            HFDMRG._bench_timing_reset!()
            win = make_window(HFDMRG.SlicedBasisBackendCached(layout, raw),
                1:1, 2:5, 6:6, ones(1, 1), ones(1, 1))
            win, HFDMRG._bench_timing_snapshot()
        end
        split_outputs = focks(split_proj), focks(split_cached)
        @test maximum(relerr(split_outputs[2][i], split_outputs[1][i])
            for i = 1:3) <= 1e-12
        @test (split_routes[:window_local_calls],
            split_routes[:window_projection_calls]) == (0.0, 1.0)
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

    @testset "Cached sliced incremental absorption" begin
        rng = MersenneTwister(92)
        relerr(A, B) = norm(A - B, Inf) / max(1.0, norm(B, Inf))

        function grow(state, side, cra, knew, backend)
            mold, dc = size(state.phi, 2), length(cra)
            O = orthonormal_cols(rng, mold + dc, knew)
            rotation = orthonormal_cols(rng, knew, knew)
            if side === :left
                Phi_old, Phi_C = O[1:mold, :], O[mold + 1:end, :]
                phi = vcat(state.phi * Phi_old, Phi_C) * rotation
                raV = (last(cra) + 1):backend.layout.offs[end]
            else
                Phi_C, Phi_old = O[1:dc, :], O[dc + 1:end, :]
                phi = vcat(Phi_C, state.phi * Phi_old) * rotation
                raV = 1:(first(cra) - 1)
            end
            HFDMRG.vee_absorb_block(side, state, cra, Phi_old, Phi_C, phi, raV,
                backend)
        end

        function fock_error(wins, w)
            rho = randn(rng, w, w); rho = (rho + rho') / 2
            rhoup = randn(rng, w, w); rhoup = (rhoup + rhoup') / 2
            rhodn = randn(rng, w, w); rhodn = (rhodn + rhodn') / 2
            Fr = [zeros(w, w) for _ = 1:3]
            Fu = [zeros(w, w) for _ = 1:3]
            Fd = [zeros(w, w) for _ = 1:3]
            for i = 1:3
                HFDMRG.vee_add_fock_r!(Fr[i], rho, wins[i])
                HFDMRG.vee_add_fock!(Fu[i], Fd[i], rhoup, rhodn, wins[i])
            end
            maximum((relerr(Fr[1], Fr[2]), relerr(Fr[1], Fr[3]),
                relerr(Fu[1], Fu[2]), relerr(Fu[1], Fu[3]),
                relerr(Fd[1], Fd[2]), relerr(Fd[1], Fd[3])))
        end

        withenv("HFDMRG_BENCH_TIMING" => "1") do
            HFDMRG._bench_timing_reset!()
            for dims in ([2, 2, 2, 2, 2, 2, 2], [2, 1, 3, 2, 3, 1, 2])
                layout = HFDMRG.SliceLayout(dims)
                ns, N = length(dims), layout.offs[end]
                V = if all(==(dims[1]), dims)
                    d = dims[1]
                    randn(rng, d, d, d, d, ns, ns)
                else
                    [[randn(rng, dims[n], dims[n], dims[m], dims[m]) for m = 1:ns]
                     for n = 1:ns]
                end
                cached = HFDMRG.SlicedBasisBackendCached(layout, V)
                projection = HFDMRG.SlicedBasisBackend(layout, V)
                L = HFDMRG.vee_init_block(:left, HFDMRG.orb_range(layout, 1),
                    (dims[1] + 1):N, orthonormal_cols(rng, dims[1], 2), cached)
                R = HFDMRG.vee_init_block(:right, HFDMRG.orb_range(layout, 7),
                    1:(N - dims[7]), orthonormal_cols(rng, dims[7], 2), cached)
                L = grow(L, :left, HFDMRG.orb_range(layout, 2), 3, cached)
                L = grow(L, :left, HFDMRG.orb_range(layout, 3), 2, cached)
                R = grow(R, :right, HFDMRG.orb_range(layout, 6), 3, cached)
                R = grow(R, :right, HFDMRG.orb_range(layout, 5), 2, cached)
                Ldirect = HFDMRG.vee_init_block(:left, L.ra, (last(L.ra) + 1):N,
                    L.phi, cached)
                Rdirect = HFDMRG.vee_init_block(:right, R.ra, 1:(first(R.ra) - 1),
                    R.phi, cached)
                Lproj = HFDMRG.vee_init_block(:left, L.ra, (last(L.ra) + 1):N,
                    L.phi, projection)
                Rproj = HFDMRG.vee_init_block(:right, R.ra, 1:(first(R.ra) - 1),
                    R.phi, projection)
                Cra = HFDMRG.orb_range(layout, 4)
                wins = (HFDMRG.vee_window(L, R, Cra, cached),
                    HFDMRG.vee_window(Ldirect, Rdirect, Cra, cached),
                    HFDMRG.vee_window(Lproj, Rproj, Cra, projection))
                @test fock_error(wins, 4 + length(Cra)) <= 1e-11
            end
            counts = HFDMRG._bench_timing_snapshot()
            @test (counts[:cache_init_calls], counts[:cache_incremental_calls],
                counts[:cache_fallback_calls]) == (8.0, 8.0, 0.0)

            HFDMRG._bench_timing_reset!()
            layout = HFDMRG.SliceLayout(fill(3, 5))
            N = layout.offs[end]
            V = randn(rng, 3, 3, 3, 3, 5, 5)
            cached = HFDMRG.SlicedBasisBackendCached(layout, V)
            projection = HFDMRG.SlicedBasisBackend(layout, V)
            L0 = HFDMRG.vee_init_block(:left, 1:1, 2:N, ones(1, 1), cached)
            R0 = HFDMRG.vee_init_block(:right, N:N, 1:(N - 1), ones(1, 1), cached)
            L0proj = HFDMRG.vee_init_block(:left, 1:1, 2:N, ones(1, 1), projection)
            R0proj = HFDMRG.vee_init_block(:right, N:N, 1:(N - 1), ones(1, 1), projection)
            split = HFDMRG.vee_window(L0, R0, 2:(N - 1), cached)
            splitproj = HFDMRG.vee_window(L0proj, R0proj, 2:(N - 1), projection)
            @test fock_error((split, split, splitproj), N) <= 1e-11
            Lphi, Rphi = orthonormal_cols(rng, 3, 2), orthonormal_cols(rng, 3, 2)
            L = HFDMRG.vee_absorb_block(:left, L0, 2:3, zeros(1, 2), zeros(2, 2),
                Lphi, 4:N, cached)
            R = HFDMRG.vee_absorb_block(:right, R0, 13:14, zeros(1, 2), zeros(2, 2),
                Rphi, 1:12, cached)
            Ldirect = HFDMRG.vee_init_block(:left, 1:3, 4:N, Lphi, cached)
            Rdirect = HFDMRG.vee_init_block(:right, 13:15, 1:12, Rphi, cached)
            Lproj = HFDMRG.vee_init_block(:left, 1:3, 4:N, Lphi, projection)
            Rproj = HFDMRG.vee_init_block(:right, 13:15, 1:12, Rphi, projection)
            partial_wins = (HFDMRG.vee_window(L, R, 4:12, cached),
                HFDMRG.vee_window(Ldirect, Rdirect, 4:12, cached),
                HFDMRG.vee_window(Lproj, Rproj, 4:12, projection))
            @test fock_error(partial_wins, 13) <= 1e-11

            oldphi = orthonormal_cols(rng, 3, 1)
            old = HFDMRG.vee_init_block(:left, 1:3, 4:N, oldphi, cached)
            offphi = orthonormal_cols(rng, 6, 2)
            off = HFDMRG.vee_absorb_block(:left, old, 4:6, zeros(1, 2), zeros(3, 2),
                offphi, 7:N, cached)
            offdirect = HFDMRG.vee_init_block(:left, 1:6, 7:N, offphi, cached)
            offproj = HFDMRG.vee_init_block(:left, 1:6, 7:N, offphi, projection)
            off_wins = (HFDMRG.vee_window(off, R, 7:12, cached),
                HFDMRG.vee_window(offdirect, Rdirect, 7:12, cached),
                HFDMRG.vee_window(offproj, Rproj, 7:12, projection))
            @test fock_error(off_wins, 10) <= 1e-11

            offRphi = orthonormal_cols(rng, 6, 2)
            offR = HFDMRG.vee_absorb_block(:right, Rdirect, 10:12, zeros(2, 2),
                zeros(3, 2), offRphi, 1:9, cached)
            offRdirect = HFDMRG.vee_init_block(:right, 10:15, 1:9, offRphi, cached)
            offRproj = HFDMRG.vee_init_block(:right, 10:15, 1:9, offRphi, projection)
            offR_wins = (HFDMRG.vee_window(L, offR, 4:9, cached),
                HFDMRG.vee_window(Ldirect, offRdirect, 4:9, cached),
                HFDMRG.vee_window(Lproj, offRproj, 4:9, projection))
            @test fock_error(offR_wins, 10) <= 1e-11

            internalL = HFDMRG.vee_init_block(:left, 4:6, 7:N,
                orthonormal_cols(rng, 3, 1), cached)
            internalR = HFDMRG.vee_init_block(:right, 10:12, 1:9,
                orthonormal_cols(rng, 3, 1), cached)
            grow(internalL, :left, 7:9, 2, cached)
            grow(internalR, :right, 7:9, 2, cached)
            counts = HFDMRG._bench_timing_snapshot()
            @test (counts[:cache_incremental_calls], counts[:cache_fallback_calls]) ==
                (0.0, 6.0)
        end
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

    @testset "Cached sliced backend vs projection" begin
        nj = 2
        ns = 6
        N = nj * ns
        layout = HFDMRG.SliceLayout(fill(nj, ns))
        Q = [sqrt(2 / (N + 1)) * sinpi(i * j / (N + 1)) for i = 1:N, j = 1:N]
        H = Q * Diagonal(0.0:N - 1) * Q'
        V6 = zeros(nj, nj, nj, nj, ns, ns)
        for n = 1:ns, m = 1:ns, a = 1:nj, c = 1:nj
            p, q = (n - 1) * nj + a, (m - 1) * nj + c
            V6[a, a, c, c, n, m] = 0.005 / (1 + abs(p - q))
        end
        backend_proj = HFDMRG.SlicedBasisBackend(layout, V6)
        backend_cached = HFDMRG.SlicedBasisBackendCached(layout, V6)
        psi0 = Q[:, 1:2]
        psiup0, psidn0 = Q[:, 1:1], Q[:, 3:3]
        projector_error(A, B) = norm(A * A' - B * B')
        energy_error(a, b) = abs(a - b) / max(1.0, abs(a), abs(b))
        _, _, e_proj = solve_hfdmrg(H, backend_proj, psi0;
            maxiter = 2, blocksize = 2, cutoff = 1e-8, verbose = false)
        _, _, e_cached = solve_hfdmrg(H, backend_cached, psi0;
            maxiter = 2, blocksize = 2, cutoff = 1e-8, verbose = false)
        @test energy_error(e_proj, e_cached) <= 1e-10

        aligned = (; maxiter = 1, blocksize = 0, block_partition = layout,
            cutoff = 0.0, scf_cutoff = 0.0, verbose = false)
        result_proj = solve_hfdmrg(H, backend_proj, psi0; aligned...)
        result_cached = solve_hfdmrg(H, backend_cached, psi0; aligned...)
        @test energy_error(result_proj[3], result_cached[3]) <= 1e-10
        @test max(projector_error(result_proj[1], result_cached[1]),
            projector_error(result_proj[2], result_cached[2])) <= 1e-9
        route_counts = withenv("HFDMRG_BENCH_TIMING" => "1") do
            HFDMRG._bench_timing_reset!()
            solve_hfdmrg(H, backend_cached, psi0; aligned...)
            HFDMRG._bench_timing_snapshot()
        end
        @test (route_counts[:cache_init_calls], route_counts[:cache_incremental_calls],
            route_counts[:cache_fallback_calls]) == (2.0, 9.0, 0.0)
        @test (route_counts[:window_local_calls],
            route_counts[:window_projection_calls]) == (6.0, 0.0)
        @test result_proj == solve_hfdmrg(H, backend_proj, psi0;
            aligned..., blocksize = 999)

        Hup = H
        Hdn = Q * Diagonal([2.0, 1.0, 0.0, collect(3.0:N - 1)...]) * Q'
        split_proj = solve_hfdmrg(Hup, Hdn, backend_proj, psiup0, psidn0; aligned...)
        split_cached = solve_hfdmrg(Hup, Hdn, backend_cached, psiup0, psidn0; aligned...)
        @test energy_error(split_proj[3], split_cached[3]) <= 1e-10
        @test max(projector_error(split_proj[1], split_cached[1]),
            projector_error(split_proj[2], split_cached[2])) <= 1e-9

        legacy = solve_hfdmrg(H, backend_proj, psiup0, psidn0;
            maxiter = 1, blocksize = 2, cutoff = 1e-8, verbose = false)
        @test legacy == solve_hfdmrg(H, backend_proj, psiup0, psidn0;
            maxiter = 1, blocksize = 2, block_partition = nothing,
            cutoff = 1e-8, verbose = false)

        Vdensity = [0.01 / (1 + abs(i - j)) for i = 1:N, j = 1:N]
        density_legacy = solve_hfdmrg(H, Vdensity, psiup0, psidn0;
            maxiter = 1, blocksize = 2, cutoff = 0.0, verbose = false)
        @test density_legacy == solve_hfdmrg(H, Vdensity, psiup0, psidn0;
            maxiter = 1, blocksize = 999, block_partition = layout,
            cutoff = 0.0, verbose = false)

        be_layout = HFDMRG.SliceLayout(fill(9, 52))
        nb, bs, ranges = HFDMRG.getblocksizes(468, 4, 0, 1, be_layout, nothing)
        @test nb == 52 && bs == fill(9, 52) && ranges == [9b - 8:9b for b = 1:52]
        @test 2nb - 2 - 4 == 98

        @test_throws ErrorException solve_hfdmrg(H, backend_proj, psiup0, psidn0;
            aligned..., nblockcenter = 0)
        bad_layout = HFDMRG.SliceLayout(fill(2, 6)); bad_layout.offs[2] += 1
        @test_throws ErrorException solve_hfdmrg(H, backend_proj, psiup0, psidn0;
            aligned..., block_partition = bad_layout)
        @test_throws ErrorException solve_hfdmrg(H, backend_proj, psiup0, psidn0;
            aligned..., block_partition = HFDMRG.SliceLayout(fill(2, 5)))
        @test_throws ErrorException solve_hfdmrg(H, backend_proj, psiup0, psidn0;
            aligned..., nblockcenter = 3)
        mismatch = HFDMRG.SliceLayout([1, 3, 2, 2, 2, 2])
        @test_throws ErrorException solve_hfdmrg(H, backend_proj, psiup0, psidn0;
            aligned..., block_partition = mismatch)
    end

    @testset "Structured three-block center" begin
        N = 17
        layout = HFDMRG.SliceLayout([2, 3, 1, 5, 1, 3, 2])
        Q = [sqrt(2 / (N + 1)) * sinpi(i * j / (N + 1)) for i = 1:N, j = 1:N]
        V = [0.2 / (1 + abs(i - j)) for i = 1:N, j = 1:N]
        function density_energy(Hup, Hdn, up, dn)
            rhoup, rhodn = up * up', dn * dn'
            occupation = diag(rhoup) + diag(rhodn)
            sum(Hup .* rhoup) + sum(Hdn .* rhodn) +
            0.5 * dot(occupation, V * occupation) -
            0.5 * sum(V .* (rhoup .^ 2 + rhodn .^ 2))
        end
        kw = (; maxiter = 1, blocksize = 999, block_partition = layout,
            nblockcenter = 3, cutoff = 0.0, scf_cutoff = 1e-13, verbose = false)

        psi = Q[:, 1:1]
        rho = psi * psi'
        G = 2Diagonal(V * diag(rho)) - V .* rho
        H = Matrix(Symmetric(Q * Diagonal(2.0 .* (0:N - 1)) * Q' - G))
        result = solve_hfdmrg(H, V, psi; kw...)
        @test norm(result[1] * result[1]' - rho) <= 1e-12
        @test isapprox(result[3], density_energy(H, H, result[1], result[2]);
            atol = 1e-12, rtol = 0)

        up, dn = Q[:, 1:1], Q[:, 2:2]
        rhoup, rhodn = up * up', dn * dn'
        occupation = diag(rhoup) + diag(rhodn)
        Hup = Matrix(Symmetric(Q * Diagonal(2.0 .* (0:N - 1)) * Q' -
            Diagonal(V * occupation) + V .* rhoup))
        spectrum_dn = [1, 0, collect(2:N - 1)...]
        Hdn = Matrix(Symmetric(Q * Diagonal(2.0 .* spectrum_dn) * Q' -
            Diagonal(V * occupation) + V .* rhodn))
        result = solve_hfdmrg(Hup, Hdn, V, up, dn; kw...)
        @test max(norm(result[1] * result[1]' - rhoup),
            norm(result[2] * result[2]' - rhodn)) <= 1e-12
        @test isapprox(result[3], density_energy(Hup, Hdn, result[1], result[2]);
            atol = 1e-12, rtol = 0)
    end

    @testset "Complete fixed boundary spans" begin
        N = 12
        eye = Matrix{Float64}(I, N, N)
        orbital(i) = eye[:, i:i]
        projector(C) = C * C'
        layout = HFDMRG.SliceLayout(fill(2, 6))
        V = zeros(N, N)
        V6 = zeros(2, 2, 2, 2, 6, 6)
        projection = HFDMRG.SlicedBasisBackend(layout, V6)
        cached = HFDMRG.SlicedBasisBackendCached(layout, V6)
        kw = (; maxiter = 1, blocksize = 2, cutoff = 0.0,
            scf_cutoff = 0.0, verbose = false)
        partitions = ((;), (; block_partition = layout))

        H = Matrix(Diagonal(fill(4.0, N)))
        H[2, 2] = H[11, 11] = -2.0
        H[2, 11] = H[11, 2] = -0.4
        rhf0 = hcat(orbital(1), orbital(12))
        Pref = projector(hcat(orbital(2), orbital(11)))
        q = (orbital(2) + orbital(11)) / sqrt(2)
        up0, dn0 = orbital(1), orbital(12)
        Hup, Hdn = copy(H), copy(H)
        Hup[11, 11] = Hdn[2, 2] = 2.0
        Hup[2, 11] = Hup[11, 2] = 0.0
        Hdn[2, 11] = Hdn[11, 2] = 0.0

        for part in partitions
            rhf = solve_hfdmrg(H, V, rhf0; kw..., part...)
            @test isapprox(rhf[3], -8.0; atol = 1e-12, rtol = 0)
            @test norm(projector(rhf[1]) - Pref) <= 1e-12
            uhf = solve_hfdmrg(H, V, up0, dn0; kw..., part...)
            @test isapprox(uhf[3], -4.8; atol = 1e-12, rtol = 0)
            @test max(norm(projector(uhf[1]) - projector(q)),
                norm(projector(uhf[2]) - projector(q))) <= 1e-12
            split = solve_hfdmrg(Hup, Hdn, V, up0, dn0; kw..., part...)
            @test isapprox(split[3], -4.0; atol = 1e-12, rtol = 0)
            @test max(norm(projector(split[1]) - projector(orbital(2))),
                norm(projector(split[2]) - projector(orbital(11)))) <= 1e-12

            for backend in (projection, cached)
                rhf = solve_hfdmrg(H, backend, rhf0; kw..., part...)
                @test isapprox(rhf[3], -8.0; atol = 1e-12, rtol = 0)
                @test norm(projector(rhf[1]) - Pref) <= 1e-12
                uhf = solve_hfdmrg(H, backend, up0, dn0; kw..., part...)
                @test isapprox(uhf[3], -4.8; atol = 1e-12, rtol = 0)
                @test max(norm(projector(uhf[1]) - projector(q)),
                    norm(projector(uhf[2]) - projector(q))) <= 1e-12
                split = solve_hfdmrg(
                    Hup, Hdn, backend, up0, dn0; kw..., part...)
                @test isapprox(split[3], -4.0; atol = 1e-12, rtol = 0)
                @test max(norm(projector(split[1]) - projector(orbital(2))),
                    norm(projector(split[2]) - projector(orbital(11)))) <= 1e-12
            end
        end

        routes = withenv("HFDMRG_BENCH_TIMING" => "1") do
            HFDMRG._bench_timing_reset!()
            solve_hfdmrg(Hup, Hdn, cached, up0, dn0;
                kw..., block_partition = layout)
            HFDMRG._bench_timing_snapshot()
        end
        @test (routes[:cache_init_calls], routes[:cache_incremental_calls],
            routes[:cache_fallback_calls]) == (2.0, 9.0, 0.0)
        @test (routes[:window_local_calls], routes[:window_projection_calls]) ==
              (6.0, 0.0)

        complete = hcat((orbital(1) + orbital(11)) / sqrt(2),
            (orbital(2) + orbital(12)) / sqrt(2))
        Hcomplete = Matrix(Diagonal(fill(4.0, N)))
        Hcomplete[1, 1] = Hcomplete[11, 11] = 0.0
        Hcomplete[1, 11] = Hcomplete[11, 1] = -3.0
        Hcomplete[2, 2] = Hcomplete[12, 12] = 0.0
        Hcomplete[2, 12] = Hcomplete[12, 2] = -2.0
        unchanged = solve_hfdmrg(Hcomplete, V, complete; kw...)
        @test isapprox(unchanged[3], -10.0; atol = 1e-12, rtol = 0)
        @test norm(projector(unchanged[1]) - projector(complete)) <= 1e-12

        dims = [3, 1, 2, 1, 4]
        ragged_layout = HFDMRG.SliceLayout(dims)
        Nragged = sum(dims)
        Qragged = [sqrt(2 / (Nragged + 1)) *
                   sinpi(i * j / (Nragged + 1))
                   for i = 1:Nragged, j = 1:Nragged]
        levels = [-4.0, -3.0, -2.0, -1.0, collect(1.0:7.0)...]
        Hragged = Qragged * Diagonal(levels) * Qragged'
        Vragged = [[zeros(dims[n], dims[n], dims[m], dims[m])
                    for m = 1:length(dims)] for n = 1:length(dims)]
        ragged_backend =
            HFDMRG.SlicedBasisBackendCached(ragged_layout, Vragged)
        ragged_eye = Matrix{Float64}(I, Nragged, Nragged)
        ragged_start = ragged_eye[:, [1, 4, 5, 6, 7, 8, 11]]
        ragged_kw = (; maxiter = 1, blocksize = 999,
            block_partition = ragged_layout, cutoff = 0.0, scf_cutoff = 0.0,
            verbose = false)
        ragged_result, ragged_routes =
            withenv("HFDMRG_BENCH_TIMING" => "1") do
                HFDMRG._bench_timing_reset!()
                result = solve_hfdmrg(
                    Hragged, ragged_backend, ragged_start; ragged_kw...)
                result, HFDMRG._bench_timing_snapshot()
            end
        ragged_oracle = Qragged[:, 1:7]
        @test isapprox(ragged_result[3], 2sum(levels[1:7]);
            atol = 1e-12, rtol = 0)
        @test norm(projector(ragged_result[1]) -
                   projector(ragged_oracle)) <= 1e-12
        @test (count(>(1e-10), svdvals(ragged_result[1][1:3, :])),
            count(>(1e-10), svdvals(ragged_result[1][8:11, :]))) == (3, 4)
        @test (ragged_routes[:cache_init_calls],
            ragged_routes[:cache_incremental_calls],
            ragged_routes[:cache_fallback_calls]) == (2.0, 6.0, 0.0)
        @test (ragged_routes[:window_local_calls],
            ragged_routes[:window_projection_calls]) == (4.0, 0.0)
    end

    @testset "Diagonal one-body coupling ranges" begin
        N = 10
        eye = Matrix{Float64}(I, N, N)
        orbital(i) = eye[:, i:i]
        projector(C) = C * C'
        levels = [0.0, -4.0, fill(5.0, 6)..., -3.0, 1.0]
        H = Matrix(Diagonal(levels))
        Hdn = Matrix(Diagonal(reverse(levels)))
        V = zeros(N, N)
        layout = HFDMRG.SliceLayout(fill(2, 5))
        rhf0 = hcat(orbital(1), orbital(10))
        up0, dn0 = orbital(1), orbital(10)
        kw = (; maxiter = 1, blocksize = 2, cutoff = 0.0,
            scf_cutoff = 0.0, verbose = false)
        for part in ((;), (; block_partition = layout))
            rhf = solve_hfdmrg(H, V, rhf0; kw..., part...)
            @test isapprox(rhf[3], -14.0; atol = 1e-12, rtol = 0)
            @test norm(projector(rhf[1]) -
                       projector(hcat(orbital(2), orbital(9)))) <= 1e-12
            uhf = solve_hfdmrg(H, V, up0, dn0; kw..., part...)
            @test isapprox(uhf[3], -8.0; atol = 1e-12, rtol = 0)
            @test max(norm(projector(uhf[1]) - projector(orbital(2))),
                norm(projector(uhf[2]) - projector(orbital(2)))) <= 1e-12
            split = solve_hfdmrg(H, Hdn, V, up0, dn0; kw..., part...)
            @test isapprox(split[3], -8.0; atol = 1e-12, rtol = 0)
            @test max(norm(projector(split[1]) - projector(orbital(2))),
                norm(projector(split[2]) - projector(orbital(9)))) <= 1e-12
        end
    end
finally
    empty!(LOAD_PATH)
    append!(LOAD_PATH, old_load_path)
end
