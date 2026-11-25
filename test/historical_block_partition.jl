@testset "Historical automatic block partition" begin
    hopping(n) = SymTridiagonal(zeros(n), fill(-1.0, n - 1))
    band_energy(n, occupied) = sum(-2cospi(k / (n + 1)) for k in 1:occupied)

    for (n, occupied, blocksize, nblockcenter) in ((8, 4, 1, 1), (10, 2, 200, 3))
        H = Matrix(hopping(n))
        V = zeros(n, n)
        seed = eigen(hopping(n)).vectors[:, 1:occupied]
        up, dn, energy = solve_hfdmrg(H, V, seed;
            blocksize, nblockcenter, maxiter=2)
        @test energy ≈ 2band_energy(n, occupied) atol=2e-12
        @test energy ≈ dot(up, H * up) + dot(dn, H * dn) atol=2e-12
        @test up == dn
        @test up' * up ≈ I atol=2e-12
    end

    n = 8
    H = Matrix(hopping(n))
    V = zeros(n, n)
    vectors = eigen(hopping(n)).vectors
    seed_up, seed_dn = vectors[:, 1:4], vectors[:, 1:3]
    for split in (false, true)
        Hdn = split ? 0.75H : H
        args = split ? (H, Hdn, V, seed_up, seed_dn) : (H, V, seed_up, seed_dn)
        up, dn, energy = solve_hfdmrg(args...; blocksize=1, maxiter=2)
        @test energy ≈ band_energy(n, 4) +
            (split ? 0.75 : 1.0) * band_energy(n, 3) atol=2e-12
        @test energy ≈ dot(up, H * up) + dot(dn, Hdn * dn) atol=2e-12
        @test up' * up ≈ I atol=2e-12
        @test dn' * dn ≈ I atol=2e-12
    end

    @test HFDMRG.getblocksizes(12, 2, 2) ==
        (6, fill(2, 6), [1:2, 3:4, 5:6, 7:8, 9:10, 11:12])
    @test_throws ArgumentError solve_hfdmrg(H, V, seed_up; blocksize=0)
    @test_throws ArgumentError solve_hfdmrg(H, V, seed_up; nblockcenter=0)
    @test_throws ArgumentError solve_hfdmrg(H, H, V, seed_up, seed_dn; nblockcenter=5)
    @test_throws ArgumentError solve_hfdmrg(H[1:4, 1:4], V[1:4, 1:4],
        reshape([1.0, 0.0, 0.0, 0.0], 4, 1))
end
