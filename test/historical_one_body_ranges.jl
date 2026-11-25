@testset "Historical structural one-body ranges" begin
    up = zeros(5, 5)
    dn = zeros(5, 5)
    up[1, 3] = up[3, 1] = -1e-14
    dn[3, 5] = dn[5, 3] = -2e-14
    # All-zero rows carry the empty-range convention (last=1, first=N).
    @test HFDMRG.getfinal(up, 5) == ([3, 1, 1, 1, 1], [3, 5, 1, 5, 5])
    @test HFDMRG.getfinal(dn, 5) == ([1, 1, 5, 1, 3], [5, 5, 5, 5, 3])
    @test HFDMRG.getfinal(up, dn, 5) == ([3, 1, 5, 1, 3], [3, 5, 1, 5, 3])

    n = 12
    chain = diagm(1 => fill(-1.0, n-1), -1 => fill(-1.0, n-1))
    V = zeros(n, n)
    C = Matrix(qr(randn(MersenneTwister(7), n, 2)).Q)
    beta = C[:, 1:1]
    levels = [-2cospi(k/(n+1)) for k = 1:2]
    expected = 2sum(levels)
    for t in (1.0, 1e-14)
        H = t*chain
        kw = (; blocksize=2, maxiter=10, cutoff=1e-12*t, scf_cutoff=1e-12*t)
        @test last(solve_hfdmrg(H, V, C; kw...))/t ≈ expected atol=1e-10 rtol=0
        @test last(solve_hfdmrg(H, V, C, C; kw...))/t ≈ expected atol=1e-10 rtol=0
        split = solve_hfdmrg(H, 2H, V, C, beta; kw...)
        @test split[3]/t ≈ sum(levels)+2levels[1] atol=1e-10 rtol=0
        @test split[3]/t ≈ dot(split[1], chain*split[1]) +
            2dot(split[2], chain*split[2]) atol=1e-10 rtol=0
    end
end
