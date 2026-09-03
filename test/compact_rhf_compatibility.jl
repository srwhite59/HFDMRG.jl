@testset "dual RHF/UHF backend dispatch boundaries" begin
    sites = 12
    h = Matrix(Diagonal(vcat(-2.0, collect(1.0:sites-1))))
    v = zeros(sites, sites)
    up = reshape(vcat(1.0, zeros(sites-1)), sites, 1)
    dn = reshape(vcat(0.0, 1.0, zeros(sites-2)), sites, 1)
    layout = HFDMRG.SliceLayout(fill(1, sites))
    sliced_data = zeros(1, 1, 1, 1, sites, sites)
    sliced = HFDMRG.SlicedBasisBackend(layout, sliced_data)
    cached = HFDMRG.SlicedBasisBackendCached(layout, sliced_data)
    cached_coulomb = HFDMRG._SlicedBasisBackendCachedCoulomb(layout,
        sliced_data)
    target = HFDMRG.DensityDensityTargetResidualBackend(v, up,
        reshape([0.01], 1, 1))
    solver = (; maxiter=1, blocksize=2, cutoff=0.0,
        scf_cutoff=Inf, verbose=false)
    sliced_solver = (; solver..., block_partition=layout)

    restricted = (
        () -> HFDMRG.solve_hfdmrg(h, v, up; solver...),
        () -> HFDMRG.solve_hfdmrg(h, sliced, up; sliced_solver...),
        () -> HFDMRG.solve_hfdmrg(h, layout, sliced_data, up; sliced_solver...),
        () -> HFDMRG.solve_hfdmrg(h, cached, up; sliced_solver...),
        () -> HFDMRG.solve_hfdmrg(h, cached_coulomb, up; sliced_solver...),
        () -> HFDMRG.solve_hfdmrg(h, target, up; solver...),
    )
    for run in restricted
        result = run()
        @test result isa Tuple{Matrix{Float64},Matrix{Float64},Float64}
        @test result[1] == result[2]
        @test size(result[1]) == (sites, 1)
        @test isfinite(result[3])
    end

    unrestricted = (
        () -> HFDMRG.solve_hfdmrg(h, v, up, dn; solver...),
        () -> HFDMRG.solve_hfdmrg(h, sliced, up, dn; sliced_solver...),
        () -> HFDMRG.solve_hfdmrg(h, layout, sliced_data, up, dn;
            sliced_solver...),
        () -> HFDMRG.solve_hfdmrg(h, cached, up, dn; sliced_solver...),
        () -> HFDMRG.solve_hfdmrg(h, cached_coulomb, up, dn;
            sliced_solver...),
        () -> HFDMRG.solve_hfdmrg(h, target, up, dn; solver...),
    )
    for run in unrestricted
        result = run()
        @test result isa Tuple{Matrix{Float64},Matrix{Float64},Float64}
        @test size(result[1]) == size(result[2]) == (sites, 1)
        @test isfinite(result[3])
    end

    one_body, operator = nucleus_fixture(2)
    compact = HFDMRG.solve_hfdmrg(one_body, operator;
        state_maxdim=4, state_cutoff=1e-12,
        energy_tolerance=1e-9, maximum_half_sweeps=4, exact=true)
    @test compact isa HFDMRG.HFDMRGResult
    @test compact.state isa HFDMRG.RHFCompactState
    @test !(compact isa Tuple)
    compact_uhf = HFDMRG.solve_hfdmrg(one_body, operator; spin=:uhf,
        initialization=:atomic_neel, state_maxdim=4, state_cutoff=1e-12,
        energy_tolerance=1e-9, maximum_half_sweeps=4, exact=true)
    @test compact_uhf isa HFDMRG.UHFDMRGResult
    @test compact_uhf.state isa HFDMRG.UHFCompactState
    @test !(compact_uhf isa Tuple)
    @test which(HFDMRG.solve_hfdmrg,
        (typeof(one_body), typeof(operator))) !==
        which(HFDMRG.solve_hfdmrg, (typeof(h), typeof(v), typeof(up)))
    @test which(HFDMRG.solve_hfdmrg,
        (typeof(one_body), typeof(operator))) !==
        which(HFDMRG.solve_hfdmrg,
            (typeof(h), typeof(v), typeof(up), typeof(dn)))
    @test isempty(Test.detect_ambiguities(HFDMRG; recursive=true))
end
