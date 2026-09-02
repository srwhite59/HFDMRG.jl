@testset "superseded RHF compatibility boundaries" begin
    sites = 4
    h = Matrix(Diagonal(collect(1.0:sites)))
    v = zeros(sites, sites)
    occupied = reshape([1.0, 0.0, 0.0, 0.0], sites, 1)
    layout = HFDMRG.SliceLayout(fill(1, sites))
    sliced_data = zeros(1, 1, 1, 1, sites, sites)
    sliced = HFDMRG.SlicedBasisBackend(layout, sliced_data)
    cached = HFDMRG.SlicedBasisBackendCached(layout, sliced_data)
    cached_coulomb = HFDMRG._SlicedBasisBackendCachedCoulomb(layout,
        sliced_data)
    target = HFDMRG.DensityDensityTargetResidualBackend(v, occupied,
        zeros(1, 1))

    routes = (
        () -> HFDMRG.solve_hfdmrg(h, v, occupied),
        () -> HFDMRG.solve_hfdmrg(h, sliced, occupied),
        () -> HFDMRG.solve_hfdmrg(h, layout, sliced_data, occupied),
        () -> HFDMRG.solve_hfdmrg(h, cached, occupied),
        () -> HFDMRG.solve_hfdmrg(h, cached_coulomb, occupied),
        () -> HFDMRG.solve_hfdmrg(h, target, occupied),
    )
    for route in routes
        error = try
            route()
            nothing
        catch exception
            exception
        end
        @test error isa ArgumentError
        @test occursin("RHF solver route was removed", error.msg)
        @test occursin("BandedOneBody", error.msg)
    end
end
