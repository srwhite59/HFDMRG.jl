@testset "Compact RHF uses two partial groups without padding" begin
    for (prefix, suffix) in ((1, 1), (14, 15))
        sites = prefix + suffix
        phase = 19 - prefix
        atoms = [1:prefix, prefix+1:sites]
        bands = [zeros(1, sites); transpose([fill(-1.0, sites-1); 0.0])]
        one_body = BandedOneBody(bands)
        operator = UnitCellInteraction(sites, phase, [0.0], zeros(0, 0),
            zeros(18, 1), zeros(1, 18), zeros(171))
        dense = SymTridiagonal(zeros(sites), fill(-1.0, sites-1))
        ground_energy = -4cospi(1 / (sites + 1))
        for exact in (true, false), starting_orientation in (:physical, :reflected)
            result = solve_hfdmrg(one_body, operator; starting_orientation,
                atom_intervals=atoms, exact, state_maxdim=2, state_cutoff=1e-12,
                energy_tolerance=1e-10, maximum_half_sweeps=4)
            lifecycle = result.state.lifecycle
            @test lifecycle.intervals == atoms
            @test lifecycle.root.rank == sites
            @test size(lifecycle.root.covariance) == (sites, sites)
            @test length(lifecycle.collection) == 1
            @test lifecycle.original_rows_read == 0
            @test result.converged
            @test result.represented_energy ≈ ground_energy atol=2e-12 rtol=0
            @test abs(result.replay_energy_change) <= 2e-12
            occupied = zeros(sites, 1)
            workspace = HFDMRG.RHFLifecycleWorkspace(one_body, operator, 2)
            HFDMRG.reconstruct_terminal!(occupied, lifecycle, workspace)
            @test occupied' * occupied ≈ I atol=2e-12
            @test result.represented_energy ≈ 2dot(occupied, dense * occupied) atol=2e-12
        end
    end
end
