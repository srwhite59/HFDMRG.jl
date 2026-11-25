@testset "Compact local interactions with zero intercell channels" begin
    onsite = 0.125
    for sites in (36, 72), spin in (:rhf, :uhf)
        diagonal = [-1.0 + 0.05i for i in 1:sites]
        one_body = BandedOneBody(reshape(diagonal, 1, sites))
        triangle = [i == j ? onsite : 0.0 for j in 1:18 for i in 1:j]
        results = map((0, 1)) do channels
            operator = UnitCellInteraction(sites, 1, zeros(channels),
                zeros(0, 0), zeros(18, channels), zeros(channels, 18), triangle)
            occupied = sites ÷ 36
            result = solve_hfdmrg(one_body, operator; spin, exact=true,
                state_maxdim=occupied, state_cutoff=0.0, energy_tolerance=1e-10,
                maximum_half_sweeps=4)
            lifecycle = result.state.lifecycle
            alpha, beta = zeros(sites, occupied), zeros(sites, occupied)
            if spin === :rhf
                workspace = HFDMRG.RHFLifecycleWorkspace(one_body, operator, occupied)
                HFDMRG.reconstruct_terminal!(alpha, lifecycle, workspace)
                copyto!(beta, alpha)
            else
                workspace = HFDMRG.UHFLifecycleWorkspace(one_body, operator,
                    occupied, occupied)
                HFDMRG.reconstruct_uhf_terminal!(alpha, beta, lifecycle, workspace)
            end
            da, db = vec(sum(abs2, alpha; dims=2)), vec(sum(abs2, beta; dims=2))
            # A purely onsite interaction contributes only opposite-spin pairs.
            reference = dot(diagonal, da + db) + onsite * dot(da, db)
            @test abs(reference) > 0.1
            @test result.represented_energy ≈ reference atol=2e-12 rtol=0
            @test alpha' * alpha ≈ I atol=2e-12
            @test beta' * beta ≈ I atol=2e-12
            (; energy=result.represented_energy, pa=alpha * alpha', pb=beta * beta')
        end
        @test results[1].energy ≈ results[2].energy atol=2e-12 rtol=0
        @test results[1].pa ≈ results[2].pa atol=2e-12 rtol=0
        @test results[1].pb ≈ results[2].pb atol=2e-12 rtol=0
    end
end
