@testset "Compact finite UHF one-spin discarded probability" begin
    cutoff = 0.015
    control = HFDMRG.RHFStateControl(1, cutoff)
    capped = HFDMRG.RHFStateControl(0, cutoff)
    restriction = [0.99;;]
    uhf = HFDMRG._factor_uhf_state_link(restriction, 0, 1, control)
    rhf = HFDMRG._factor_state_link(restriction, 0, 1, control)
    rhf_capped = HFDMRG._factor_state_link(restriction, 0, 1, capped)
    @test uhf.selection.completed == 1 && uhf.active_rank == 0
    @test uhf.selection.discarded_squared_weight ≈ 0.01 atol=2e-15
    @test uhf.selection.limiting_reason === :cutoff
    @test rhf.active_rank == 1 && iszero(rhf.selection.discarded_squared_weight)
    @test rhf_capped.selection.completed == 1
    @test rhf_capped.selection.discarded_squared_weight ≈ 0.0199 atol=2e-15
    @test rhf_capped.selection.limiting_reason === :maximum_rank

    mixed = HFDMRG._factor_uhf_state_link(
        Diagonal([0.02, 0.99]), 0, 2, capped)
    @test mixed.selection.completed == mixed.selection.empty == 1
    @test mixed.selection.discarded_squared_weight ≈ 1 - 0.99*0.98 atol=2e-15

    exact = HFDMRG.RHFStateControl(1, cutoff; exact=true)
    for factor in (HFDMRG._factor_state_link, HFDMRG._factor_uhf_state_link)
        link = factor(restriction, 0, 1, exact)
        @test link.active_rank == 1
        @test link.selection.limiting_reason === :exact
        @test iszero(link.selection.discarded_squared_weight)
    end

    # The alpha bond has probability 0.99 on the right of the final cut.
    # There is also one localized alpha orbital and one localized beta orbital.
    for cells in (2, 3), spin_swapped in (false, true)
        sites = 18cells
        i = sites - 18
        bands = zeros(2, sites)
        bands[1, :] .= 1.0
        bands[1, i] = 1 - 2*0.01
        bands[1, i+1] = 1 - 2*0.99
        bands[2, i] = -2sqrt(0.01*0.99)
        bands[1, i+2] = -0.5
        bands[1, i+3] = -0.3
        one_body = BandedOneBody(bands)
        operator = UnitCellInteraction(sites, 1, [0.0], zeros(0, 0),
            zeros(18, 1), zeros(1, 18), zeros(171))
        atoms = [1:i+1, i+2:i+2, i+3:sites]
        workspace = HFDMRG.UHFLifecycleWorkspace(one_body, operator, 1, 1)
        state_control = HFDMRG.UHFStateControl(control)
        lifecycle = HFDMRG.initialize_uhf_neel_lifecycle(one_body, operator,
            state_control, workspace; atom_intervals=atoms, spin_swapped)
        @test lifecycle.root.alpha.total_occupied == (spin_swapped ? 1 : 2)
        @test lifecycle.root.beta.total_occupied == (spin_swapped ? 2 : 1)
        expected = spin_swapped ? [0.0, 0.01] : [0.01, 0.0]
        if cells == 3
            suffix = lifecycle.collection[end]
            @test suffix.alpha.rank == suffix.beta.rank == 0
            @test [suffix.alpha.link.selection.discarded_squared_weight,
                suffix.beta.link.selection.discarded_squared_weight] ≈
                expected atol=2e-15
        else
            # With only two cells the initial root contains both physical
            # intervals. The first truncation occurs in fixed-state transport.
            evidence = only(HFDMRG.run_uhf_fixed_half_sweep!(
                lifecycle, one_body, operator, state_control, workspace))
            @test evidence.alpha_active == evidence.beta_active == 0
            @test [evidence.alpha_discarded, evidence.beta_discarded] ≈
                expected atol=2e-15
        end
    end
end
