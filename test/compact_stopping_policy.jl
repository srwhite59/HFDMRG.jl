# Packet 10A7D: compact extensive-energy stopping policy.

@testset "compact extensive-energy threshold arithmetic" begin
    for requested in (1.0e-7, 1.0e-9, 1.0e-11)
        threshold = HFDMRG._energy_stopping_threshold(requested, -10.0, 4)
        @test threshold.requested == requested
        @test threshold.effective == requested
        @test threshold.reason == :requested_energy_tolerance
    end

    requested = 1.0e-11
    energy = -1.109e5
    rank = 5
    expected = 256eps(Float64) * abs(energy) * sqrt(rank)
    extensive = HFDMRG._energy_stopping_threshold(requested, energy, rank)
    @test extensive.requested == requested
    @test extensive.effective == expected
    @test extensive.reason == :float64_extensive_energy_floor
    @test HFDMRG._energy_change_is_quiet(expected, extensive.effective)
    @test !HFDMRG._energy_change_is_quiet(nextfloat(expected),
        extensive.effective)
    @test !HFDMRG._energy_change_is_quiet(NaN, extensive.effective)

    low = HFDMRG._energy_stopping_threshold(1.0e-7, energy, rank)
    @test low.effective == 1.0e-7
    @test low.reason == :requested_energy_tolerance
    @test_throws ArgumentError HFDMRG._energy_stopping_threshold(0.0,
        energy, rank)
    @test_throws ArgumentError HFDMRG._energy_stopping_threshold(requested,
        Inf, rank)
    @test_throws ArgumentError HFDMRG._energy_stopping_threshold(requested,
        energy, -1)
end

@testset "compact RHF and UHF stopping diagnostics" begin
    requested = 1.0e-11
    rhf, one_body, operator = solve_fixture(2, :physical; sweeps=30,
        tolerance=requested)
    @test rhf.converged
    @test rhf.requested_energy_tolerance == requested
    @test rhf.effective_energy_tolerance >= requested
    @test rhf.stopping_tolerance_reason in
        (:requested_energy_tolerance, :float64_extensive_energy_floor)
    @test length(rhf.full_sweep_changes) >= 2
    @test all(abs.(rhf.full_sweep_changes[end-1:end]) .<=
        rhf.effective_energy_tolerance)

    uhf, _, _ = solve_uhf_fixture(2, :physical; sweeps=30,
        tolerance=requested)
    @test uhf.converged
    @test uhf.requested_energy_tolerance == requested
    @test uhf.effective_energy_tolerance >= requested
    @test uhf.stopping_tolerance_reason in
        (:requested_energy_tolerance, :float64_extensive_energy_floor)
    @test length(uhf.full_sweep_changes) >= 2
    @test all(abs.(uhf.full_sweep_changes[end-1:end]) .<=
        uhf.effective_energy_tolerance)

    rhf_limited = HFDMRG.solve_hfdmrg(one_body, operator;
        state_maxdim=1, state_cutoff=1.0e-12,
        energy_tolerance=requested, maximum_half_sweeps=2)
    @test !rhf_limited.converged
    @test rhf_limited.reason == :maximum_half_sweeps
    @test rhf_limited.requested_energy_tolerance == requested
    @test isfinite(rhf_limited.effective_energy_tolerance)

    uhf_limited = HFDMRG.solve_hfdmrg(one_body, operator; spin=:uhf,
        state_maxdim=1, state_cutoff=1.0e-12,
        energy_tolerance=requested, maximum_half_sweeps=2)
    @test !uhf_limited.converged
    @test uhf_limited.reason == :maximum_half_sweeps
    @test uhf_limited.requested_energy_tolerance == requested
    @test isfinite(uhf_limited.effective_energy_tolerance)

    exact, _, _ = solve_fixture(2, :reflected; exact=true, sweeps=30,
        tolerance=requested)
    @test exact.converged
    @test exact.maximum_discarded_squared_weight == 0.0
    @test abs(exact.replay_energy_change) <=
        HFDMRG._rhf_energy_envelope(exact.represented_energy,
            exact.state.lifecycle.root.rank)
end
