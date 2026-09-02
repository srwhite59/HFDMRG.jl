function solve_fixture(atoms, orientation; exact=false, sweeps=60,
        tolerance=1.0e-7)
    one_body, operator = nucleus_fixture(atoms)
    result = HFDMRG.solve_hfdmrg(one_body, operator;
        state_maxdim=exact ? atoms ÷ 2 : min(4, atoms ÷ 2),
        state_cutoff=1.0e-12, energy_tolerance=tolerance,
        maximum_half_sweeps=sweeps, starting_orientation=orientation, exact)
    result, one_body, operator
end

function compact_projector(result, workspace)
    lifecycle = result.state.lifecycle
    output = zeros(lifecycle.cells * 18, lifecycle.root.total_occupied)
    HFDMRG.reconstruct_terminal!(output, lifecycle, workspace)
    output * transpose(output)
end

function warmed_publication_allocation(lifecycle, candidate)
    HFDMRG._publish_advance_candidate!(lifecycle, candidate)
    @allocated HFDMRG._publish_advance_candidate!(lifecycle, candidate)
end

@testset "standalone monotone RHF facade" begin
    energies = Dict{Tuple{Int,Bool,Symbol},Float64}()
    for atoms in (2, 10), exact in (false, true),
            orientation in (:physical, :reflected)
        result, one_body, operator = solve_fixture(atoms, orientation; exact)
        @test result.converged
        @test result.reason == :energy_converged
        @test result.half_sweeps <= 60
        @test result.accepted_full + result.accepted_fallback > 0
        @test result.rejected_updates >= 0
        @test result.maximum_state_rank <= (exact ? atoms ÷ 2 : 4)
        @test result.maximum_pair_rank ==
            HFDMRG.pair_dimension(result.maximum_state_rank)
        @test result.maximum_discarded_squared_weight >= 0.0
        @test abs(result.replay_energy_change) <= 1.0e-10
        @test result.projector_envelope >=
            4sqrt(result.maximum_discarded_squared_weight)
        lifecycle = result.state.lifecycle
        workspace = HFDMRG.RHFLifecycleWorkspace(one_body, operator,
            exact ? atoms ÷ 2 : min(4, atoms ÷ 2))
        coefficients = zeros(18atoms, atoms ÷ 2)
        HFDMRG.reconstruct_terminal!(coefficients, lifecycle, workspace)
        @test norm(transpose(coefficients)*coefficients-I, Inf) <= 2e-11
        density, energy, _ = dense_rhf_state(one_body, operator, coefficients)
        @test energy ≈ result.represented_energy atol=2e-10
        @test norm(density*density-density, Inf) <= 2e-11
        energies[(atoms, exact, orientation)] = result.represented_energy
    end
    for atoms in (2, 10), exact in (false, true)
        # Representation exactness leaves the caller's 1e-7 nonlinear
        # stopping target unchanged.  Two independently started solves have
        # the sum of those stopping envelopes; each exact replay above remains
        # constrained at roundoff scale.
        @test abs(energies[(atoms, exact, :physical)] -
            energies[(atoms, exact, :reflected)]) <= 2.0e-7
    end
end

@testset "H20 finite convergence and matched publication" begin
    results = HFDMRG.HFDMRGResult[]
    for orientation in (:physical, :reflected)
        result, one_body, operator = solve_fixture(20, orientation;
            sweeps=60, tolerance=1.0e-7)
        @test result.converged
        @test result.replay_energy_change <= 1.0e-7
        push!(results, result)
        lifecycle = result.state.lifecycle
        @test lifecycle.original_rows_read == 0
        @test length(lifecycle.collection) == 19
        @test lifecycle.root.total_occupied == 10
    end
    # Each independently started finite solve certifies its own 1e-7 working-
    # energy envelope.  Their crossed comparison therefore has the sum of the
    # two certified envelopes; exact-mode orientation parity is tested below.
    @test abs(results[1].represented_energy-results[2].represented_energy) <=
        2.0e-7

    one_body, operator = nucleus_fixture(10)
    control = HFDMRG.RHFStateControl(4, 1.0e-12)
    workspace = HFDMRG.RHFNonlinearWorkspace(one_body, operator, 4)
    lifecycle = HFDMRG.initialize_rhf_dimer_lifecycle(one_body, operator,
        control, workspace.lifecycle)
    evidence, _ = HFDMRG._nonlinear_half_sweep!(lifecycle, one_body,
        operator, control, workspace)
    @test all(item.published_energy <= item.baseline_energy +
        HFDMRG._rhf_energy_envelope(item.baseline_energy, item.center_rank)
        for item in evidence)
    @test all(item.status in (:accepted_full, :accepted_geodesic,
        :rejected_no_safe_trial) for item in evidence)
end

@testset "H20 exact reference orientations" begin
    energies = Float64[]
    for orientation in (:physical, :reflected)
        result, _, _ = solve_fixture(20, orientation; exact=true,
            sweeps=60, tolerance=1.0e-7)
        @test result.converged
        @test result.reason == :energy_converged
        @test abs(result.replay_energy_change) <= 1.0e-10
        @test result.maximum_discarded_squared_weight == 0.0
        push!(energies, result.represented_energy)
    end
    @test abs(energies[1]-energies[2]) <= 1.0e-7
end

@testset "stationary, perturbed, insufficient, and atomic paths" begin
    result, one_body, operator = solve_fixture(10, :physical)
    control = HFDMRG.RHFStateControl(4, 1.0e-12)
    workspace = HFDMRG.RHFNonlinearWorkspace(one_body, operator, 4)
    stationary = result.state.lifecycle
    energy_before = HFDMRG._invariant_energy!(stationary, one_body, operator,
        workspace.lifecycle)
    replay = HFDMRG._solver_loop!(stationary, one_body, operator, control,
        1.0e-7, 8, workspace)
    @test replay[1]
    @test abs(last(replay[3]) - energy_before) <= 1.0e-7

    perturbed = HFDMRG.initialize_rhf_dimer_lifecycle(one_body, operator,
        control, workspace.lifecycle)
    rank = perturbed.root.rank
    incoming = copy(@view perturbed.root.covariance[1:rank, 1:rank])
    factors = eigen(Symmetric(incoming))
    basis = Matrix(@view factors.vectors[:, rank-perturbed.root.occupied+1:rank])
    virtual = factors.vectors[:, 1]
    angle = 0.04
    basis[:, 1] .= cos(angle)*basis[:, 1] + sin(angle)*virtual
    mul!(@view(perturbed.root.covariance[1:rank, 1:rank]), basis,
        transpose(basis))
    HFDMRG._check_moving_root(perturbed.root)
    perturbed_outcome = HFDMRG._solver_loop!(perturbed, one_body, operator,
        control, 1.0e-7, 60, workspace)
    @test perturbed_outcome[1]
    @test last(perturbed_outcome[3]) <= energy_before + 1.0e-7

    insufficient = HFDMRG.solve_hfdmrg(one_body, operator;
        state_maxdim=4, state_cutoff=1.0e-12,
        energy_tolerance=1.0e-11, maximum_half_sweeps=2)
    @test !insufficient.converged
    @test insufficient.reason == :maximum_half_sweeps

    invalid_workspace = HFDMRG.RHFNonlinearWorkspace(one_body, operator, 4)
    invalid_lifecycle = HFDMRG.initialize_rhf_dimer_lifecycle(one_body,
        operator, control, invalid_workspace.lifecycle)
    covariance_before = copy(invalid_lifecycle.root.covariance)
    collection_before = copy(invalid_lifecycle.collection)
    wrong = HFDMRG.RHFStateControl(3, 1.0e-12)
    @test_throws ArgumentError HFDMRG._monotone_center_transaction!(
        invalid_lifecycle, one_body, operator, wrong, invalid_workspace)
    @test isequal(covariance_before, invalid_lifecycle.root.covariance)
    @test collection_before == invalid_lifecycle.collection
    HFDMRG._preflight_nonlinear(invalid_lifecycle, one_body, operator,
        control, invalid_workspace)
    @test @allocated(HFDMRG._preflight_nonlinear(invalid_lifecycle, one_body,
        operator, control, invalid_workspace)) == 0
    left, right = HFDMRG._lifecycle_blocks(invalid_lifecycle,
        invalid_workspace.lifecycle)
    HFDMRG._copy_moving_root!(invalid_workspace.trial_root,
        invalid_lifecycle.root)
    candidate = HFDMRG._build_advance_candidate!(invalid_workspace.trial_root,
        invalid_lifecycle, left, right, one_body, operator, control,
        invalid_workspace.lifecycle, invalid_workspace.trial_work)
    @test warmed_publication_allocation(invalid_lifecycle, candidate) == 0
    @test_throws ArgumentError HFDMRG.solve_hfdmrg(one_body, operator;
        state_maxdim=4, state_cutoff=1.0e-12,
        energy_tolerance=1.0e-7, maximum_half_sweeps=2,
        starting_orientation=:sideways)
end
