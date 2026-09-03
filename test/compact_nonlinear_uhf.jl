function solve_uhf_fixture(atoms, orientation; exact=false, sweeps=60,
        tolerance=1.0e-7, spin_swapped=false)
    one_body, operator = nucleus_fixture(atoms)
    result = HFDMRG.solve_hfdmrg(one_body, operator; spin=:uhf,
        initialization=:atomic_neel,
        state_maxdim=exact ? (atoms + 1) ÷ 2 : min(4, (atoms + 1) ÷ 2),
        state_cutoff=1.0e-12, energy_tolerance=tolerance,
        maximum_half_sweeps=sweeps, starting_orientation=orientation, exact,
        spin_swapped)
    result, one_body, operator
end

function reconstructed_uhf_projectors(result, one_body, operator)
    lifecycle = result.state.lifecycle
    capacity = max(result.maximum_alpha_rank, result.maximum_beta_rank, 1)
    workspace = HFDMRG.UHFLifecycleWorkspace(one_body, operator, capacity,
        capacity)
    alpha = zeros(operator.sites, lifecycle.root.alpha.total_occupied)
    beta = zeros(operator.sites, lifecycle.root.beta.total_occupied)
    HFDMRG.reconstruct_uhf_terminal!(alpha, beta, lifecycle, workspace)
    alpha * transpose(alpha), beta * transpose(beta), alpha, beta
end

function scaled_interaction(operator, scale)
    HFDMRG.UnitCellInteraction(operator.sites, operator.phase_at_first,
        operator.tail_transfer, operator.residual_transfer,
        scale * operator.phase_source, operator.phase_sink,
        scale * operator.local_triangle)
end

@testset "standalone monotone UHF facade" begin
    energies = Dict{Tuple{Int,Bool,Symbol},Float64}()
    projectors = Dict{Tuple{Int,Bool,Symbol},Tuple{Matrix{Float64},Matrix{Float64}}}()
    for atoms in (2, 10), exact in (false, true),
            orientation in (:physical, :reflected)
        result, one_body, operator = solve_uhf_fixture(atoms, orientation;
            exact)
        @test result isa HFDMRG.UHFDMRGResult
        @test result.converged
        @test result.reason == :energy_converged
        @test result.accepted_full + result.accepted_fallback > 0
        @test result.maximum_alpha_pair_rank ==
            HFDMRG.pair_dimension(result.maximum_alpha_rank)
        @test result.maximum_beta_pair_rank ==
            HFDMRG.pair_dimension(result.maximum_beta_rank)
        @test result.alpha_projector_envelope >=
            4sqrt(result.maximum_alpha_discarded_squared_weight)
        @test result.beta_projector_envelope >=
            4sqrt(result.maximum_beta_discarded_squared_weight)
        Pa, Pb, alpha, beta = reconstructed_uhf_projectors(result, one_body,
            operator)
        @test transpose(alpha) * alpha ≈ I atol=3e-11
        @test transpose(beta) * beta ≈ I atol=3e-11
        expected, _, _, _, _ = dense_uhf_oracle(one_body, operator,
            Matrix{Float64}(I, operator.sites, operator.sites),
            Matrix{Float64}(I, operator.sites, operator.sites), Pa, Pb)
        @test result.represented_energy ≈ expected atol=3e-10
        energies[(atoms, exact, orientation)] = result.represented_energy
        projectors[(atoms, exact, orientation)] = (Pa, Pb)
    end
    for atoms in (2, 10), exact in (false, true)
        @test abs(energies[(atoms, exact, :physical)] -
            energies[(atoms, exact, :reflected)]) <= 2.0e-7
        if exact
            # Exact representation removes truncation; it does not tighten
            # the caller's 1e-7 nonlinear energy target.  Independently
            # started determinants are therefore compared on its square-root
            # projector scale, while each exact replay remains roundoff-gated.
            @test maximum(abs, projectors[(atoms, exact, :physical)][1] -
                projectors[(atoms, exact, :reflected)][1]) <= 1e-5
            @test maximum(abs, projectors[(atoms, exact, :physical)][2] -
                projectors[(atoms, exact, :reflected)][2]) <= 1e-5
        end
    end
end

@testset "H20 finite/exact UHF convergence and spin covariance" begin
    for exact in (false, true)
        physical, one_body, operator = solve_uhf_fixture(20, :physical;
            exact)
        reflected, _, _ = solve_uhf_fixture(20, :reflected; exact)
        @test physical.converged && reflected.converged
        @test abs(physical.represented_energy-reflected.represented_energy) <=
            2.0e-7
        @test physical.state.lifecycle.original_rows_read == 0
        @test reflected.state.lifecycle.original_rows_read == 0
        @test length(physical.state.lifecycle.collection) == 19
        exact && @test max(physical.maximum_alpha_discarded_squared_weight,
            physical.maximum_beta_discarded_squared_weight) == 0.0
        Pa, Pb, _, _ = reconstructed_uhf_projectors(physical, one_body,
            operator)
        @test norm(Pa*Pa-Pa, Inf) <= 3e-10
        @test norm(Pb*Pb-Pb, Inf) <= 3e-10
    end

    ordinary, one_body, operator = solve_uhf_fixture(10, :physical)
    swapped, _, _ = solve_uhf_fixture(10, :physical; spin_swapped=true)
    Pa, Pb, _, _ = reconstructed_uhf_projectors(ordinary, one_body, operator)
    Sa, Sb, _, _ = reconstructed_uhf_projectors(swapped, one_body, operator)
    @test ordinary.represented_energy ≈ swapped.represented_energy atol=2e-10
    @test Pa ≈ Sb atol=3e-8
    @test Pb ≈ Sa atol=3e-8

    lifecycle = ordinary.state.lifecycle
    alpha_before = lifecycle.root.alpha
    beta_before = lifecycle.root.beta
    HFDMRG._spin_swap_uhf_lifecycle!(lifecycle, operator)
    @test lifecycle.root.alpha === beta_before
    @test lifecycle.root.beta === alpha_before
    @test all(block.alpha.left_core === block.beta.left_core
        for block in lifecycle.collection)
    HFDMRG._spin_swap_uhf_lifecycle!(lifecycle, operator)
    @test lifecycle.root.alpha === alpha_before
    @test lifecycle.root.beta === beta_before
end

@testset "UHF full, common-geodesic, rejection, and atomic paths" begin
    one_body, base = nucleus_fixture(2)
    state = HFDMRG.RHFStateControl(2, 1e-12)
    control = HFDMRG.UHFStateControl(state)
    rows = [1:18, 19:36]

    full_workspace = HFDMRG.UHFNonlinearWorkspace(one_body, base, 2, 2)
    full_lifecycle = HFDMRG.initialize_uhf_neel_lifecycle(one_body, base,
        control, full_workspace.lifecycle; atom_intervals=rows)
    full = HFDMRG._monotone_uhf_center_transaction!(full_lifecycle,
        one_body, base, control, full_workspace)
    @test full.status == :accepted_full
    @test full.step == 1.0
    @test full.published_energy <= full.baseline_energy +
        HFDMRG._rhf_energy_envelope(full.baseline_energy,
            full.alpha_center_rank + full.beta_center_rank)

    hostile = scaled_interaction(base, -5.0)
    fallback_workspace = HFDMRG.UHFNonlinearWorkspace(one_body, hostile, 2, 2)
    fallback_lifecycle = HFDMRG.initialize_uhf_neel_lifecycle(one_body,
        hostile, control, fallback_workspace.lifecycle; atom_intervals=rows)
    fallback = HFDMRG._monotone_uhf_center_transaction!(fallback_lifecycle,
        one_body, hostile, control, fallback_workspace)
    @test fallback.status == :accepted_geodesic
    @test 0.0 < fallback.step < 1.0
    swapped_workspace = HFDMRG.UHFNonlinearWorkspace(one_body, hostile, 2, 2)
    swapped_lifecycle = HFDMRG.initialize_uhf_neel_lifecycle(one_body,
        hostile, control, swapped_workspace.lifecycle; atom_intervals=rows,
        spin_swapped=true)
    swapped_fallback = HFDMRG._monotone_uhf_center_transaction!(
        swapped_lifecycle, one_body, hostile, control, swapped_workspace)
    @test swapped_fallback.status == :accepted_geodesic
    @test swapped_fallback.step == fallback.step

    rejected_workspace = HFDMRG.UHFNonlinearWorkspace(one_body, hostile, 2, 2)
    rejected_lifecycle = HFDMRG.initialize_uhf_neel_lifecycle(one_body,
        hostile, control, rejected_workspace.lifecycle; atom_intervals=rows)
    rejected = HFDMRG._monotone_uhf_center_transaction!(rejected_lifecycle,
        one_body, hostile, control, rejected_workspace; geodesic_steps=())
    @test rejected.status == :rejected_no_safe_trial
    @test rejected.step == 0.0
    @test rejected.published_energy == rejected.baseline_energy

    atomic_workspace = HFDMRG.UHFNonlinearWorkspace(one_body, base, 2, 2)
    atomic = HFDMRG.initialize_uhf_neel_lifecycle(one_body, base, control,
        atomic_workspace.lifecycle; atom_intervals=rows)
    alpha_before = copy(atomic.root.alpha.covariance)
    beta_before = copy(atomic.root.beta.covariance)
    collection_before = copy(atomic.collection)
    wrong = HFDMRG.UHFStateControl(HFDMRG.RHFStateControl(1, 1e-12))
    @test_throws ArgumentError HFDMRG._monotone_uhf_center_transaction!(
        atomic, one_body, base, wrong, atomic_workspace)
    @test isequal(alpha_before, atomic.root.alpha.covariance)
    @test isequal(beta_before, atomic.root.beta.covariance)
    @test all(atomic.collection[index] === collection_before[index]
        for index in eachindex(collection_before))
    HFDMRG._preflight_uhf_nonlinear(atomic, one_body, base, control,
        atomic_workspace)
    @test @allocated(HFDMRG._preflight_uhf_nonlinear(atomic, one_body, base,
        control, atomic_workspace)) == 0
    @test @allocated(HFDMRG._copy_uhf_root!(atomic_workspace.trial_root,
        atomic.root)) == 0
    left, right = HFDMRG._uhf_lifecycle_blocks(atomic,
        atomic_workspace.lifecycle)
    HFDMRG._copy_uhf_root!(atomic_workspace.trial_root, atomic.root)
    candidate = HFDMRG._uhf_build_candidate!(atomic_workspace.trial_root,
        atomic, left, right, one_body, base, control,
        atomic_workspace.lifecycle, atomic_workspace.trial_alpha_work,
        atomic_workspace.trial_beta_work)
    HFDMRG._publish_uhf_advance_candidate!(atomic, candidate)
    @test @allocated(HFDMRG._publish_uhf_advance_candidate!(atomic,
        candidate)) == 0
end

@testset "UHF stationary, perturbed, insufficient, and fragment routes" begin
    result, one_body, operator = solve_uhf_fixture(10, :physical)
    state = HFDMRG.RHFStateControl(4, 1e-12)
    control = HFDMRG.UHFStateControl(state)
    workspace = HFDMRG.UHFNonlinearWorkspace(one_body, operator, 4, 4)
    stationary = result.state.lifecycle
    before = HFDMRG._uhf_invariant_energy!(stationary, one_body, operator,
        workspace.lifecycle)
    replay = HFDMRG._uhf_solver_loop!(stationary, one_body, operator,
        control, 1e-7, 8, workspace)
    @test replay.converged
    @test abs(last(replay.energies)-before) <= 1e-7

    perturbed = HFDMRG.initialize_uhf_neel_lifecycle(one_body, operator,
        control, workspace.lifecycle)
    for (root, angle) in ((perturbed.root.alpha, 0.04),
            (perturbed.root.beta, -0.03))
        rank = root.rank
        factors = eigen(Symmetric(Matrix(@view(root.covariance[1:rank,
            1:rank]))))
        basis = Matrix(@view(factors.vectors[:, rank-root.occupied+1:rank]))
        virtual = factors.vectors[:, 1]
        basis[:, 1] .= cos(angle)*basis[:, 1] + sin(angle)*virtual
        mul!(@view(root.covariance[1:rank, 1:rank]), basis, transpose(basis))
    end
    outcome = HFDMRG._uhf_solver_loop!(perturbed, one_body, operator,
        control, 1e-7, 60, workspace)
    @test outcome.converged

    insufficient = HFDMRG.solve_hfdmrg(one_body, operator; spin=:uhf,
        state_maxdim=4, state_cutoff=1e-12, energy_tolerance=1e-11,
        maximum_half_sweeps=2)
    @test !insufficient.converged
    @test insufficient.reason == :maximum_half_sweeps

    fragment_h1, fragment_operator, atoms = finite_fragment_fixture()
    fragment = HFDMRG.solve_hfdmrg(fragment_h1, fragment_operator; spin=:uhf,
        state_maxdim=4, state_cutoff=1e-12, energy_tolerance=1e-7,
        maximum_half_sweeps=30, atom_intervals=atoms)
    @test fragment.converged
    @test fragment.state.lifecycle.intervals[end] == 37:37
    @test_throws ArgumentError HFDMRG.solve_hfdmrg(one_body, operator;
        spin=:uhf, initialization=:dimer, state_maxdim=4,
        state_cutoff=1e-12, energy_tolerance=1e-7,
        maximum_half_sweeps=2)
end
