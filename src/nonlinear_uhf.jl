mutable struct UHFNonlinearWorkspace
    lifecycle::UHFLifecycleWorkspace
    trial_root::UHFMovingRoot
    baseline_root::UHFMovingRoot
    trial_alpha_work::RHFRootWorkspace
    trial_beta_work::RHFRootWorkspace
    alpha_incoming::Matrix{Float64}
    beta_incoming::Matrix{Float64}
    alpha_proposal::Matrix{Float64}
    beta_proposal::Matrix{Float64}
    alpha_geodesic::Matrix{Float64}
    beta_geodesic::Matrix{Float64}
    alpha_fock::Matrix{Float64}
    beta_fock::Matrix{Float64}
    transaction_calls::Int
    eigensolution_calls::Int
    trial_energy_calls::Int
    publication_calls::Int
end

function UHFNonlinearWorkspace(one_body::BandedOneBody,
        operator::UnitCellInteraction, maximum_alpha::Int,
        maximum_beta::Int)
    lifecycle = UHFLifecycleWorkspace(one_body, operator, maximum_alpha,
        maximum_beta)
    alpha_capacity = lifecycle.maximum_alpha
    beta_capacity = lifecycle.maximum_beta
    matrix(capacity) = Matrix{Float64}(undef, capacity, capacity)
    moving(capacity) = RHFMovingRoot(matrix(capacity), 0, 0, 0, 0, 0, 0)
    root_work(capacity) = RHFRootWorkspace(matrix(capacity), matrix(capacity))
    UHFNonlinearWorkspace(lifecycle,
        UHFMovingRoot(moving(alpha_capacity), moving(beta_capacity)),
        UHFMovingRoot(moving(alpha_capacity), moving(beta_capacity)),
        root_work(alpha_capacity), root_work(beta_capacity),
        matrix(alpha_capacity), matrix(beta_capacity),
        matrix(alpha_capacity), matrix(beta_capacity),
        matrix(alpha_capacity), matrix(beta_capacity),
        matrix(alpha_capacity), matrix(beta_capacity), 0, 0, 0, 0)
end

struct UHFCompactState
    lifecycle::UHFFixedLifecycle
end

struct UHFDMRGResult
    converged::Bool
    reason::Symbol
    represented_energy::Float64
    half_sweeps::Int
    accepted_full::Int
    accepted_fallback::Int
    rejected_updates::Int
    maximum_alpha_rank::Int
    maximum_beta_rank::Int
    maximum_alpha_pair_rank::Int
    maximum_beta_pair_rank::Int
    maximum_alpha_discarded_squared_weight::Float64
    maximum_beta_discarded_squared_weight::Float64
    full_sweep_changes::Vector{Float64}
    replay_energy_change::Float64
    replay_alpha_projector_change::Float64
    replay_beta_projector_change::Float64
    alpha_projector_envelope::Float64
    beta_projector_envelope::Float64
    state::UHFCompactState
end

function _copy_uhf_root!(destination::UHFMovingRoot, source::UHFMovingRoot)
    _copy_moving_root!(destination.alpha, source.alpha)
    _copy_moving_root!(destination.beta, source.beta)
    destination
end

function _uhf_candidate_blocks(candidate, lifecycle::UHFFixedLifecycle,
        workspace::UHFLifecycleWorkspace)
    workspace.empty.provenance = lifecycle.provenance
    workspace.empty.alpha.provenance = lifecycle.provenance
    workspace.empty.beta.provenance = lifecycle.provenance
    cell = candidate.next_center
    if lifecycle.cells == 2
        return workspace.empty, workspace.empty
    end
    left = if cell == 1
        workspace.empty
    elseif cell - 1 == candidate.cell
        candidate.block
    else
        lifecycle.collection[cell-1]
    end
    right = if cell == lifecycle.cells - 1
        workspace.empty
    elseif cell + 1 == candidate.cell
        candidate.block
    else
        lifecycle.collection[cell+1]
    end
    left, right
end

function _uhf_candidate_energy!(candidate, lifecycle::UHFFixedLifecycle,
        one_body::BandedOneBody, operator::UnitCellInteraction,
        workspace::UHFLifecycleWorkspace)
    left, right = _uhf_candidate_blocks(candidate, lifecycle, workspace)
    root = candidate.root
    cell = candidate.next_center
    alpha_expected = left.alpha.rank + _center_width_uhf(lifecycle, cell) +
        right.alpha.rank
    beta_expected = left.beta.rank + _center_width_uhf(lifecycle, cell) +
        right.beta.rank
    root.alpha.rank == alpha_expected && root.beta.rank == beta_expected ||
        throw(DimensionMismatch("candidate UHF root does not span the destination center"))
    prepare_uhf_lifecycle_center!(workspace.prepared, left, right, cell,
        lifecycle.intervals[cell], lifecycle.intervals[cell+1], one_body,
        operator)
    uhf_center_true_energy!(workspace.prepared,
        @view(root.alpha.covariance[1:root.alpha.rank, 1:root.alpha.rank]),
        @view(root.beta.covariance[1:root.beta.rank, 1:root.beta.rank]),
        root.alpha.occupied, root.beta.occupied).total
end

function _uhf_trial_advance!(trial::UHFMovingRoot, source::UHFMovingRoot,
        lifecycle::UHFFixedLifecycle, left::UHFEntryBlock,
        right::UHFEntryBlock, alpha_covariance, beta_covariance,
        one_body::BandedOneBody, operator::UnitCellInteraction,
        control::UHFStateControl, workspace::UHFNonlinearWorkspace)
    _copy_uhf_root!(trial, source)
    copyto!(@view(trial.alpha.covariance[1:source.alpha.rank,
        1:source.alpha.rank]), alpha_covariance)
    copyto!(@view(trial.beta.covariance[1:source.beta.rank,
        1:source.beta.rank]), beta_covariance)
    candidate = _uhf_build_candidate!(trial, lifecycle, left, right,
        one_body, operator, control, workspace.lifecycle,
        workspace.trial_alpha_work, workspace.trial_beta_work)
    energy = _uhf_candidate_energy!(candidate, lifecycle, one_body, operator,
        workspace.lifecycle)
    workspace.trial_energy_calls += 1
    candidate, energy
end

@inline function _uhf_provenance_matches(lifecycle, operator, control)
    common = (operator.provenance, lifecycle.intervals,
        lifecycle.atom_intervals, control.alpha.maximum_active,
        control.alpha.cutoff, control.alpha.exact,
        control.beta.maximum_active, control.beta.cutoff,
        control.beta.exact)
    lifecycle.provenance == hash((common..., false)) ||
        lifecycle.provenance == hash((common..., true))
end

function _preflight_uhf_nonlinear(lifecycle::UHFFixedLifecycle,
        one_body::BandedOneBody, operator::UnitCellInteraction,
        control::UHFStateControl, workspace::UHFNonlinearWorkspace)
    _uhf_provenance_matches(lifecycle, operator, control) ||
        throw(ArgumentError("UHF lifecycle controls or operator provenance changed"))
    h1_sites(one_body) == operator.sites || throw(DimensionMismatch(
        "H1 and interaction lengths disagree"))
    lifecycle.root.alpha.rank <= size(workspace.alpha_incoming, 1) &&
        lifecycle.root.beta.rank <= size(workspace.beta_incoming, 1) ||
        throw(DimensionMismatch("UHF nonlinear workspace capacity is too small"))
    _check_uhf_root(lifecycle.root)
    nothing
end

function _monotone_uhf_center_transaction!(lifecycle::UHFFixedLifecycle,
        one_body::BandedOneBody, operator::UnitCellInteraction,
        control::UHFStateControl, workspace::UHFNonlinearWorkspace;
        geodesic_steps=_RHF_GEODESIC_STEPS)
    _preflight_uhf_nonlinear(lifecycle, one_body, operator, control, workspace)
    workspace.transaction_calls += 1
    localwork = workspace.lifecycle
    left, right = _uhf_lifecycle_blocks(lifecycle, localwork)
    cell, forward = lifecycle.next_center, lifecycle.forward
    alpha_root, beta_root = lifecycle.root.alpha, lifecycle.root.beta
    ma, mb = alpha_root.rank, beta_root.rank
    prepare_uhf_lifecycle_center!(localwork.prepared, left, right, cell,
        lifecycle.intervals[cell], lifecycle.intervals[cell+1], one_body,
        operator)
    incoming = uhf_center_energy_fock!(localwork.prepared,
        @view(alpha_root.covariance[1:ma, 1:ma]),
        @view(beta_root.covariance[1:mb, 1:mb]), alpha_root.occupied,
        beta_root.occupied)
    copyto!(@view(workspace.alpha_incoming[1:ma, 1:ma]),
        @view(alpha_root.covariance[1:ma, 1:ma]))
    copyto!(@view(workspace.beta_incoming[1:mb, 1:mb]),
        @view(beta_root.covariance[1:mb, 1:mb]))
    copyto!(@view(workspace.alpha_fock[1:ma, 1:ma]),
        @view(localwork.prepared.alpha.fock[1:ma, 1:ma]))
    copyto!(@view(workspace.beta_fock[1:mb, 1:mb]),
        @view(localwork.prepared.beta.fock[1:mb, 1:mb]))
    _lowest_projector!(workspace.alpha_proposal, workspace.alpha_fock, ma,
        alpha_root.occupied)
    _lowest_projector!(workspace.beta_proposal, workspace.beta_fock, mb,
        beta_root.occupied)
    workspace.eigensolution_calls += 2
    alpha_change = _projector_distance(workspace.alpha_incoming,
        workspace.alpha_proposal, ma)
    beta_change = _projector_distance(workspace.beta_incoming,
        workspace.beta_proposal, mb)

    full, proposal_energy = _uhf_trial_advance!(workspace.trial_root,
        lifecycle.root, lifecycle, left, right,
        @view(workspace.alpha_proposal[1:ma, 1:ma]),
        @view(workspace.beta_proposal[1:mb, 1:mb]), one_body, operator,
        control, workspace)
    envelope = _rhf_energy_envelope(incoming.total, ma + mb)
    if proposal_energy <= incoming.total + envelope
        _publish_uhf_advance_candidate!(lifecycle, full)
        workspace.publication_calls += 1
        lifecycle.center_visits += 1
        return _uhf_update_record(cell, forward, incoming.total,
            incoming.total, proposal_energy, proposal_energy, 1.0,
            :accepted_full, alpha_change, beta_change, full, ma, mb)
    end

    baseline, baseline_energy = _uhf_trial_advance!(workspace.baseline_root,
        lifecycle.root, lifecycle, left, right,
        @view(workspace.alpha_incoming[1:ma, 1:ma]),
        @view(workspace.beta_incoming[1:mb, 1:mb]), one_body, operator,
        control, workspace)
    envelope = _rhf_energy_envelope(baseline_energy, ma + mb)
    for step in geodesic_steps
        _geodesic_projector!(workspace.alpha_geodesic,
            workspace.alpha_incoming, workspace.alpha_proposal, ma,
            alpha_root.occupied, step)
        _geodesic_projector!(workspace.beta_geodesic,
            workspace.beta_incoming, workspace.beta_proposal, mb,
            beta_root.occupied, step)
        trial, energy = _uhf_trial_advance!(workspace.trial_root,
            lifecycle.root, lifecycle, left, right,
            @view(workspace.alpha_geodesic[1:ma, 1:ma]),
            @view(workspace.beta_geodesic[1:mb, 1:mb]), one_body, operator,
            control, workspace)
        if energy <= baseline_energy + envelope
            _publish_uhf_advance_candidate!(lifecycle, trial)
            workspace.publication_calls += 1
            lifecycle.center_visits += 1
            return _uhf_update_record(cell, forward, incoming.total,
                baseline_energy, proposal_energy, energy, step,
                :accepted_geodesic, alpha_change, beta_change, trial, ma, mb)
        end
    end
    _publish_uhf_advance_candidate!(lifecycle, baseline)
    workspace.publication_calls += 1
    lifecycle.center_visits += 1
    _uhf_update_record(cell, forward, incoming.total, baseline_energy,
        proposal_energy, baseline_energy, 0.0, :rejected_no_safe_trial,
        alpha_change, beta_change, baseline, ma, mb)
end

@inline function _uhf_update_record(cell, forward, incoming, baseline,
        proposal, published, step, status, alpha_change, beta_change,
        candidate, alpha_center_rank, beta_center_rank)
    (cell=cell, forward=forward, incoming_energy=incoming,
        baseline_energy=baseline, proposal_energy=proposal,
        published_energy=published, step=step, status=status,
        alpha_self_change=alpha_change, beta_self_change=beta_change,
        alpha_discarded=candidate.alpha_selection.discarded_squared_weight,
        beta_discarded=candidate.beta_selection.discarded_squared_weight,
        alpha_rank=candidate.block.alpha.rank,
        beta_rank=candidate.block.beta.rank,
        alpha_center_rank=alpha_center_rank,
        beta_center_rank=beta_center_rank,
        alpha_pair_rank=pair_dimension(candidate.block.alpha.rank),
        beta_pair_rank=pair_dimension(candidate.block.beta.rank))
end

function _uhf_invariant_energy!(lifecycle::UHFFixedLifecycle,
        one_body::BandedOneBody, operator::UnitCellInteraction,
        workspace::UHFLifecycleWorkspace)
    left, right = _uhf_lifecycle_blocks(lifecycle, workspace)
    cell = lifecycle.next_center
    prepare_uhf_lifecycle_center!(workspace.prepared, left, right, cell,
        lifecycle.intervals[cell], lifecycle.intervals[cell+1], one_body,
        operator)
    root = lifecycle.root
    uhf_center_energy_fock!(workspace.prepared,
        @view(root.alpha.covariance[1:root.alpha.rank, 1:root.alpha.rank]),
        @view(root.beta.covariance[1:root.beta.rank, 1:root.beta.rank]),
        root.alpha.occupied, root.beta.occupied).total
end

function _uhf_nonlinear_half_sweep!(lifecycle, one_body, operator, control,
        workspace)
    orientation = lifecycle.forward
    evidence = NamedTuple[]
    while lifecycle.forward == orientation
        push!(evidence, _monotone_uhf_center_transaction!(lifecycle,
            one_body, operator, control, workspace))
    end
    lifecycle.half_sweeps += 1
    evidence, _uhf_invariant_energy!(lifecycle, one_body, operator,
        workspace.lifecycle)
end

function _uhf_solver_loop!(lifecycle, one_body, operator, state_control,
        energy_tolerance::Float64, maximum_half_sweeps::Int,
        workspace::UHFNonlinearWorkspace)
    energies = Float64[]; changes = Float64[]; all_updates = NamedTuple[]
    full = fallback = rejected = 0
    max_alpha = max_beta = max_alpha_pair = max_beta_pair = 0
    max_alpha_discarded = max_beta_discarded = 0.0
    quiet = 0; converged = false; reason = :maximum_half_sweeps
    for half = 1:maximum_half_sweeps
        updates, energy = _uhf_nonlinear_half_sweep!(lifecycle, one_body,
            operator, state_control, workspace)
        append!(all_updates, updates); push!(energies, energy)
        for update in updates
            full += update.status === :accepted_full
            fallback += update.status === :accepted_geodesic
            rejected += update.status === :rejected_no_safe_trial
            max_alpha = max(max_alpha, update.alpha_rank)
            max_beta = max(max_beta, update.beta_rank)
            max_alpha_pair = max(max_alpha_pair, update.alpha_pair_rank)
            max_beta_pair = max(max_beta_pair, update.beta_pair_rank)
            max_alpha_discarded = max(max_alpha_discarded,
                update.alpha_discarded)
            max_beta_discarded = max(max_beta_discarded,
                update.beta_discarded)
            update.published_energy <= update.baseline_energy +
                _rhf_energy_envelope(update.baseline_energy,
                    update.alpha_center_rank + update.beta_center_rank) ||
                error("published UHF optimizer update violates monotonicity")
        end
        if half >= 3
            change = energy - energies[end-2]
            push!(changes, change)
            quiet = abs(change) <= energy_tolerance ? quiet + 1 : 0
            if quiet >= 2
                converged = true; reason = :energy_converged
                break
            end
        end
    end
    (converged=converged, reason=reason, energies=energies, changes=changes,
        full=full, fallback=fallback, rejected=rejected,
        max_alpha=max_alpha, max_beta=max_beta,
        max_alpha_pair=max_alpha_pair, max_beta_pair=max_beta_pair,
        max_alpha_discarded=max_alpha_discarded,
        max_beta_discarded=max_beta_discarded, updates=all_updates)
end

function _solve_uhf_hfdmrg(one_body::BandedOneBody,
        operator::UnitCellInteraction; state_maxdim::Int,
        state_cutoff::Real, energy_tolerance::Real,
        maximum_half_sweeps::Int, starting_orientation::Symbol,
        exact::Bool, atom_intervals, spin_swapped::Bool)
    control = UHFStateControl(RHFStateControl(state_maxdim,
        Float64(state_cutoff); exact))
    atoms = atom_intervals === nothing ?
        length(_traversal_intervals(operator)) : length(atom_intervals)
    alpha_total = spin_swapped ? atoms ÷ 2 : (atoms + 1) ÷ 2
    beta_total = atoms - alpha_total
    alpha_capacity = exact ? max(state_maxdim, alpha_total) : state_maxdim
    beta_capacity = exact ? max(state_maxdim, beta_total) : state_maxdim
    workspace = UHFNonlinearWorkspace(one_body, operator, alpha_capacity,
        beta_capacity)
    lifecycle = initialize_uhf_neel_lifecycle(one_body, operator, control,
        workspace.lifecycle; forward=starting_orientation === :physical,
        atom_intervals, spin_swapped)
    outcome = _uhf_solver_loop!(lifecycle, one_body, operator, control,
        Float64(energy_tolerance), maximum_half_sweeps, workspace)

    replay_before = _uhf_invariant_energy!(lifecycle, one_body, operator,
        workspace.lifecycle)
    alpha_rank = lifecycle.root.alpha.rank
    beta_rank = lifecycle.root.beta.rank
    copyto!(@view(workspace.alpha_incoming[1:alpha_rank, 1:alpha_rank]),
        @view(lifecycle.root.alpha.covariance[1:alpha_rank, 1:alpha_rank]))
    copyto!(@view(workspace.beta_incoming[1:beta_rank, 1:beta_rank]),
        @view(lifecycle.root.beta.covariance[1:beta_rank, 1:beta_rank]))
    for _ = 1:2
        run_uhf_fixed_half_sweep!(lifecycle, one_body, operator, control,
            workspace.lifecycle)
    end
    replay_after = _uhf_invariant_energy!(lifecycle, one_body, operator,
        workspace.lifecycle)
    alpha_projector = lifecycle.root.alpha.rank == alpha_rank ?
        _projector_distance(workspace.alpha_incoming,
            lifecycle.root.alpha.covariance, alpha_rank) : Inf
    beta_projector = lifecycle.root.beta.rank == beta_rank ?
        _projector_distance(workspace.beta_incoming,
            lifecycle.root.beta.covariance, beta_rank) : Inf
    alpha_envelope = 4sqrt(outcome.max_alpha_discarded) +
        8192eps(Float64)*max(lifecycle.root.alpha.rank, 1)
    beta_envelope = 4sqrt(outcome.max_beta_discarded) +
        8192eps(Float64)*max(lifecycle.root.beta.rank, 1)
    converged, reason = outcome.converged, outcome.reason
    replay_energy = replay_after - replay_before
    envelope = _rhf_energy_envelope(replay_before,
        lifecycle.root.alpha.rank + lifecycle.root.beta.rank)
    if exact
        if converged && abs(replay_energy) > envelope
            converged = false; reason = :exact_replay_failed
        end
    elseif converged && abs(replay_energy) > Float64(energy_tolerance) +
            envelope
        converged = false; reason = :finite_replay_floor
    end
    UHFDMRGResult(converged, reason, replay_after, length(outcome.energies),
        outcome.full, outcome.fallback, outcome.rejected,
        outcome.max_alpha, outcome.max_beta, outcome.max_alpha_pair,
        outcome.max_beta_pair, outcome.max_alpha_discarded,
        outcome.max_beta_discarded, outcome.changes, replay_energy,
        alpha_projector, beta_projector, alpha_envelope, beta_envelope,
        UHFCompactState(lifecycle))
end
