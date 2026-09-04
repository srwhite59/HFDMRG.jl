const _RHF_GEODESIC_STEPS = (0.9375, 0.875, 0.75, 0.5, 0.25,
    0.125, 0.0625, 0.03125, 0.015625, 0.0078125)

mutable struct RHFNonlinearWorkspace
    lifecycle::RHFLifecycleWorkspace
    trial_root::RHFMovingRoot
    baseline_root::RHFMovingRoot
    trial_work::RHFRootWorkspace
    incoming::Matrix{Float64}
    proposal::Matrix{Float64}
    geodesic::Matrix{Float64}
    first_fock::Matrix{Float64}
    transaction_calls::Int
    eigensolution_calls::Int
    trial_energy_calls::Int
    publication_calls::Int
end

function RHFNonlinearWorkspace(one_body::BandedOneBody,
        operator::UnitCellInteraction, maximum_active::Int)
    lifecycle = RHFLifecycleWorkspace(one_body, operator, maximum_active)
    capacity = lifecycle.maximum_rank
    matrix() = Matrix{Float64}(undef, capacity, capacity)
    trial = RHFMovingRoot(matrix(), 0, 0, 0, 0, 0, 0)
    baseline = RHFMovingRoot(matrix(), 0, 0, 0, 0, 0, 0)
    root_work = RHFRootWorkspace(matrix(), matrix())
    RHFNonlinearWorkspace(lifecycle, trial, baseline, root_work, matrix(),
        matrix(), matrix(), matrix(), 0, 0, 0, 0)
end

struct RHFCompactState
    lifecycle::RHFFixedLifecycle
end


struct HFDMRGResult
    converged::Bool
    reason::Symbol
    represented_energy::Float64
    half_sweeps::Int
    accepted_full::Int
    accepted_fallback::Int
    rejected_updates::Int
    maximum_state_rank::Int
    maximum_pair_rank::Int
    maximum_discarded_squared_weight::Float64
    full_sweep_changes::Vector{Float64}
    replay_energy_change::Float64
    replay_projector_change::Float64
    projector_envelope::Float64
    state::RHFCompactState
end

@inline function _rhf_energy_envelope(energy::Float64, dimension::Int)
    256eps(Float64) * max(abs(energy), 1.0) * sqrt(max(dimension, 1))
end

function _copy_moving_root!(destination::RHFMovingRoot,
        source::RHFMovingRoot)
    copyto!(@view(destination.covariance[1:source.rank, 1:source.rank]),
        @view(source.covariance[1:source.rank, 1:source.rank]))
    destination.rank = source.rank
    destination.occupied = source.occupied
    destination.total_occupied = source.total_occupied
    destination.closes = source.closes
    destination.opens = source.opens
    destination.turnarounds = source.turnarounds
    destination
end

function _candidate_blocks(candidate, lifecycle::RHFFixedLifecycle,
        workspace::RHFLifecycleWorkspace)
    workspace.empty.provenance = lifecycle.provenance
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

function _candidate_energy!(candidate, lifecycle::RHFFixedLifecycle,
        one_body::BandedOneBody, operator::UnitCellInteraction,
        workspace::RHFLifecycleWorkspace)
    left, right = _candidate_blocks(candidate, lifecycle, workspace)
    root = candidate.root
    expected = left.rank + _center_width(lifecycle, candidate.next_center) +
        right.rank
    root.rank == expected || throw(DimensionMismatch(
        "candidate root does not span the destination center"))
    _prepare_lifecycle_center!(workspace.prepared, left, right,
        candidate.next_center, lifecycle.intervals[candidate.next_center],
        lifecycle.intervals[candidate.next_center + 1], one_body, operator)
    _center_true_energy!(workspace.prepared,
        @view(root.covariance[1:root.rank, 1:root.rank]), root.occupied).total
end

function _trial_advance!(trial::RHFMovingRoot, source::RHFMovingRoot,
        lifecycle::RHFFixedLifecycle, left::RHFOuterBlock,
        right::RHFOuterBlock, covariance, one_body::BandedOneBody,
        operator::UnitCellInteraction, control::RHFStateControl,
        workspace::RHFNonlinearWorkspace)
    _copy_moving_root!(trial, source)
    copyto!(@view(trial.covariance[1:source.rank, 1:source.rank]), covariance)
    candidate = _build_advance_candidate!(trial, lifecycle, left, right,
        one_body, operator, control, workspace.lifecycle,
        workspace.trial_work)
    energy = _candidate_energy!(candidate, lifecycle, one_body, operator,
        workspace.lifecycle)
    workspace.trial_energy_calls += 1
    candidate, energy
end

function _lowest_projector!(output, fock, rank::Int, occupied::Int)
    factors = eigen(Symmetric(Matrix(@view(fock[1:rank, 1:rank]))))
    active = @view output[1:rank, 1:rank]
    fill!(active, 0.0)
    if occupied > 0
        basis = @view factors.vectors[:, 1:occupied]
        mul!(active, basis, transpose(basis))
    end
    active
end

function _projector_basis(projector, rank::Int, occupied::Int)
    occupied == 0 && return zeros(rank, 0)
    factors = eigen(Symmetric(Matrix(@view(projector[1:rank, 1:rank]))))
    Matrix(@view(factors.vectors[:, rank-occupied+1:rank]))
end

function _geodesic_projector!(output, incoming, proposal, rank::Int,
        occupied::Int, step::Float64)
    active = @view output[1:rank, 1:rank]
    occupied == 0 && (fill!(active, 0.0); return active)
    x = _projector_basis(incoming, rank, occupied)
    y = _projector_basis(proposal, rank, occupied)
    factors = svd(transpose(x) * y)
    xa = x * factors.U
    ya = y * transpose(factors.Vt)
    z = similar(xa)
    @inbounds for column = 1:occupied
        cosine = clamp(factors.S[column], 0.0, 1.0)
        angle = acos(cosine)
        sine = sin(angle)
        if sine > 64eps(Float64)
            for row = 1:rank
                tangent = (ya[row, column] - cosine*xa[row, column]) / sine
                z[row, column] = cos(step*angle)*xa[row, column] +
                    sin(step*angle)*tangent
            end
        else
            copyto!(@view(z[:, column]), @view(xa[:, column]))
        end
    end
    basis = Matrix(qr(z).Q)[:, 1:occupied]
    mul!(active, basis, transpose(basis))
    active
end

@inline function _projector_distance(first, second, rank::Int)
    maximum(abs(first[row, column] - second[row, column])
        for column = 1:rank, row = 1:rank; init=0.0)
end

function _preflight_nonlinear(lifecycle, one_body, operator, control,
        workspace)
    lifecycle.provenance == hash((operator.provenance, lifecycle.intervals,
        lifecycle.atom_intervals,
        control.maximum_active, control.cutoff, control.exact)) ||
        throw(ArgumentError("lifecycle controls or operator provenance changed"))
    h1_sites(one_body) == operator.sites || throw(DimensionMismatch(
        "H1 and interaction lengths disagree"))
    lifecycle.root.rank <= size(workspace.incoming, 1) ||
        throw(DimensionMismatch("nonlinear workspace rank capacity is too small"))
    _check_moving_root(lifecycle.root)
    nothing
end

function _monotone_center_transaction!(lifecycle::RHFFixedLifecycle,
        one_body::BandedOneBody, operator::UnitCellInteraction,
        control::RHFStateControl, workspace::RHFNonlinearWorkspace)
    _preflight_nonlinear(lifecycle, one_body, operator, control, workspace)
    workspace.transaction_calls += 1
    localwork = workspace.lifecycle
    left, right = _lifecycle_blocks(lifecycle, localwork)
    cell, forward = lifecycle.next_center, lifecycle.forward
    root = lifecycle.root
    rank = root.rank
    _prepare_lifecycle_center!(localwork.prepared, left, right, cell,
        lifecycle.intervals[cell], lifecycle.intervals[cell + 1], one_body,
        operator)
    incoming_energy = _center_energy_fock!(localwork.prepared,
        @view(root.covariance[1:rank, 1:rank]), root.occupied).total
    copyto!(@view(workspace.incoming[1:rank, 1:rank]),
        @view(root.covariance[1:rank, 1:rank]))
    copyto!(@view(workspace.first_fock[1:rank, 1:rank]),
        @view(localwork.prepared.fock[1:rank, 1:rank]))
    _lowest_projector!(workspace.proposal, workspace.first_fock, rank,
        root.occupied)
    workspace.eigensolution_calls += 1
    self_change = _projector_distance(workspace.incoming,
        workspace.proposal, rank)

    full, proposal_energy = _trial_advance!(workspace.trial_root, root,
        lifecycle, left, right, @view(workspace.proposal[1:rank, 1:rank]),
        one_body, operator, control, workspace)
    envelope = _rhf_energy_envelope(incoming_energy, rank)
    if proposal_energy <= incoming_energy + envelope
        _publish_advance_candidate!(lifecycle, full)
        workspace.publication_calls += 1
        lifecycle.center_visits += 1
        return (cell=cell, forward=forward, incoming_energy=incoming_energy,
            baseline_energy=incoming_energy, proposal_energy=proposal_energy,
            published_energy=proposal_energy, step=1.0, status=:accepted_full,
            self_change=self_change, discarded=full.selection.discarded_squared_weight,
            rank=full.block.rank, center_rank=rank,
            pair_rank=pair_dimension(full.block.rank))
    end

    baseline, baseline_energy = _trial_advance!(workspace.baseline_root, root,
        lifecycle, left, right, @view(workspace.incoming[1:rank, 1:rank]),
        one_body, operator, control, workspace)
    envelope = _rhf_energy_envelope(baseline_energy, rank)
    for step in _RHF_GEODESIC_STEPS
        _geodesic_projector!(workspace.geodesic, workspace.incoming,
            workspace.proposal, rank, root.occupied, step)
        trial, energy = _trial_advance!(workspace.trial_root, root, lifecycle,
            left, right, @view(workspace.geodesic[1:rank, 1:rank]), one_body,
            operator, control, workspace)
        if energy <= baseline_energy + envelope
            _publish_advance_candidate!(lifecycle, trial)
            workspace.publication_calls += 1
            lifecycle.center_visits += 1
            return (cell=cell, forward=forward,
                incoming_energy=incoming_energy, baseline_energy=baseline_energy,
                proposal_energy=proposal_energy, published_energy=energy,
                step=step, status=:accepted_geodesic,
                self_change=self_change,
                discarded=trial.selection.discarded_squared_weight,
                rank=trial.block.rank, center_rank=rank,
                pair_rank=pair_dimension(trial.block.rank))
        end
    end
    _publish_advance_candidate!(lifecycle, baseline)
    workspace.publication_calls += 1
    lifecycle.center_visits += 1
    (cell=cell, forward=forward, incoming_energy=incoming_energy,
        baseline_energy=baseline_energy, proposal_energy=proposal_energy,
        published_energy=baseline_energy, step=0.0,
        status=:rejected_no_safe_trial, self_change=self_change,
        discarded=baseline.selection.discarded_squared_weight,
        rank=baseline.block.rank, center_rank=rank,
        pair_rank=pair_dimension(baseline.block.rank))
end

function _invariant_energy!(lifecycle, one_body, operator, workspace)
    left, right = _lifecycle_blocks(lifecycle, workspace)
    _prepare_lifecycle_center!(workspace.prepared, left, right,
        lifecycle.next_center, lifecycle.intervals[lifecycle.next_center],
        lifecycle.intervals[lifecycle.next_center + 1], one_body, operator)
    root = lifecycle.root
    _center_energy_fock!(workspace.prepared,
        @view(root.covariance[1:root.rank, 1:root.rank]), root.occupied).total
end

function _nonlinear_half_sweep!(lifecycle, one_body, operator, control,
        workspace)
    orientation = lifecycle.forward
    evidence = NamedTuple[]
    while lifecycle.forward == orientation
        push!(evidence, _monotone_center_transaction!(lifecycle, one_body,
            operator, control, workspace))
    end
    lifecycle.half_sweeps += 1
    evidence, _invariant_energy!(lifecycle, one_body, operator,
        workspace.lifecycle)
end

function _solver_loop!(lifecycle, one_body, operator, state_control,
        energy_tolerance::Float64, maximum_half_sweeps::Int,
        workspace::RHFNonlinearWorkspace)
    energies = Float64[]
    changes = Float64[]
    full = fallback = rejected = 0
    maximum_rank = maximum_pair = 0
    maximum_discarded = 0.0
    quiet = 0
    converged = false
    reason = :maximum_half_sweeps
    all_updates = NamedTuple[]
    for half = 1:maximum_half_sweeps
        updates, energy = _nonlinear_half_sweep!(lifecycle, one_body,
            operator, state_control, workspace)
        append!(all_updates, updates); push!(energies, energy)
        for update in updates
            full += update.status === :accepted_full
            fallback += update.status === :accepted_geodesic
            rejected += update.status === :rejected_no_safe_trial
            maximum_rank = max(maximum_rank, update.rank)
            maximum_pair = max(maximum_pair, update.pair_rank)
            maximum_discarded = max(maximum_discarded, update.discarded)
            update.published_energy <= update.baseline_energy +
                _rhf_energy_envelope(update.baseline_energy,
                    update.center_rank) ||
                error("published optimizer update violates monotonicity")
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
    (converged, reason, energies, changes, full, fallback, rejected,
        maximum_rank, maximum_pair, maximum_discarded, all_updates)
end

function solve_hfdmrg(one_body::BandedOneBody,
        operator::UnitCellInteraction; state_maxdim::Int,
        state_cutoff::Real, energy_tolerance::Real,
        maximum_half_sweeps::Int, starting_orientation::Symbol=:physical,
        exact::Bool=false, atom_intervals=nothing, spin::Symbol=:rhf,
        initialization::Symbol=(spin === :uhf ? :atomic_neel : :dimer),
        spin_swapped::Bool=false)
    isfinite(energy_tolerance) && energy_tolerance > 0 || throw(ArgumentError(
        "energy convergence tolerance must be finite and positive"))
    maximum_half_sweeps >= 1 || throw(ArgumentError(
        "maximum half-sweeps must be positive"))
    starting_orientation in (:physical, :reflected) || throw(ArgumentError(
        "starting orientation must be physical or reflected"))
    spin in (:rhf, :uhf) || throw(ArgumentError("spin must be rhf or uhf"))
    if spin === :uhf
        initialization === :atomic_neel || throw(ArgumentError(
            "UHF currently requires atomic_neel initialization"))
        return _solve_uhf_hfdmrg(one_body, operator; state_maxdim,
            state_cutoff, energy_tolerance, maximum_half_sweeps,
            starting_orientation, exact, atom_intervals, spin_swapped)
    end
    initialization === :dimer || throw(ArgumentError(
        "RHF currently requires dimer initialization"))
    !spin_swapped || throw(ArgumentError(
        "spin-swapped initialization applies only to UHF"))
    control = RHFStateControl(state_maxdim, Float64(state_cutoff); exact)
    atoms = atom_intervals === nothing ?
        length(_traversal_intervals(operator)) : length(atom_intervals)
    capacity = exact ? max(state_maxdim, atoms ÷ 2) : state_maxdim
    workspace = RHFNonlinearWorkspace(one_body, operator, capacity)
    lifecycle = initialize_rhf_dimer_lifecycle(one_body, operator, control,
        workspace.lifecycle; forward=starting_orientation === :physical,
        atom_intervals)
    outcome = _solver_loop!(lifecycle, one_body, operator, control,
        Float64(energy_tolerance), maximum_half_sweeps, workspace)

    replay_before = _invariant_energy!(lifecycle, one_body, operator,
        workspace.lifecycle)
    rank_before = lifecycle.root.rank
    copyto!(@view(workspace.incoming[1:rank_before, 1:rank_before]),
        @view(lifecycle.root.covariance[1:rank_before, 1:rank_before]))
    for _ = 1:2
        run_fixed_half_sweep!(lifecycle, one_body, operator, control,
            workspace.lifecycle)
    end
    replay_after = _invariant_energy!(lifecycle, one_body, operator,
        workspace.lifecycle)
    replay_projector = lifecycle.root.rank == rank_before ?
        _projector_distance(workspace.incoming, lifecycle.root.covariance,
            rank_before) : Inf
    projector_envelope = 4sqrt(outcome[10]) +
        8192eps(Float64)*max(lifecycle.root.rank, 1)
    converged = outcome[1]
    reason = outcome[2]
    replay_energy = replay_after - replay_before
    if exact
        envelope = _rhf_energy_envelope(replay_before, lifecycle.root.rank)
        # Root-coordinate matrices before and after a complete replay live in
        # independently rotated outer bases.  Exact physical projector parity
        # is certified by the reversible links and is checked by the small
        # streamed reconstruction oracle, not by this coordinate diagnostic.
        if converged && abs(replay_energy) > envelope
            converged = false; reason = :exact_replay_failed
        end
    elseif converged && abs(replay_energy) > Float64(energy_tolerance) +
            _rhf_energy_envelope(replay_before, lifecycle.root.rank)
        converged = false; reason = :finite_replay_floor
    end
    HFDMRGResult(converged, reason, replay_after, length(outcome[3]),
        outcome[5], outcome[6], outcome[7], outcome[8], outcome[9],
        outcome[10], outcome[4], replay_energy, replay_projector,
        projector_envelope, RHFCompactState(lifecycle))
end
