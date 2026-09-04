struct UHFStateControl
    alpha::RHFStateControl
    beta::RHFStateControl
end

UHFStateControl(control::RHFStateControl) = UHFStateControl(control, control)

mutable struct UHFMovingRoot
    alpha::RHFMovingRoot
    beta::RHFMovingRoot
end

function _check_uhf_root(root::UHFMovingRoot)
    _check_moving_root(root.alpha); _check_moving_root(root.beta)
    nothing
end

mutable struct UHFFixedLifecycle
    collection::Vector{UHFOuterBlock}
    root::UHFMovingRoot
    intervals::Vector{UnitRange{Int}}
    atom_intervals::Vector{UnitRange{Int}}
    cells::Int
    forward::Bool
    next_center::Int
    provenance::UInt
    generation::Int
    half_sweeps::Int
    fragment_ingestions::Int
    original_rows_read::Int
    replacements::Int
    center_visits::Int
end

mutable struct UHFLifecycleWorkspace
    prepared::FastUHFPreparedCenter
    builder::UHFBlockBuildWorkspace
    alpha_root::RHFRootWorkspace
    beta_root::RHFRootWorkspace
    empty::UHFOuterBlock
    maximum_alpha::Int
    maximum_beta::Int
end

function UHFLifecycleWorkspace(one_body::BandedOneBody,
        operator::UnitCellInteraction, maximum_alpha::Int,
        maximum_beta::Int; maximum_rhs::Int=3)
    min(maximum_alpha, maximum_beta) >= 0 || throw(ArgumentError(
        "UHF state capacities must be nonnegative"))
    center_alpha = 2UNIT_CELL_WIDTH + 2maximum_alpha
    center_beta = 2UNIT_CELL_WIDTH + 2maximum_beta
    channels = _maximum_channel_rank(operator)
    prepared = FastUHFPreparedCenter(center_alpha, center_beta, operator;
        maximum_rhs=max(maximum_rhs, channels),
        maximum_alpha_side=maximum_alpha,
        maximum_beta_side=maximum_beta)
    builder = UHFBlockBuildWorkspace(max(UNIT_CELL_WIDTH + maximum_alpha,
            UNIT_CELL_WIDTH), max(UNIT_CELL_WIDTH + maximum_beta,
            UNIT_CELL_WIDTH), channels)
    provenance = UInt(0)
    UHFLifecycleWorkspace(prepared, builder,
        RHFRootWorkspace(Matrix{Float64}(undef, center_alpha, center_alpha),
            Matrix{Float64}(undef, center_alpha, center_alpha)),
        RHFRootWorkspace(Matrix{Float64}(undef, center_beta, center_beta),
            Matrix{Float64}(undef, center_beta, center_beta)),
        _empty_uhf_block(operator, provenance), center_alpha, center_beta)
end

function _uhf_neel_fragments(one_body::BandedOneBody, atom_intervals)
    atoms = Vector{UnitRange{Int}}(atom_intervals)
    !isempty(atoms) || throw(ArgumentError("atomic product requires atoms"))
    previous = 0
    for interval in atoms
        !isempty(interval) && first(interval) == previous + 1 ||
            throw(ArgumentError("atom intervals must be ordered and contiguous"))
        previous = last(interval)
    end
    previous == h1_sites(one_body) || throw(ArgumentError(
        "atom intervals must partition the complete one-body range"))
    alpha = _RHFLocalFragment[]; beta = _RHFLocalFragment[]
    for (index, rows) in enumerate(atoms)
        fragment = _RHFLocalFragment(rows, _dimer_orbital(one_body, rows))
        push!(isodd(index) ? alpha : beta, fragment)
    end
    alpha, beta
end

function _uhf_physical_block(one_body::BandedOneBody,
        operator::UnitCellInteraction, cell::Int, rows::UnitRange{Int},
        provenance::UInt; reflected::Bool=false)
    seed_uhf_interval_block(one_body, operator, cell, rows; provenance,
        reflected)
end

function initialize_uhf_neel_lifecycle(one_body::BandedOneBody,
        operator::UnitCellInteraction, control::UHFStateControl,
        workspace::UHFLifecycleWorkspace; forward::Bool=true,
        atom_intervals=nothing, spin_swapped::Bool=false)
    h1_sites(one_body) == operator.sites || throw(DimensionMismatch(
        "H1 and interaction lengths disagree"))
    intervals = _check_traversal_intervals(_traversal_intervals(operator),
        operator.sites)
    atoms = atom_intervals === nothing ? intervals :
        Vector{UnitRange{Int}}(atom_intervals)
    alpha_fragments, beta_fragments = _uhf_neel_fragments(one_body, atoms)
    spin_swapped && ((alpha_fragments, beta_fragments) =
        (beta_fragments, alpha_fragments))
    alpha_cursor = _RHFFragmentCursor(alpha_fragments)
    beta_cursor = _RHFFragmentCursor(beta_fragments)
    cells = length(intervals)
    provenance = hash((operator.provenance, intervals, atoms,
        control.alpha.maximum_active, control.alpha.cutoff,
        control.alpha.exact, control.beta.maximum_active,
        control.beta.cutoff, control.beta.exact, spin_swapped))
    empty = _empty_uhf_block(operator, provenance)
    collection = Vector{UHFOuterBlock}(undef, cells - 1)
    current = empty
    if forward
        collection[1] = empty
        for cell = cells:-1:3
            physical = _uhf_physical_block(one_body, operator, cell,
                intervals[cell], provenance)
            alpha_ids, alpha_parent, alpha_cross = _fragment_parent(
                alpha_cursor, intervals[cell], current.alpha.rank, true)
            beta_ids, beta_parent, beta_cross = _fragment_parent(beta_cursor,
                intervals[cell], current.beta.rank, true)
            alpha_link = _factor_state_link(alpha_parent *
                transpose(alpha_parent), length(intervals[cell]),
                current.alpha.rank, control.alpha; cross=alpha_cross,
                maximum_occupied=length(alpha_ids))
            beta_link = _factor_state_link(beta_parent * transpose(beta_parent),
                length(intervals[cell]), current.beta.rank, control.beta;
                cross=beta_cross, maximum_occupied=length(beta_ids))
            current = build_uhf_outer_block(physical, current, alpha_link,
                beta_link, one_body, operator, workspace.builder)
            collection[cell-1] = current
            _advance_fragment_cursor!(alpha_cursor, alpha_ids, alpha_parent,
                alpha_link)
            _advance_fragment_cursor!(beta_cursor, beta_ids, beta_parent,
                beta_link)
        end
    else
        collection[cells-1] = empty
        for cell = 1:cells-2
            physical = _uhf_physical_block(one_body, operator, cell,
                intervals[cell], provenance; reflected=true)
            alpha_ids, alpha_parent, alpha_cross = _fragment_parent(
                alpha_cursor, intervals[cell], current.alpha.rank, false)
            beta_ids, beta_parent, beta_cross = _fragment_parent(beta_cursor,
                intervals[cell], current.beta.rank, false)
            alpha_link = _factor_state_link(alpha_parent *
                transpose(alpha_parent), current.alpha.rank,
                length(intervals[cell]), control.alpha; cross=alpha_cross,
                maximum_occupied=length(alpha_ids))
            beta_link = _factor_state_link(beta_parent * transpose(beta_parent),
                current.beta.rank, length(intervals[cell]), control.beta;
                cross=beta_cross, maximum_occupied=length(beta_ids))
            current = build_uhf_outer_block(current, physical, alpha_link,
                beta_link, one_body, operator, workspace.builder)
            current.alpha.reflected = true; current.beta.reflected = true
            collection[cell] = current
            _advance_fragment_cursor!(alpha_cursor, alpha_ids, alpha_parent,
                alpha_link)
            _advance_fragment_cursor!(beta_cursor, beta_ids, beta_parent,
                beta_link)
        end
    end
    center_rows = first(intervals[forward ? 1 : cells-1]):last(
        intervals[forward ? 2 : cells])
    alpha_projector, alpha_occupied = _root_fragment_projector(alpha_cursor,
        center_rows, current.alpha.rank, forward, length(alpha_fragments))
    beta_projector, beta_occupied = _root_fragment_projector(beta_cursor,
        center_rows, current.beta.rank, forward, length(beta_fragments))
    alpha_covariance = fill(NaN, workspace.maximum_alpha,
        workspace.maximum_alpha)
    beta_covariance = fill(NaN, workspace.maximum_beta,
        workspace.maximum_beta)
    copyto!(@view(alpha_covariance[1:size(alpha_projector, 1),
        1:size(alpha_projector, 2)]), alpha_projector)
    copyto!(@view(beta_covariance[1:size(beta_projector, 1),
        1:size(beta_projector, 2)]), beta_projector)
    root = UHFMovingRoot(RHFMovingRoot(alpha_covariance,
            size(alpha_projector, 1), alpha_occupied,
            length(alpha_fragments), 0, 0, 0),
        RHFMovingRoot(beta_covariance, size(beta_projector, 1), beta_occupied,
            length(beta_fragments), 0, 0, 0))
    lifecycle = UHFFixedLifecycle(collection, root, intervals, atoms, cells,
        forward, forward ? 1 : cells-1, provenance, 0, 0,
        alpha_cursor.ingestions + beta_cursor.ingestions, 0, 0, 0)
    _check_uhf_root(root)
    lifecycle
end

function _uhf_lifecycle_blocks(lifecycle::UHFFixedLifecycle,
        workspace::UHFLifecycleWorkspace)
    cell = lifecycle.next_center
    workspace.empty.provenance = lifecycle.provenance
    workspace.empty.alpha.provenance = lifecycle.provenance
    workspace.empty.beta.provenance = lifecycle.provenance
    left = cell == 1 ? workspace.empty : lifecycle.collection[cell-1]
    right = cell == lifecycle.cells-1 ? workspace.empty :
        lifecycle.collection[cell+1]
    expected_alpha = left.alpha.rank + _center_width_uhf(lifecycle, cell) +
        right.alpha.rank
    expected_beta = left.beta.rank + _center_width_uhf(lifecycle, cell) +
        right.beta.rank
    lifecycle.root.alpha.rank == expected_alpha &&
        lifecycle.root.beta.rank == expected_beta || throw(DimensionMismatch(
        "UHF collection and root ranks disagree"))
    lifecycle.root.alpha.occupied + left.alpha.completed_count +
        right.alpha.completed_count == lifecycle.root.alpha.total_occupied ||
        throw(ArgumentError("alpha occupied count is inconsistent"))
    lifecycle.root.beta.occupied + left.beta.completed_count +
        right.beta.completed_count == lifecycle.root.beta.total_occupied ||
        throw(ArgumentError("beta occupied count is inconsistent"))
    left, right
end

@inline _center_width_uhf(lifecycle::UHFFixedLifecycle, cell::Int) =
    length(lifecycle.intervals[cell]) + length(lifecycle.intervals[cell+1])

function _restore_subroot!(root::RHFMovingRoot, backup, rank, occupied,
        closes, opens, turnarounds)
    _restore_root!(root, backup, rank, occupied, closes, opens, turnarounds)
end

function _uhf_build_candidate!(root::UHFMovingRoot,
        lifecycle::UHFFixedLifecycle,
        left::UHFOuterBlock, right::UHFOuterBlock,
        one_body::BandedOneBody, operator::UnitCellInteraction,
        control::UHFStateControl, workspace::UHFLifecycleWorkspace,
        alpha_workspace::RHFRootWorkspace=workspace.alpha_root,
        beta_workspace::RHFRootWorkspace=workspace.beta_root)
    cell = lifecycle.next_center; forward = lifecycle.forward
    terminal = forward ? cell == lifecycle.cells-1 : cell == 1
    alpha_root, beta_root = root.alpha, root.beta
    if forward && !terminal
        width = length(lifecycle.intervals[cell])
        alpha_local = left.alpha.rank + width
        beta_local = left.beta.rank + width
        alpha_link = _factor_state_link(
            @view(alpha_root.covariance[1:alpha_local, 1:alpha_local]),
            left.alpha.rank, width, control.alpha;
            cross=@view(alpha_root.covariance[1:alpha_local,
                alpha_local+1:alpha_root.rank]),
            maximum_occupied=alpha_root.occupied)
        beta_link = _factor_state_link(
            @view(beta_root.covariance[1:beta_local, 1:beta_local]),
            left.beta.rank, width, control.beta;
            cross=@view(beta_root.covariance[1:beta_local,
                beta_local+1:beta_root.rank]),
            maximum_occupied=beta_root.occupied)
        physical = _uhf_physical_block(one_body, operator, cell,
            lifecycle.intervals[cell], lifecycle.provenance)
        block = build_uhf_outer_block(left, physical, alpha_link, beta_link,
            one_body, operator, workspace.builder)
        _close_root!(alpha_root, alpha_link, alpha_local, true,
            alpha_workspace; exact=control.alpha.exact)
        _close_root!(beta_root, beta_link, beta_local, true,
            beta_workspace; exact=control.beta.exact)
        _open_root!(alpha_root, right.alpha.link, false, alpha_workspace)
        _open_root!(beta_root, right.beta.link, false, beta_workspace)
        next_forward, next_center = true, cell + 1
    elseif forward
        width = length(lifecycle.intervals[cell+1])
        alpha_local = width + right.alpha.rank
        beta_local = width + right.beta.rank
        alpha_start = alpha_root.rank-alpha_local+1
        beta_start = beta_root.rank-beta_local+1
        alpha_link = _factor_state_link(@view(alpha_root.covariance[
                alpha_start:alpha_root.rank, alpha_start:alpha_root.rank]),
            width, right.alpha.rank, control.alpha;
            cross=transpose(@view(alpha_root.covariance[1:alpha_start-1,
                alpha_start:alpha_root.rank])),
            maximum_occupied=alpha_root.occupied)
        beta_link = _factor_state_link(@view(beta_root.covariance[
                beta_start:beta_root.rank, beta_start:beta_root.rank]),
            width, right.beta.rank, control.beta;
            cross=transpose(@view(beta_root.covariance[1:beta_start-1,
                beta_start:beta_root.rank])),
            maximum_occupied=beta_root.occupied)
        physical = _uhf_physical_block(one_body, operator, cell+1,
            lifecycle.intervals[cell+1], lifecycle.provenance; reflected=true)
        block = build_uhf_outer_block(physical, right, alpha_link, beta_link,
            one_body, operator, workspace.builder)
        block.alpha.reflected = true; block.beta.reflected = true
        _close_root!(alpha_root, alpha_link, alpha_local, false,
            alpha_workspace; exact=control.alpha.exact)
        _close_root!(beta_root, beta_link, beta_local, false,
            beta_workspace; exact=control.beta.exact)
        if lifecycle.cells == 2
            _open_root!(alpha_root, block.alpha.link, false,
                alpha_workspace)
            _open_root!(beta_root, block.beta.link, false, beta_workspace)
        else
            _open_root!(alpha_root, left.alpha.link, true,
                alpha_workspace)
            _open_root!(beta_root, left.beta.link, true, beta_workspace)
        end
        next_forward, next_center = false,
            lifecycle.cells == 2 ? 1 : cell - 1
        alpha_root.turnarounds += 1; beta_root.turnarounds += 1
    elseif !terminal
        width = length(lifecycle.intervals[cell+1])
        alpha_local = width + right.alpha.rank
        beta_local = width + right.beta.rank
        alpha_start = alpha_root.rank-alpha_local+1
        beta_start = beta_root.rank-beta_local+1
        alpha_link = _factor_state_link(@view(alpha_root.covariance[
                alpha_start:alpha_root.rank, alpha_start:alpha_root.rank]),
            width, right.alpha.rank, control.alpha;
            cross=transpose(@view(alpha_root.covariance[1:alpha_start-1,
                alpha_start:alpha_root.rank])),
            maximum_occupied=alpha_root.occupied)
        beta_link = _factor_state_link(@view(beta_root.covariance[
                beta_start:beta_root.rank, beta_start:beta_root.rank]),
            width, right.beta.rank, control.beta;
            cross=transpose(@view(beta_root.covariance[1:beta_start-1,
                beta_start:beta_root.rank])),
            maximum_occupied=beta_root.occupied)
        physical = _uhf_physical_block(one_body, operator, cell+1,
            lifecycle.intervals[cell+1], lifecycle.provenance; reflected=true)
        block = build_uhf_outer_block(physical, right, alpha_link, beta_link,
            one_body, operator, workspace.builder)
        block.alpha.reflected = true; block.beta.reflected = true
        _close_root!(alpha_root, alpha_link, alpha_local, false,
            alpha_workspace; exact=control.alpha.exact)
        _close_root!(beta_root, beta_link, beta_local, false,
            beta_workspace; exact=control.beta.exact)
        _open_root!(alpha_root, left.alpha.link, true, alpha_workspace)
        _open_root!(beta_root, left.beta.link, true, beta_workspace)
        next_forward, next_center = false, cell - 1
    else
        width = length(lifecycle.intervals[cell])
        alpha_local = left.alpha.rank + width
        beta_local = left.beta.rank + width
        alpha_link = _factor_state_link(
            @view(alpha_root.covariance[1:alpha_local, 1:alpha_local]),
            left.alpha.rank, width, control.alpha;
            cross=@view(alpha_root.covariance[1:alpha_local,
                alpha_local+1:alpha_root.rank]),
            maximum_occupied=alpha_root.occupied)
        beta_link = _factor_state_link(
            @view(beta_root.covariance[1:beta_local, 1:beta_local]),
            left.beta.rank, width, control.beta;
            cross=@view(beta_root.covariance[1:beta_local,
                beta_local+1:beta_root.rank]),
            maximum_occupied=beta_root.occupied)
        physical = _uhf_physical_block(one_body, operator, cell,
            lifecycle.intervals[cell], lifecycle.provenance)
        block = build_uhf_outer_block(left, physical, alpha_link, beta_link,
            one_body, operator, workspace.builder)
        _close_root!(alpha_root, alpha_link, alpha_local, true,
            alpha_workspace; exact=control.alpha.exact)
        _close_root!(beta_root, beta_link, beta_local, true,
            beta_workspace; exact=control.beta.exact)
        if lifecycle.cells == 2
            _open_root!(alpha_root, block.alpha.link, true,
                alpha_workspace)
            _open_root!(beta_root, block.beta.link, true, beta_workspace)
        else
            _open_root!(alpha_root, right.alpha.link, false,
                alpha_workspace)
            _open_root!(beta_root, right.beta.link, false, beta_workspace)
        end
        next_forward, next_center = true,
            lifecycle.cells == 2 ? 1 : cell + 1
        alpha_root.turnarounds += 1; beta_root.turnarounds += 1
    end
    (block=block, root=root, cell=cell, forward=next_forward,
        next_center=next_center,
        turned=next_forward != forward, alpha_selection=alpha_link.selection,
        beta_selection=beta_link.selection,
        published_bytes=uhf_block_storage_bytes(block))
end

function _publish_uhf_advance_candidate!(lifecycle::UHFFixedLifecycle,
        candidate)
    if candidate.root !== lifecycle.root
        _copy_moving_root!(lifecycle.root.alpha, candidate.root.alpha)
        _copy_moving_root!(lifecycle.root.beta, candidate.root.beta)
    end
    lifecycle.collection[candidate.cell] = candidate.block
    lifecycle.forward = candidate.forward
    lifecycle.next_center = candidate.next_center
    candidate.turned && (lifecycle.generation += 1)
    lifecycle.replacements += 1
    candidate
end

function _transactional_uhf_advance!(lifecycle::UHFFixedLifecycle,
        left::UHFOuterBlock, right::UHFOuterBlock,
        one_body::BandedOneBody, operator::UnitCellInteraction,
        control::UHFStateControl, workspace::UHFLifecycleWorkspace)
    alpha, beta = lifecycle.root.alpha, lifecycle.root.beta
    alpha_state = (alpha.rank, alpha.occupied, alpha.closes, alpha.opens,
        alpha.turnarounds)
    beta_state = (beta.rank, beta.occupied, beta.closes, beta.opens,
        beta.turnarounds)
    copyto!(@view(workspace.alpha_root.second[1:alpha.rank, 1:alpha.rank]),
        @view(alpha.covariance[1:alpha.rank, 1:alpha.rank]))
    copyto!(@view(workspace.beta_root.second[1:beta.rank, 1:beta.rank]),
        @view(beta.covariance[1:beta.rank, 1:beta.rank]))
    candidate = try
        _uhf_build_candidate!(lifecycle.root, lifecycle, left, right,
            one_body, operator, control, workspace)
    catch
        _restore_subroot!(alpha, workspace.alpha_root.second, alpha_state...)
        _restore_subroot!(beta, workspace.beta_root.second, beta_state...)
        rethrow()
    end
    _publish_uhf_advance_candidate!(lifecycle, candidate)
end

function run_uhf_fixed_half_sweep!(lifecycle::UHFFixedLifecycle,
        one_body::BandedOneBody, operator::UnitCellInteraction,
        control::UHFStateControl, workspace::UHFLifecycleWorkspace)
    _check_uhf_root(lifecycle.root)
    lifecycle.provenance == hash((operator.provenance, lifecycle.intervals,
        lifecycle.atom_intervals, control.alpha.maximum_active,
        control.alpha.cutoff, control.alpha.exact,
        control.beta.maximum_active, control.beta.cutoff, control.beta.exact,
        false)) || lifecycle.provenance == hash((operator.provenance,
        lifecycle.intervals, lifecycle.atom_intervals,
        control.alpha.maximum_active, control.alpha.cutoff,
        control.alpha.exact, control.beta.maximum_active,
        control.beta.cutoff, control.beta.exact, true)) || throw(ArgumentError(
        "UHF lifecycle controls or operator provenance changed"))
    starting = lifecycle.forward; evidence = NamedTuple[]
    while lifecycle.forward == starting
        left, right = _uhf_lifecycle_blocks(lifecycle, workspace)
        cell = lifecycle.next_center
        prepare_uhf_lifecycle_center!(workspace.prepared, left, right, cell,
            lifecycle.intervals[cell], lifecycle.intervals[cell+1], one_body,
            operator)
        energy = uhf_center_energy_fock!(workspace.prepared,
            @view(lifecycle.root.alpha.covariance[1:lifecycle.root.alpha.rank,
                1:lifecycle.root.alpha.rank]),
            @view(lifecycle.root.beta.covariance[1:lifecycle.root.beta.rank,
                1:lifecycle.root.beta.rank]), lifecycle.root.alpha.occupied,
            lifecycle.root.beta.occupied)
        candidate = _transactional_uhf_advance!(lifecycle, left, right,
            one_body, operator, control, workspace)
        lifecycle.center_visits += 1
        push!(evidence, (cell=cell, forward=starting, energy=energy.total,
            alpha_rank=workspace.prepared.alpha.active_rank,
            beta_rank=workspace.prepared.beta.active_rank,
            alpha_active=candidate.alpha_selection.active,
            beta_active=candidate.beta_selection.active,
            alpha_discarded=candidate.alpha_selection.discarded_squared_weight,
            beta_discarded=candidate.beta_selection.discarded_squared_weight,
            cross_rank=cross_hartree_rank(candidate.block.cross_hartree),
            published_bytes=candidate.published_bytes))
    end
    lifecycle.half_sweeps += 1
    evidence
end

@inline function _uhf_spin_block(lifecycle::UHFFixedLifecycle, slot::Int,
        alpha::Bool)
    block = lifecycle.collection[slot]
    alpha ? block.alpha : block.beta
end

function _expand_uhf_prefix!(output, column::Int, initial, slot::Int,
        lifecycle::UHFFixedLifecycle, workspace::RHFRootWorkspace,
        alpha::Bool)
    current = initial
    first_buffer = @view workspace.first[:, 1]
    second_buffer = @view workspace.second[:, 1]
    for index = slot:-1:1
        block = _uhf_spin_block(lifecycle, index, alpha)
        parent = block.link.old_rank + block.link.physical_width
        mul!(@view(first_buffer[1:parent]), block.link.close_map, current)
        rows = lifecycle.intervals[index]
        copyto!(@view(output[rows, column]),
            @view(first_buffer[block.link.old_rank+1:parent]))
        copyto!(@view(second_buffer[1:block.link.old_rank]),
            @view(first_buffer[1:block.link.old_rank]))
        current = @view second_buffer[1:block.link.old_rank]
    end
    nothing
end

function _expand_uhf_suffix!(output, column::Int, initial, slot::Int,
        lifecycle::UHFFixedLifecycle, workspace::RHFRootWorkspace,
        alpha::Bool)
    current = initial
    first_buffer = @view workspace.first[:, 1]
    second_buffer = @view workspace.second[:, 1]
    for index = slot:lifecycle.cells-1
        block = _uhf_spin_block(lifecycle, index, alpha)
        parent = block.link.old_rank + block.link.physical_width
        mul!(@view(first_buffer[1:parent]), block.link.close_map, current)
        rows = lifecycle.intervals[index+1]; width = length(rows)
        copyto!(@view(output[rows, column]), @view(first_buffer[1:width]))
        child = block.link.physical_width
        copyto!(@view(second_buffer[1:child]),
            @view(first_buffer[width+1:parent]))
        current = @view second_buffer[1:child]
    end
    nothing
end

function _reconstruct_uhf_spin!(output, lifecycle::UHFFixedLifecycle,
        workspace::RHFRootWorkspace, alpha::Bool)
    root = alpha ? lifecycle.root.alpha : lifecycle.root.beta
    size(output) == (last(lifecycle.intervals[end]), root.total_occupied) ||
        throw(DimensionMismatch("UHF reconstruction output dimensions disagree"))
    cell = lifecycle.next_center; left_slot = cell-1; right_slot = cell+1
    left_rank = left_slot == 0 ? 0 :
        _uhf_spin_block(lifecycle, left_slot, alpha).rank
    right_rank = right_slot >= lifecycle.cells ? 0 :
        _uhf_spin_block(lifecycle, right_slot, alpha).rank
    factors = eigen(Symmetric(Matrix(@view(root.covariance[1:root.rank,
        1:root.rank]))))
    fill!(output, 0.0); column = 0
    for mode = root.rank-root.occupied+1:root.rank
        column += 1; vector = @view factors.vectors[:, mode]
        left_rank > 0 && _expand_uhf_prefix!(output, column,
            @view(vector[1:left_rank]), left_slot, lifecycle, workspace, alpha)
        rows = first(lifecycle.intervals[cell]):last(lifecycle.intervals[cell+1])
        copyto!(@view(output[rows, column]),
            @view(vector[left_rank+1:left_rank+length(rows)]))
        right_rank > 0 && _expand_uhf_suffix!(output, column,
            @view(vector[end-right_rank+1:end]), right_slot, lifecycle,
            workspace, alpha)
    end
    for slot = 1:left_slot
        link = _uhf_spin_block(lifecycle, slot, alpha).link
        for mode in axes(link.completed_map, 2)
            column += 1; rows = lifecycle.intervals[slot]
            copyto!(@view(output[rows, column]),
                @view(link.completed_map[link.old_rank+1:end, mode]))
            link.old_rank > 0 && _expand_uhf_prefix!(output, column,
                @view(link.completed_map[1:link.old_rank, mode]), slot-1,
                lifecycle, workspace, alpha)
        end
    end
    for slot = right_slot:lifecycle.cells-1
        link = _uhf_spin_block(lifecycle, slot, alpha).link
        for mode in axes(link.completed_map, 2)
            column += 1; rows = lifecycle.intervals[slot+1]
            width = length(rows)
            copyto!(@view(output[rows, column]),
                @view(link.completed_map[1:width, mode]))
            link.physical_width > 0 && _expand_uhf_suffix!(output, column,
                @view(link.completed_map[width+1:end, mode]), slot+1,
                lifecycle, workspace, alpha)
        end
    end
    column == root.total_occupied || throw(ArgumentError(
        "UHF reconstruction omitted an occupied spin orbital"))
    output
end

function reconstruct_uhf_terminal!(alpha_output::AbstractMatrix{Float64},
        beta_output::AbstractMatrix{Float64},
        lifecycle::UHFFixedLifecycle, workspace::UHFLifecycleWorkspace)
    Base.mightalias(alpha_output, beta_output) && throw(ArgumentError(
        "spin reconstruction outputs must not alias"))
    _reconstruct_uhf_spin!(alpha_output, lifecycle, workspace.alpha_root, true)
    _reconstruct_uhf_spin!(beta_output, lifecycle, workspace.beta_root, false)
    alpha_output, beta_output
end

function spin_swap_uhf_block(block::UHFOuterBlock,
        operator::UnitCellInteraction)
    alpha = deepcopy(block.beta); beta = deepcopy(block.alpha)
    left_core = alpha.left_core; right_core = alpha.right_core
    beta.left_core = left_core; beta.right_core = right_core
    record = UHFCrossHartree(copy(block.cross_hartree.beta),
        copy(block.cross_hartree.alpha), copy(block.cross_hartree.weight),
        alpha.rank, beta.rank)
    output = UHFOuterBlock(alpha, beta, record, block.provenance)
    _validate_uhf_block(output, operator)
    output
end

function _spin_swap_uhf_lifecycle!(lifecycle::UHFFixedLifecycle,
        operator::UnitCellInteraction)
    _check_uhf_root(lifecycle.root)
    for block in lifecycle.collection
        _validate_uhf_block(block, operator)
    end
    for block in lifecycle.collection
        block.alpha, block.beta = block.beta, block.alpha
        record = block.cross_hartree
        block.cross_hartree = UHFCrossHartree(record.beta, record.alpha,
            record.weight, record.beta_rank, record.alpha_rank)
    end
    lifecycle.root.alpha, lifecycle.root.beta = lifecycle.root.beta,
        lifecycle.root.alpha
    _check_uhf_root(lifecycle.root)
    lifecycle
end

function uhf_lifecycle_storage_bytes(lifecycle::UHFFixedLifecycle,
        operator::UnitCellInteraction, workspace::UHFLifecycleWorkspace)
    (collection=sum(uhf_block_storage_bytes, lifecycle.collection),
        root=sizeof(lifecycle.root.alpha.covariance) +
            sizeof(lifecycle.root.beta.covariance),
        static_operator=sizeof(operator.tail_transfer) +
            sizeof(operator.residual_transfer) + sizeof(operator.phase_source) +
            sizeof(operator.phase_sink) + sizeof(operator.local_triangle),
        reusable_workspace=Base.summarysize(workspace))
end
