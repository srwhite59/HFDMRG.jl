mutable struct RHFOuterBlock
    cells::UnitRange{Int}
    rows::UnitRange{Int}
    rank::Int
    completed_count::Int
    constant_energy::Float64
    h1::Matrix{Float64}
    pair_field::Matrix{Float64}
    left_feature::Matrix{Float64}
    right_feature::Matrix{Float64}
    left_core::Vector{Float64}
    right_core::Vector{Float64}
    left_edge::Matrix{Float64}
    right_edge::Matrix{Float64}
    link::RHFStateLink
    reflected::Bool
    provenance::UInt
end


function _prepare_physical_block!(block::RHFOuterBlock, cell::Int,
        rows::UnitRange{Int},
        one_body::BandedOneBody, operator::UnitCellInteraction,
        provenance::UInt; reflected::Bool=false)
    width = length(rows)
    block.rank == width && block.completed_count == 0 ||
        throw(ArgumentError("physical block workspace is malformed"))
    1 <= first(rows) <= last(rows) <= operator.sites || throw(BoundsError(
        1:operator.sites, rows))
    @inbounds for column = 1:width, row = 1:width
        block.h1[row, column] = h1_entry(one_body, first(rows) + row - 1,
            first(rows) + column - 1)
    end
    block.cells = cell:cell
    block.rows = rows
    block.reflected = reflected
    block.provenance = provenance
    block
end

function _prepare_physical_block!(block::RHFOuterBlock, cell::Int,
        one_body::BandedOneBody, operator::UnitCellInteraction,
        provenance::UInt; reflected::Bool=false)
    rows = (cell - 1) * UNIT_CELL_WIDTH + 1:cell * UNIT_CELL_WIDTH
    _prepare_physical_block!(block, cell, rows, one_body, operator,
        provenance; reflected)
end

mutable struct BlockBuildWorkspace
    maximum_rank::Int
    maximum_pairs::Int
    channel_rank::Int
    h1::Matrix{Float64}
    pair_field::Matrix{Float64}
    feature_left::Matrix{Float64}
    feature_right::Matrix{Float64}
    pair_map::Matrix{Float64}
    pair_scratch::Matrix{Float64}
    h1_scratch::Matrix{Float64}
    transfer::Matrix{Float64}
    transfer_scratch::Matrix{Float64}
    left_core::Vector{Float64}
    right_core::Vector{Float64}
    channel_transfer_calls::Int
    residual_transfer_steps::Int
    maximum_transfer_power::Int
end

function BlockBuildWorkspace(maximum_rank::Int, channel_rank::Int)
    maximum_rank >= 1 || throw(ArgumentError("maximum rank must be positive"))
    channel_rank >= 0 || throw(ArgumentError("channel rank must be nonnegative"))
    pairs = pair_dimension(maximum_rank)
    BlockBuildWorkspace(maximum_rank, pairs, channel_rank,
        zeros(maximum_rank, maximum_rank), zeros(pairs, pairs),
        zeros(pairs, channel_rank), zeros(pairs, channel_rank),
        zeros(pairs, pairs), zeros(pairs, pairs),
        zeros(maximum_rank, maximum_rank), zeros(channel_rank, channel_rank),
        zeros(channel_rank, channel_rank), zeros(channel_rank),
        zeros(channel_rank), 0, 0, 0)
end

@inline function _record_channel_transfer!(workspace::BlockBuildWorkspace,
        operator::UnitCellInteraction, power::Int)
    workspace.channel_transfer_calls += 1
    size(operator.residual_transfer, 1) > 0 &&
        (workspace.residual_transfer_steps += power)
    workspace.maximum_transfer_power = max(workspace.maximum_transfer_power,
        power)
    nothing
end

function _channel_transfer!(output, scratch, operator::UnitCellInteraction,
        power::Int)
    rank = _maximum_channel_rank(operator)
    size(output, 1) >= rank && size(output, 2) >= rank ||
        throw(DimensionMismatch("transfer capacity is too small"))
    fill!(@view(output[1:rank, 1:rank]), 0.0)
    tail = length(operator.tail_transfer)
    @inbounds for channel = 1:tail
        output[channel, channel] = operator.tail_transfer[channel]^power
    end
    residual = size(operator.residual_transfer, 1)
    if residual > 0
        block = @view output[tail+1:tail+residual, tail+1:tail+residual]
        work = @view scratch[tail+1:tail+residual, tail+1:tail+residual]
        fill!(block, 0.0)
        @inbounds for channel = 1:residual
            block[channel, channel] = 1.0
        end
        for _ = 1:power
            mul!(work, block, operator.residual_transfer)
            copyto!(block, work)
        end
    end
    @view output[1:rank, 1:rank]
end

function _validate_block(block::RHFOuterBlock, operator::UnitCellInteraction)
    rank = block.rank
    pairs = pair_dimension(rank)
    channels = _maximum_channel_rank(operator)
    block.completed_count >= 0 || throw(ArgumentError(
        "completed count must be nonnegative"))
    isfinite(block.constant_energy) || throw(ArgumentError(
        "block constant energy must be finite"))
    size(block.h1) == (rank, rank) || throw(DimensionMismatch("H1 block size"))
    size(block.pair_field) == (pairs, pairs) || throw(DimensionMismatch(
        "pair block size"))
    size(block.left_feature) == (pairs, channels) || throw(DimensionMismatch(
        "left feature size"))
    size(block.right_feature) == (pairs, channels) || throw(DimensionMismatch(
        "right feature size"))
    length(block.left_core) == channels && length(block.right_core) == channels ||
        throw(DimensionMismatch("completed-core channel size"))
    size(block.left_edge, 2) == rank && size(block.right_edge, 2) == rank ||
        throw(DimensionMismatch("edge projection size"))
    block.link.active_rank == rank || throw(DimensionMismatch(
        "block/link rank mismatch"))
    all(isfinite, block.h1) && all(isfinite, block.pair_field) &&
        all(isfinite, block.left_feature) && all(isfinite, block.right_feature) &&
        all(isfinite, block.left_core) && all(isfinite, block.right_core) ||
        throw(ArgumentError("block payload must be finite"))
    nothing
end

@inline function _edge_value(block::RHFOuterBlock, physical::Int,
        column::Int, from_left::Bool)
    if from_left
        row = physical - first(block.rows) + 1
        1 <= row <= size(block.left_edge, 1) || return 0.0
        return block.left_edge[row, column]
    end
    row = physical - (last(block.rows) - size(block.right_edge, 1) + 1) + 1
    1 <= row <= size(block.right_edge, 1) || return 0.0
    block.right_edge[row, column]
end

function _empty_outer_block(operator::UnitCellInteraction,
        provenance::UInt; reflected::Bool=false)
    channels = _maximum_channel_rank(operator)
    selection = StateSelection(0, 0, 0, 0.0, :structural_boundary)
    link = RHFStateLink(0, 0, zeros(0, 0), zeros(0, 0), selection)
    block = RHFOuterBlock(1:0, 1:0, 0, 0, 0.0, zeros(0, 0), zeros(0, 0),
        zeros(0, channels), zeros(0, channels), zeros(channels),
        zeros(channels), zeros(0, 0), zeros(0, 0), link, reflected,
        provenance)
    _validate_block(block, operator)
    block
end

function seed_cell_block(one_body::BandedOneBody,
        operator::UnitCellInteraction, cell::Int; provenance::UInt=UInt(0))
    width = UNIT_CELL_WIDTH
    first_site = (cell - 1) * width + 1
    last_site = cell * width
    1 <= first_site <= last_site <= operator.sites || throw(BoundsError(
        1:operator.sites, first_site:last_site))
    h1 = Matrix{Float64}(undef, width, width)
    @inbounds for column = 1:width, row = 1:width
        h1[row, column] = h1_entry(one_body, first_site + row - 1,
            first_site + column - 1)
    end
    pairs = pair_dimension(width)
    channels = _maximum_channel_rank(operator)
    field = zeros(pairs, pairs)
    left = zeros(pairs, channels)
    right = zeros(pairs, channels)
    @inbounds for second = 1:width, first = 1:width
        field[pair_index(first, first), pair_index(second, second)] =
            operator.local_triangle[_local_index(first, second)]
    end
    @inbounds for phase = 1:width, channel = 1:channels
        diagonal = pair_index(phase, phase)
        right[diagonal, channel] = operator.phase_source[phase, channel]
        left[diagonal, channel] = operator.phase_sink[channel, phase]
    end
    selection = StateSelection(0, width, 0, 0.0, :exact)
    link = RHFStateLink(0, width, Matrix{Float64}(I, width, width), selection)
    edge_width = min(width, max(h1_bandwidth(one_body), 1))
    left_edge = zeros(edge_width, width)
    right_edge = zeros(edge_width, width)
    @inbounds for index = 1:edge_width
        left_edge[index, index] = 1.0
        right_edge[index, width - edge_width + index] = 1.0
    end
    block = RHFOuterBlock(cell:cell, first_site:last_site, width, 0, 0.0,
        h1, field, left, right,
        zeros(channels), zeros(channels),
        left_edge, right_edge, link, false, provenance)
    _validate_block(block, operator)
    block
end


function seed_interval_block(one_body::BandedOneBody,
        operator::UnitCellInteraction, cell::Int, rows::UnitRange{Int};
        provenance::UInt=UInt(0), reflected::Bool=false)
    width = length(rows)
    1 <= width <= UNIT_CELL_WIDTH || throw(DimensionMismatch(
        "physical interval width must lie in 1:18"))
    1 <= first(rows) <= last(rows) <= operator.sites || throw(BoundsError(
        1:operator.sites, rows))
    h1 = Matrix{Float64}(undef, width, width)
    @inbounds for column = 1:width, row = 1:width
        h1[row, column] = h1_entry(one_body, first(rows) + row - 1,
            first(rows) + column - 1)
    end
    pairs = pair_dimension(width)
    channels = _maximum_channel_rank(operator)
    field = zeros(pairs, pairs)
    left = zeros(pairs, channels)
    right = zeros(pairs, channels)
    scratch1 = zeros(size(operator.residual_transfer, 1))
    scratch2 = similar(scratch1)
    @inbounds for second_index = 1:width, first_index = 1:width
        i = first(rows) + first_index - 1
        j = first(rows) + second_index - 1
        field[pair_index(first_index, first_index),
            pair_index(second_index, second_index)] =
            _interaction_entry(operator, i, j, scratch1, scratch2)
    end
    @inbounds for local_index = 1:width, channel = 1:channels
        phase = _phase(operator, first(rows) + local_index - 1)
        diagonal = pair_index(local_index, local_index)
        right[diagonal, channel] = operator.phase_source[phase, channel]
        left[diagonal, channel] = operator.phase_sink[channel, phase]
    end
    selection = StateSelection(0, width, 0, 0.0, :exact)
    link = RHFStateLink(0, width, Matrix{Float64}(I, width, width), selection)
    edge_width = min(width, max(h1_bandwidth(one_body), 1))
    left_edge = zeros(edge_width, width)
    right_edge = zeros(edge_width, width)
    @inbounds for index = 1:edge_width
        left_edge[index, index] = 1.0
        right_edge[index, width - edge_width + index] = 1.0
    end
    block = RHFOuterBlock(cell:cell, rows, width, 0, 0.0, h1, field,
        left, right, zeros(channels), zeros(channels), left_edge, right_edge,
        link, reflected, provenance)
    _validate_block(block, operator)
    block
end

function _fill_merged!(workspace::BlockBuildWorkspace, left::RHFOuterBlock,
        right::RHFOuterBlock, one_body::BandedOneBody,
        operator::UnitCellInteraction, required_face::Symbol=:both)
    # Algorithmic donor: PPP 0bcd4f11831a98ae4bf4fb3d09eac5968df9d96c,
    # PeriodicInteractionStateBlocks.jl `_absorb!` and
    # PeriodicInteractionCenters.jl `_transfer_vector!`. PPP's side-specific
    # state owner advances its sole inward channel by one absorbed cell.
    # HFDMRG adapts that ownership to its bidirectional block record: compact
    # lifecycle publication requests one face, while :both remains the dense
    # block oracle used to verify the same arithmetic.
    required_face in (:left, :right, :both) || throw(ArgumentError(
        "required interaction face must be :left, :right, or :both"))
    if !isempty(left.cells) && !isempty(right.cells)
        last(left.cells) + 1 == first(right.cells) || throw(ArgumentError(
            "blocks must be physically adjacent"))
    end
    rank_left, rank_right = left.rank, right.rank
    rank = rank_left + rank_right
    rank <= workspace.maximum_rank || throw(DimensionMismatch(
        "block-build rank capacity is too small"))
    pairs = pair_dimension(rank)
    fill!(@view(workspace.h1[1:rank, 1:rank]), 0.0)
    copyto!(@view(workspace.h1[1:rank_left, 1:rank_left]), left.h1)
    copyto!(@view(workspace.h1[rank_left+1:rank, rank_left+1:rank]), right.h1)
    left_boundary = min(h1_bandwidth(one_body), size(left.right_edge, 1))
    right_boundary = min(h1_bandwidth(one_body), size(right.left_edge, 1))
    @inbounds for a = 1:rank_left, b = 1:rank_right
        coupling = 0.0
        for li = 1:left_boundary, rj = 1:right_boundary
            physical_left = last(left.rows) - left_boundary + li
            physical_right = first(right.rows) + rj - 1
            coupling += left.right_edge[li, a] *
                h1_entry(one_body, physical_left, physical_right) *
                right.left_edge[rj, b]
        end
        workspace.h1[a, rank_left + b] = coupling
        workspace.h1[rank_left + b, a] = coupling
    end
    fill!(@view(workspace.pair_field[1:pairs, 1:pairs]), 0.0)
    fill!(@view(workspace.feature_left[1:pairs, :]), 0.0)
    fill!(@view(workspace.feature_right[1:pairs, :]), 0.0)
    if required_face !== :right
        for j = 1:rank_left, i = 1:j
            destination = pair_index(i, j)
            source = pair_index(i, j)
            @inbounds for channel = 1:size(left.left_feature, 2)
                workspace.feature_left[destination, channel] =
                    left.left_feature[source, channel]
            end
        end
    end
    for j = 1:rank_left, i = 1:j
        destination = pair_index(i, j)
        source = pair_index(i, j)
        for l = 1:rank_left, k = 1:l
            workspace.pair_field[destination, pair_index(k, l)] =
                left.pair_field[source, pair_index(k, l)]
        end
    end
    if required_face !== :left
        for j = 1:rank_right, i = 1:j
            destination = pair_index(rank_left + i, rank_left + j)
            source = pair_index(i, j)
            @inbounds for channel = 1:size(right.right_feature, 2)
                workspace.feature_right[destination, channel] =
                    right.right_feature[source, channel]
            end
        end
    end
    for j = 1:rank_right, i = 1:j
        destination = pair_index(rank_left + i, rank_left + j)
        source = pair_index(i, j)
        for l = 1:rank_right, k = 1:l
            workspace.pair_field[destination,
                pair_index(rank_left + k, rank_left + l)] =
                right.pair_field[source, pair_index(k, l)]
        end
    end
    channels = _maximum_channel_rank(operator)
    fill!(@view(workspace.left_core[1:channels]), 0.0)
    fill!(@view(workspace.right_core[1:channels]), 0.0)
    @inbounds for channel = 1:channels
        required_face !== :right &&
            (workspace.left_core[channel] = left.left_core[channel])
        required_face !== :left &&
            (workspace.right_core[channel] = right.right_core[channel])
    end
    constant = left.constant_energy + right.constant_energy
    @inbounds for channel = 1:channels
        constant = muladd(4.0left.right_core[channel],
            right.left_core[channel], constant)
    end
    # Completed density on one child contributes only Hartree to the active
    # basis of the disjoint child.  Exchange cross terms vanish because a
    # certified completed orbital has no support outside its own block.
    @inbounds for j = 1:rank_right, i = 1:j
        pair = pair_index(i, j)
        value = 0.0
        for channel = 1:channels
            value = muladd(left.right_core[channel],
                right.left_feature[pair, channel], value)
        end
        scale = pair_scale(i, j)
        workspace.h1[rank_left+i, rank_left+j] += 2.0value / scale
        i != j && (workspace.h1[rank_left+j, rank_left+i] += 2.0value / scale)
    end
    @inbounds for j = 1:rank_left, i = 1:j
        pair = pair_index(i, j)
        value = 0.0
        for channel = 1:channels
            value = muladd(right.left_core[channel],
                left.right_feature[pair, channel], value)
        end
        scale = pair_scale(i, j)
        workspace.h1[i, j] += 2.0value / scale
        i != j && (workspace.h1[j, i] += 2.0value / scale)
    end
    if required_face !== :left
        power = length(right.cells)
        transfer_right = _channel_transfer!(workspace.transfer,
            workspace.transfer_scratch, operator, power)
        _record_channel_transfer!(workspace, operator, power)
        @inbounds for pair = 1:pair_dimension(rank_left), outgoing = 1:channels
            value = 0.0
            for incoming = 1:channels
                value = muladd(left.right_feature[pair, incoming],
                    transfer_right[incoming, outgoing], value)
            end
            workspace.feature_right[pair, outgoing] = value
        end
        @inbounds for outgoing = 1:channels
            value = right.right_core[outgoing]
            for incoming = 1:channels
                value = muladd(left.right_core[incoming],
                    transfer_right[incoming, outgoing], value)
            end
            workspace.right_core[outgoing] = value
        end
    end
    if required_face !== :right
        power = length(left.cells)
        transfer_left = _channel_transfer!(workspace.transfer,
            workspace.transfer_scratch, operator, power)
        _record_channel_transfer!(workspace, operator, power)
        # The destination rows are noncontiguous in the enlarged packed
        # ordering, so propagate each exact row without a gathered panel.
        source = 0
        @inbounds for j = 1:rank_right, i = 1:j
            source += 1
            destination = pair_index(rank_left + i, rank_left + j)
            for channel = 1:channels
                value = 0.0
                for old_channel = 1:channels
                    value = muladd(right.left_feature[source, old_channel],
                        transfer_left[channel, old_channel], value)
                end
                workspace.feature_left[destination, channel] = value
            end
        end
        @inbounds for channel = 1:channels
            value = left.left_core[channel]
            for old_channel = 1:channels
                value = muladd(transfer_left[channel, old_channel],
                    right.left_core[old_channel], value)
            end
            workspace.left_core[channel] = value
        end
    end
    for jl = 1:rank_left, il = 1:jl, jr = 1:rank_right, ir = 1:jr
        lp = pair_index(il, jl)
        rp = pair_index(ir, jr)
        destination_r = pair_index(rank_left + ir, rank_left + jr)
        value = 0.0
        for channel = 1:channels
            value = muladd(left.right_feature[lp, channel],
                right.left_feature[rp, channel], value)
        end
        workspace.pair_field[lp, destination_r] = value
        workspace.pair_field[destination_r, lp] = value
    end
    rank, pairs, constant, left.completed_count + right.completed_count
end

function build_outer_block(left::RHFOuterBlock, right::RHFOuterBlock,
        link::RHFStateLink, one_body::BandedOneBody,
        operator::UnitCellInteraction, workspace::BlockBuildWorkspace;
        required_face::Symbol=:both)
    _validate_block(left, operator); _validate_block(right, operator)
    left.provenance == right.provenance || throw(ArgumentError(
        "block provenance is inconsistent"))
    old_rank, old_pairs, constant, completed_before = _fill_merged!(workspace,
        left, right, one_body, operator, required_face)
    old_rank == link.old_rank + link.physical_width || throw(DimensionMismatch(
        "link does not describe the merged block"))
    new_rank = link.active_rank
    newly_completed = link.selection.completed
    if newly_completed > 0
        completed_density = @view workspace.h1_scratch[1:old_rank, 1:old_rank]
        mul!(completed_density, link.completed_map,
            transpose(link.completed_map))
        base_h1 = @view workspace.h1[1:old_rank, 1:old_rank]
        field = @view workspace.pair_field[1:old_pairs, 1:old_pairs]
        one_body_energy = 2.0dot(base_h1, completed_density)
        hartree_energy = 0.0
        exchange_energy = 0.0
        @inbounds for b = 1:old_rank, a = 1:old_rank
            jvalue = 0.0
            kvalue = 0.0
            for d = 1:old_rank, c = 1:old_rank
                jvalue = muladd(completed_density[c, d],
                    _raw_integral(field, a, b, c, d), jvalue)
                kvalue = muladd(completed_density[c, d],
                    _raw_integral(field, a, c, b, d), kvalue)
            end
            hartree_energy = muladd(2.0jvalue,
                completed_density[a, b], hartree_energy)
            exchange_energy = muladd(-kvalue,
                completed_density[a, b], exchange_energy)
            base_h1[a, b] += 2.0jvalue - kvalue
        end
        constant += one_body_energy + hartree_energy + exchange_energy
        channels = _maximum_channel_rank(operator)
        @inbounds for channel = 1:channels
            left_value = workspace.left_core[channel]
            right_value = workspace.right_core[channel]
            for j = 1:old_rank, i = 1:j
                packed = pair_scale(i, j) * completed_density[i, j]
                pair = pair_index(i, j)
                required_face !== :right && (left_value = muladd(packed,
                    workspace.feature_left[pair, channel], left_value))
                required_face !== :left && (right_value = muladd(packed,
                    workspace.feature_right[pair, channel], right_value))
            end
            required_face !== :right &&
                (workspace.left_core[channel] = left_value)
            required_face !== :left &&
                (workspace.right_core[channel] = right_value)
        end
    end
    transform = _pair_transform!(workspace.pair_map, link.close_map,
        old_rank, new_rank)
    new_pairs = pair_dimension(new_rank)
    h1 = Matrix{Float64}(undef, new_rank, new_rank)
    h1_scratch = @view workspace.h1_scratch[1:old_rank, 1:new_rank]
    mul!(h1_scratch, @view(workspace.h1[1:old_rank, 1:old_rank]),
        link.close_map)
    mul!(h1, transpose(link.close_map), h1_scratch)
    scratch = @view workspace.pair_scratch[1:old_pairs, 1:new_pairs]
    mul!(scratch, @view(workspace.pair_field[1:old_pairs, 1:old_pairs]), transform)
    pair_field = Matrix{Float64}(undef, new_pairs, new_pairs)
    mul!(pair_field, transpose(transform), scratch)
    channels = _maximum_channel_rank(operator)
    left_feature = zeros(new_pairs, channels)
    right_feature = zeros(new_pairs, channels)
    required_face !== :right && mul!(left_feature, transpose(transform),
        @view(workspace.feature_left[1:old_pairs, 1:channels]))
    required_face !== :left && mul!(right_feature, transpose(transform),
        @view(workspace.feature_right[1:old_pairs, 1:channels]))
    rows = isempty(left.rows) ? right.rows : isempty(right.rows) ?
        left.rows : first(left.rows):last(right.rows)
    edge_width = min(length(rows), max(h1_bandwidth(one_body), 1))
    left_edge = zeros(edge_width, new_rank)
    right_edge = zeros(edge_width, new_rank)
    children = ((left, 0), (right, left.rank))
    @inbounds for edge = 1:edge_width, column = 1:new_rank
        left_physical = first(rows) + edge - 1
        right_physical = last(rows) - edge_width + edge
        for (child, offset) in children, old = 1:child.rank
            left_edge[edge, column] = muladd(
                _edge_value(child, left_physical, old, true),
                link.close_map[offset + old, column],
                left_edge[edge, column])
            right_edge[edge, column] = muladd(
                _edge_value(child, right_physical, old, false),
                link.close_map[offset + old, column],
                right_edge[edge, column])
        end
    end
    left_core = copy(@view workspace.left_core[1:channels])
    right_core = copy(@view workspace.right_core[1:channels])
    cells = isempty(left.cells) ? right.cells : isempty(right.cells) ?
        left.cells : first(left.cells):last(right.cells)
    block = RHFOuterBlock(cells, rows, new_rank,
        completed_before + newly_completed, constant,
        h1, pair_field, left_feature, right_feature, left_core, right_core,
        left_edge, right_edge, link,
        false, left.provenance)
    _validate_block(block, operator)
    block
end

function Base.reverse(block::RHFOuterBlock)
    RHFOuterBlock(block.cells, block.rows, block.rank, block.completed_count,
        block.constant_energy, copy(block.h1),
        copy(block.pair_field), copy(block.right_feature),
        copy(block.left_feature), copy(block.right_core), copy(block.left_core),
        copy(block.right_edge),
        copy(block.left_edge), block.link, !block.reflected,
        block.provenance)
end

function reflect_block(block::RHFOuterBlock, total_cells::Int,
        total_rows::Int)
    isempty(block.cells) && return RHFOuterBlock(1:0, 1:0, 0, 0, 0.0,
        copy(block.h1), copy(block.pair_field), copy(block.right_feature),
        copy(block.left_feature), copy(block.right_core), copy(block.left_core),
        copy(block.right_edge), copy(block.left_edge), block.link,
        !block.reflected, block.provenance)
    1 <= first(block.cells) <= last(block.cells) <= total_cells ||
        throw(BoundsError(1:total_cells, block.cells))
    reflected_cells = (total_cells - last(block.cells) + 1):(total_cells - first(block.cells) + 1)
    reflected_rows = (total_rows - last(block.rows) + 1):(total_rows -
        first(block.rows) + 1)
    RHFOuterBlock(reflected_cells, reflected_rows, block.rank,
        block.completed_count,
        block.constant_energy, copy(block.h1),
        copy(block.pair_field), copy(block.right_feature),
        copy(block.left_feature), copy(block.right_core), copy(block.left_core),
        copy(block.right_edge),
        copy(block.left_edge), block.link, !block.reflected,
        block.provenance)
end

reflect_block(block::RHFOuterBlock, total_cells::Int) = reflect_block(block,
    total_cells, total_cells * UNIT_CELL_WIDTH)

function block_storage_bytes(block::RHFOuterBlock)
    sizeof(block.h1) + sizeof(block.pair_field) +
        sizeof(block.left_feature) + sizeof(block.right_feature) +
        sizeof(block.left_core) + sizeof(block.right_core) +
        sizeof(block.left_edge) + sizeof(block.right_edge) +
        sizeof(block.link.close_map) + sizeof(block.link.completed_map)
end
