mutable struct PreparedRHFCenter
    maximum_rank::Int
    maximum_pairs::Int
    maximum_rhs::Int
    active_rank::Int
    constant_energy::Float64
    h1::Matrix{Float64}
    pair_field::Matrix{Float64}
    density::Matrix{Float64}
    coulomb::Matrix{Float64}
    exchange::Matrix{Float64}
    fock::Matrix{Float64}
    pair_density::Vector{Float64}
    action::Matrix{Float64}
    feature_scratch::Matrix{Float64}
    transfer::Matrix{Float64}
    transfer_scratch::Matrix{Float64}
    core_scratch::Vector{Float64}
end

function PreparedRHFCenter(maximum_rank::Int, channel_rank::Int;
        maximum_rhs::Int=1)
    maximum_rank >= 1 || throw(ArgumentError("maximum rank must be positive"))
    maximum_rhs >= 1 || throw(ArgumentError("RHS capacity must be positive"))
    pairs = pair_dimension(maximum_rank)
    PreparedRHFCenter(maximum_rank, pairs, maximum_rhs, 0, 0.0,
        fill(NaN, maximum_rank, maximum_rank),
        fill(NaN, pairs, pairs), fill(NaN, maximum_rank, maximum_rank),
        fill(NaN, maximum_rank, maximum_rank),
        fill(NaN, maximum_rank, maximum_rank),
        fill(NaN, maximum_rank, maximum_rank), fill(NaN, pairs),
        fill(NaN, maximum_rank, maximum_rhs),
        fill(NaN, pairs, channel_rank),
        fill(NaN, channel_rank, channel_rank),
        fill(NaN, channel_rank, channel_rank), fill(NaN, channel_rank))
end

function _add_core_hartree!(h1, block::RHFOuterBlock, offset::Int,
        core, feature)
    channels = length(core)
    @inbounds for j = 1:block.rank, i = 1:j
        value = 0.0
        pair = pair_index(i, j)
        for channel = 1:channels
            value = muladd(core[channel], feature[pair, channel], value)
        end
        value *= 2.0 / pair_scale(i, j)
        h1[offset+i, offset+j] += value
        i != j && (h1[offset+j, offset+i] += value)
    end
    nothing
end

@inline function _core_cross(first, second)
    value = 0.0
    @inbounds for channel in eachindex(first, second)
        value = muladd(first[channel], second[channel], value)
    end
    4.0value
end

function _copy_pair_block!(destination, block::RHFOuterBlock, offset::Int)
    @inbounds for j = 1:block.rank, i = 1:j
        source_first = pair_index(i, j)
        target_first = pair_index(offset + i, offset + j)
        for l = 1:block.rank, k = 1:l
            destination[target_first, pair_index(offset + k, offset + l)] =
                block.pair_field[source_first, pair_index(k, l)]
        end
    end
end

function _fill_cross_pair!(destination, left::RHFOuterBlock, left_offset::Int,
        right::RHFOuterBlock, right_offset::Int, propagated_right,
        channels::Int)
    @inbounds for jl = 1:left.rank, il = 1:jl,
            jr = 1:right.rank, ir = 1:jr
        lp = pair_index(il, jl)
        rp = pair_index(ir, jr)
        value = 0.0
        for channel = 1:channels
            value = muladd(propagated_right[lp, channel],
                right.left_feature[rp, channel], value)
        end
        first_pair = pair_index(left_offset + il, left_offset + jl)
        second_pair = pair_index(right_offset + ir, right_offset + jr)
        destination[first_pair, second_pair] = value
        destination[second_pair, first_pair] = value
    end
end

function _fill_h1_cross!(h1, left::RHFOuterBlock, left_offset::Int,
        right::RHFOuterBlock, right_offset::Int, one_body::BandedOneBody)
    left_boundary = min(h1_bandwidth(one_body), size(left.right_edge, 1))
    right_boundary = min(h1_bandwidth(one_body), size(right.left_edge, 1))
    @inbounds for a = 1:left.rank, b = 1:right.rank
        value = 0.0
        for li = 1:left_boundary, rj = 1:right_boundary
            site_left = last(left.rows) - left_boundary + li
            site_right = first(right.rows) + rj - 1
            value = muladd(left.right_edge[li, a] *
                h1_entry(one_body, site_left, site_right),
                right.left_edge[rj, b], value)
        end
        h1[left_offset + a, right_offset + b] = value
        h1[right_offset + b, left_offset + a] = value
    end
end

function prepare_two_sided_center!(prepared::PreparedRHFCenter,
        left::RHFOuterBlock, center::RHFOuterBlock, right::RHFOuterBlock,
        one_body::BandedOneBody, operator::UnitCellInteraction)
    last(left.cells) + 1 == first(center.cells) &&
        last(center.cells) + 1 == first(right.cells) || throw(ArgumentError(
            "two-sided center blocks must be adjacent"))
    _validate_block(left, operator); _validate_block(center, operator)
    _validate_block(right, operator)
    left.provenance == center.provenance &&
        center.provenance == right.provenance ||
        throw(ArgumentError("center provenance is inconsistent"))
    rank_left, rank_center, rank_right = left.rank, center.rank, right.rank
    rank = rank_left + rank_center + rank_right
    rank <= prepared.maximum_rank || throw(DimensionMismatch(
        "prepared-center rank capacity is too small"))
    pairs = pair_dimension(rank)
    channels = _maximum_channel_rank(operator)
    size(prepared.transfer, 1) >= channels || throw(DimensionMismatch(
        "prepared-center channel capacity is too small"))
    h1 = @view prepared.h1[1:rank, 1:rank]
    field = @view prepared.pair_field[1:pairs, 1:pairs]
    fill!(h1, 0.0)
    fill!(field, 0.0)
    offsets = (0, rank_left, rank_left + rank_center)
    blocks = (left, center, right)
    for (block, offset) in zip(blocks, offsets)
        copyto!(@view(h1[offset+1:offset+block.rank,
            offset+1:offset+block.rank]), block.h1)
        _copy_pair_block!(field, block, offset)
    end
    _fill_h1_cross!(h1, left, offsets[1], center, offsets[2], one_body)
    _fill_h1_cross!(h1, center, offsets[2], right, offsets[3], one_body)
    _fill_cross_pair!(field, left, offsets[1], center, offsets[2],
        left.right_feature, channels)
    _fill_cross_pair!(field, center, offsets[2], right, offsets[3],
        center.right_feature, channels)
    constant = left.constant_energy + center.constant_energy +
        right.constant_energy
    constant += _core_cross(left.right_core, center.left_core)
    constant += _core_cross(center.right_core, right.left_core)
    _add_core_hartree!(h1, center, offsets[2], left.right_core,
        center.left_feature)
    _add_core_hartree!(h1, left, offsets[1], center.left_core,
        left.right_feature)
    _add_core_hartree!(h1, right, offsets[3], center.right_core,
        right.left_feature)
    _add_core_hartree!(h1, center, offsets[2], right.left_core,
        center.right_feature)
    transfer = _channel_transfer!(prepared.transfer,
        prepared.transfer_scratch, operator,
        length(center.cells))
    # Reuse action storage as the propagated left feature; its row capacity is
    # rank rather than pair rank, so perform the modest channel contraction
    # directly when filling the nonadjacent cross term.
    @inbounds for outgoing = 1:channels
        prepared.core_scratch[outgoing] = 0.0
        for incoming = 1:channels
            prepared.core_scratch[outgoing] = muladd(
                left.right_core[incoming], transfer[incoming, outgoing],
                prepared.core_scratch[outgoing])
        end
    end
    constant += _core_cross(@view(prepared.core_scratch[1:channels]),
        right.left_core)
    _add_core_hartree!(h1, right, offsets[3],
        @view(prepared.core_scratch[1:channels]), right.left_feature)
    @inbounds for incoming = 1:channels
        prepared.core_scratch[incoming] = 0.0
        for outgoing = 1:channels
            prepared.core_scratch[incoming] = muladd(
                transfer[incoming, outgoing], right.left_core[outgoing],
                prepared.core_scratch[incoming])
        end
    end
    _add_core_hartree!(h1, left, offsets[1],
        @view(prepared.core_scratch[1:channels]), left.right_feature)
    @inbounds for jl = 1:left.rank, il = 1:jl,
            jr = 1:right.rank, ir = 1:jr
        lp = pair_index(il, jl)
        rp = pair_index(ir, jr)
        value = 0.0
        for outgoing = 1:channels, incoming = 1:channels
            value = muladd(left.right_feature[lp, incoming] *
                transfer[incoming, outgoing],
                right.left_feature[rp, outgoing], value)
        end
        first_pair = pair_index(il, jl)
        second_pair = pair_index(offsets[3] + ir, offsets[3] + jr)
        field[first_pair, second_pair] = value
        field[second_pair, first_pair] = value
    end
    prepared.active_rank = rank
    prepared.constant_energy = constant
    prepared
end

function prepare_terminal_center!(prepared::PreparedRHFCenter,
        left::RHFOuterBlock, right::RHFOuterBlock,
        one_body::BandedOneBody, operator::UnitCellInteraction)
    last(left.cells) + 1 == first(right.cells) || throw(ArgumentError(
        "terminal-center blocks must be adjacent"))
    _validate_block(left, operator); _validate_block(right, operator)
    left.provenance == right.provenance || throw(ArgumentError(
        "terminal-center provenance is inconsistent"))
    rank = left.rank + right.rank
    rank <= prepared.maximum_rank || throw(DimensionMismatch(
        "prepared-center rank capacity is too small"))
    pairs = pair_dimension(rank)
    channels = _maximum_channel_rank(operator)
    h1 = @view prepared.h1[1:rank, 1:rank]
    field = @view prepared.pair_field[1:pairs, 1:pairs]
    fill!(h1, 0.0); fill!(field, 0.0)
    copyto!(@view(h1[1:left.rank, 1:left.rank]), left.h1)
    copyto!(@view(h1[left.rank+1:rank, left.rank+1:rank]), right.h1)
    _copy_pair_block!(field, left, 0)
    _copy_pair_block!(field, right, left.rank)
    _fill_h1_cross!(h1, left, 0, right, left.rank, one_body)
    _fill_cross_pair!(field, left, 0, right, left.rank,
        left.right_feature, channels)
    constant = left.constant_energy + right.constant_energy +
        _core_cross(left.right_core, right.left_core)
    _add_core_hartree!(h1, right, left.rank, left.right_core,
        right.left_feature)
    _add_core_hartree!(h1, left, 0, right.left_core,
        left.right_feature)
    prepared.active_rank = rank
    prepared.constant_energy = constant
    prepared
end

function _physical_feature(operator::UnitCellInteraction, phase::Int,
        channel::Int, left::Bool)
    left ? operator.phase_sink[channel, phase] :
        operator.phase_source[phase, channel]
end

function _fill_block_physical_cross!(field, block::RHFOuterBlock,
        block_offset::Int, physical_offset::Int, physical_rows::UnitRange{Int},
        operator::UnitCellInteraction, propagated, block_on_left::Bool)
    channels = _maximum_channel_rank(operator)
    @inbounds for j = 1:block.rank, i = 1:j,
            local_index in eachindex(physical_rows)
        phase = _phase(operator, physical_rows[local_index])
        block_pair = pair_index(i, j)
        physical_pair = pair_index(physical_offset + local_index,
            physical_offset + local_index)
        value = 0.0
        for channel = 1:channels
            value = muladd(propagated[block_pair, channel],
                _physical_feature(operator, phase, channel, block_on_left),
                value)
        end
        target = pair_index(block_offset + i, block_offset + j)
        field[target, physical_pair] = value
        field[physical_pair, target] = value
    end
    nothing
end

function _propagate_feature!(output, input, transfer, pairs, channels)
    @inbounds for pair = 1:pairs, outgoing = 1:channels
        value = 0.0
        for incoming = 1:channels
            value = muladd(input[pair, incoming],
                transfer[incoming, outgoing], value)
        end
        output[pair, outgoing] = value
    end
    output
end

"""Prepare `[left | first rows | second rows | right]` in reusable storage."""
function _prepare_lifecycle_center!(prepared::PreparedRHFCenter,
        left::RHFOuterBlock, right::RHFOuterBlock, cell::Int,
        first_rows::UnitRange{Int}, second_rows::UnitRange{Int},
        one_body::BandedOneBody, operator::UnitCellInteraction)
    _validate_block(left, operator); _validate_block(right, operator)
    left.provenance == right.provenance || throw(ArgumentError(
        "lifecycle-center provenance is inconsistent"))
    isempty(left.cells) || last(left.cells) + 1 == cell || throw(ArgumentError(
        "left lifecycle interval is inconsistent"))
    isempty(right.cells) || cell + 2 == first(right.cells) ||
        throw(ArgumentError("right lifecycle interval is inconsistent"))
    !isempty(first_rows) && !isempty(second_rows) &&
        last(first_rows) + 1 == first(second_rows) &&
        1 <= first(first_rows) <= last(second_rows) <= operator.sites ||
        throw(BoundsError(1:operator.sites,
            first(first_rows):last(second_rows)))
    width_first, width_second = length(first_rows), length(second_rows)
    1 <= width_first <= UNIT_CELL_WIDTH && 1 <= width_second <= UNIT_CELL_WIDTH ||
        throw(DimensionMismatch("lifecycle physical groups exceed unit-cell capacity"))
    rank_left, rank_right = left.rank, right.rank
    rank = rank_left + width_first + width_second + rank_right
    rank <= prepared.maximum_rank || throw(DimensionMismatch(
        "lifecycle-center rank capacity is too small"))
    pairs = pair_dimension(rank)
    channels = _maximum_channel_rank(operator)
    h1 = @view prepared.h1[1:rank, 1:rank]
    field = @view prepared.pair_field[1:pairs, 1:pairs]
    fill!(h1, 0.0); fill!(field, 0.0)
    offset_first = rank_left
    offset_second = rank_left + width_first
    offset_right = offset_second + width_second
    rank_left > 0 && copyto!(@view(h1[1:rank_left, 1:rank_left]), left.h1)
    rank_right > 0 && copyto!(@view(h1[offset_right+1:rank,
        offset_right+1:rank]), right.h1)
    _copy_pair_block!(field, left, 0)
    _copy_pair_block!(field, right, offset_right)
    center_rows = first(first_rows):last(second_rows)
    @inbounds for local_column in eachindex(center_rows),
            local_row in eachindex(center_rows)
        h1[rank_left + local_row, rank_left + local_column] = h1_entry(
            one_body, center_rows[local_row], center_rows[local_column])
    end
    scratch1 = @view prepared.core_scratch[1:size(operator.residual_transfer, 1)]
    scratch2 = @view prepared.transfer[1:size(operator.residual_transfer, 1), 1]
    for (rows, offset) in ((first_rows, offset_first),
            (second_rows, offset_second))
        @inbounds for second in eachindex(rows), first in eachindex(rows)
            field[pair_index(offset + first, offset + first),
                pair_index(offset + second, offset + second)] =
                _interaction_entry(operator, rows[first], rows[second],
                    scratch1, scratch2)
        end
    end
    # Boundary H1 contractions are written directly to avoid constructing
    # physical block records on the production path.
    boundary = min(h1_bandwidth(one_body), size(left.right_edge, 1))
    @inbounds for a = 1:rank_left, b = 1:length(center_rows)
        value = 0.0
        for li = 1:boundary
            site_left = last(left.rows) - boundary + li
            value = muladd(left.right_edge[li, a],
                h1_entry(one_body, site_left, center_rows[b]), value)
        end
        h1[a, rank_left+b] = value; h1[rank_left+b, a] = value
    end
    boundary = min(h1_bandwidth(one_body), size(right.left_edge, 1))
    @inbounds for a = 1:length(center_rows), b = 1:rank_right
        value = 0.0
        for rj = 1:boundary
            site_right = first(right.rows) + rj - 1
            value = muladd(h1_entry(one_body,
                center_rows[a], site_right),
                right.left_edge[rj, b], value)
        end
        h1[rank_left+a, offset_right+b] = value
        h1[offset_right+b, rank_left+a] = value
    end
    _fill_block_physical_cross!(field, left, 0, offset_first, first_rows,
        operator, left.right_feature, true)
    @inbounds for first in eachindex(first_rows), second in eachindex(second_rows)
        value = 0.0
        first_phase = _phase(operator, first_rows[first])
        second_phase = _phase(operator, second_rows[second])
        for channel = 1:channels
            value = muladd(operator.phase_source[first_phase, channel],
                operator.phase_sink[channel, second_phase], value)
        end
        p = pair_index(offset_first+first, offset_first+first)
        q = pair_index(offset_second+second, offset_second+second)
        field[p, q] = value; field[q, p] = value
    end
    _fill_block_physical_cross!(field, right, offset_right, offset_second,
        second_rows, operator, right.left_feature, false)
    transfer = _channel_transfer!(prepared.transfer,
        prepared.transfer_scratch, operator, 1)
    propagated = @view prepared.feature_scratch[
        1:pair_dimension(rank_left), 1:channels]
    _propagate_feature!(propagated, left.right_feature, transfer,
        pair_dimension(rank_left), channels)
    _fill_block_physical_cross!(field, left, 0, offset_second, second_rows,
        operator, propagated, true)
    # Cell one to the right block crosses cell two's group transfer.
    @inbounds for local_index in eachindex(first_rows), jr = 1:right.rank,
            ir = 1:jr
        phase = _phase(operator, first_rows[local_index])
        value = 0.0
        rp = pair_index(ir, jr)
        for incoming = 1:channels, outgoing = 1:channels
            value = muladd(operator.phase_source[phase, incoming] *
                transfer[incoming, outgoing],
                right.left_feature[rp, outgoing], value)
        end
        p = pair_index(offset_first+local_index, offset_first+local_index)
        q = pair_index(offset_right+ir, offset_right+jr)
        field[p, q] = value; field[q, p] = value
    end
    transfer2 = _channel_transfer!(prepared.transfer,
        prepared.transfer_scratch, operator, 2)
    propagated = @view prepared.feature_scratch[
        1:pair_dimension(rank_left), 1:channels]
    _propagate_feature!(propagated, left.right_feature, transfer2,
        pair_dimension(rank_left), channels)
    _fill_cross_pair!(field, left, 0, right, offset_right, propagated, channels)
    constant = left.constant_energy + right.constant_energy
    @inbounds for outgoing = 1:channels
        prepared.core_scratch[outgoing] = 0.0
        for incoming = 1:channels
            prepared.core_scratch[outgoing] = muladd(left.right_core[incoming],
                transfer2[incoming, outgoing], prepared.core_scratch[outgoing])
        end
    end
    constant += _core_cross(@view(prepared.core_scratch[1:channels]),
        right.left_core)
    _add_core_hartree!(h1, right, offset_right,
        @view(prepared.core_scratch[1:channels]), right.left_feature)
    transfer = _channel_transfer!(prepared.transfer,
        prepared.transfer_scratch, operator, 1)
    @inbounds for local_index in eachindex(first_rows)
        phase = _phase(operator, first_rows[local_index])
        first_value = 0.0; second_value = 0.0
        for channel = 1:channels
            first_value = muladd(left.right_core[channel],
                operator.phase_sink[channel, phase], first_value)
            for outgoing = 1:channels
                second_value = muladd(left.right_core[channel] *
                    transfer[channel, outgoing],
                    operator.phase_sink[outgoing, phase], second_value)
            end
        end
        h1[offset_first+local_index, offset_first+local_index] += 2.0first_value
    end
    @inbounds for local_index in eachindex(second_rows)
        phase = _phase(operator, second_rows[local_index])
        second_value = 0.0
        for channel = 1:channels, outgoing = 1:channels
            second_value = muladd(left.right_core[channel] *
                transfer[channel, outgoing],
                operator.phase_sink[outgoing, phase], second_value)
        end
        h1[offset_second+local_index, offset_second+local_index] += 2.0second_value
    end
    transfer2 = _channel_transfer!(prepared.transfer,
        prepared.transfer_scratch, operator, 2)
    @inbounds for incoming = 1:channels
        prepared.core_scratch[incoming] = 0.0
        for outgoing = 1:channels
            prepared.core_scratch[incoming] = muladd(
                transfer2[incoming, outgoing], right.left_core[outgoing],
                prepared.core_scratch[incoming])
        end
    end
    _add_core_hartree!(h1, left, 0,
        @view(prepared.core_scratch[1:channels]), left.right_feature)
    transfer = _channel_transfer!(prepared.transfer,
        prepared.transfer_scratch, operator, 1)
    @inbounds for local_index in eachindex(second_rows)
        phase = _phase(operator, second_rows[local_index])
        second_value = 0.0; first_value = 0.0
        for channel = 1:channels
            second_value = muladd(right.left_core[channel],
                operator.phase_source[phase, channel], second_value)
            for outgoing = 1:channels
                first_value = muladd(operator.phase_source[phase, channel] *
                    transfer[channel, outgoing], right.left_core[outgoing],
                    first_value)
            end
        end
        h1[offset_second+local_index, offset_second+local_index] += 2.0second_value
    end
    @inbounds for local_index in eachindex(first_rows)
        phase = _phase(operator, first_rows[local_index])
        first_value = 0.0
        for channel = 1:channels, outgoing = 1:channels
            first_value = muladd(operator.phase_source[phase, channel] *
                transfer[channel, outgoing], right.left_core[outgoing],
                first_value)
        end
        h1[offset_first+local_index, offset_first+local_index] += 2.0first_value
    end
    prepared.active_rank = rank
    prepared.constant_energy = constant
    prepared
end

function _prepare_lifecycle_center!(prepared::PreparedRHFCenter,
        left::RHFOuterBlock, right::RHFOuterBlock, cell::Int,
        one_body::BandedOneBody, operator::UnitCellInteraction)
    first_rows = (cell - 1) * UNIT_CELL_WIDTH + 1:cell * UNIT_CELL_WIDTH
    second_rows = cell * UNIT_CELL_WIDTH + 1:(cell + 1) * UNIT_CELL_WIDTH
    _prepare_lifecycle_center!(prepared, left, right, cell, first_rows,
        second_rows, one_body, operator)
end

function _check_center_density(prepared::PreparedRHFCenter, covariance,
        occupied::Int)
    rank = prepared.active_rank
    size(covariance) == (rank, rank) || throw(DimensionMismatch(
        "root and prepared-center ranks disagree"))
    0 <= occupied <= rank || throw(DimensionMismatch(
        "root occupation exceeds the prepared center"))
    nothing
end

function center_energy_fock!(prepared::PreparedRHFCenter, root::RHFRoot)
    _center_energy_fock!(prepared, root.covariance, root.occupied)
end

function _center_energy_fock!(prepared::PreparedRHFCenter, covariance,
        occupied::Int)
    _check_center_density(prepared, covariance, occupied)
    rank = prepared.active_rank
    pairs = pair_dimension(rank)
    density = @view prepared.density[1:rank, 1:rank]
    copyto!(density, covariance)
    pair_density = @view prepared.pair_density[1:pairs]
    @inbounds for j = 1:rank, i = 1:j
        pair_density[pair_index(i, j)] = pair_scale(i, j) * density[i, j]
    end
    h1 = @view prepared.h1[1:rank, 1:rank]
    field = @view prepared.pair_field[1:pairs, 1:pairs]
    coulomb = @view prepared.coulomb[1:rank, 1:rank]
    exchange = @view prepared.exchange[1:rank, 1:rank]
    fock = @view prepared.fock[1:rank, 1:rank]
    fill!(coulomb, 0.0)
    fill!(exchange, 0.0)
    @inbounds for b = 1:rank, a = 1:rank
        jvalue = 0.0
        kvalue = 0.0
        for d = 1:rank, c = 1:rank
            jvalue = muladd(density[c, d],
                _raw_integral(field, a, b, c, d), jvalue)
            kvalue = muladd(density[c, d],
                _raw_integral(field, a, c, b, d), kvalue)
        end
        coulomb[a, b] = jvalue
        exchange[a, b] = kvalue
        fock[a, b] = h1[a, b] + 2.0jvalue - kvalue
    end
    one_body = 2.0dot(h1, density)
    hartree = 2.0dot(coulomb, density)
    exchange_energy = -dot(exchange, density)
    (total=prepared.constant_energy + one_body + hartree + exchange_energy,
        core=prepared.constant_energy, one_body=one_body,
        hartree=hartree, exchange=exchange_energy)
end

function apply_center_fock!(output::AbstractVecOrMat{Float64},
        prepared::PreparedRHFCenter, input::AbstractVecOrMat{Float64})
    rank = prepared.active_rank
    size(input, 1) == rank && size(output) == size(input) ||
        throw(DimensionMismatch("center action dimensions disagree"))
    Base.mightalias(output, input) && throw(ArgumentError(
        "center action output must not alias input"))
    mul!(output, @view(prepared.fock[1:rank, 1:rank]), input)
end

function center_storage_formula(rank::Int, channels::Int, rhs::Int)
    pairs = pair_dimension(rank)
    (h1_density_fock_matrices=5rank^2,
        pair_field=pairs^2, pair_vector=pairs,
        channel_transfer=2channels^2, channel_scratch=channels,
        feature_scratch=pairs * channels,
        rhs_panel=rank * rhs)
end

function block_storage_formula(rank::Int, channels::Int, edge_width::Int,
        link_rows::Int; completed::Int=0)
    pairs = pair_dimension(rank)
    (h1=rank^2, pair_field=pairs^2,
        boundary_features=2pairs * channels,
        completed_core_channels=2channels,
        physical_edges=2edge_width * rank,
        reversible_link=link_rows * (rank + completed))
end
