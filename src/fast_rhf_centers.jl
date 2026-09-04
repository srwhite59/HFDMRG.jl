"""
Reusable separated-record RHF center.

Algorithmic donor: PPP 0bcd4f11831a98ae4bf4fb3d09eac5968df9d96c,
`PeriodicInteractionCenters.jl`: `_panels!`, `_core_fields!`,
`_prepare_common!`, `prepare_rhf_periodic_interaction_center!`,
`_rhf_fock_core!`, and `rhf_periodic_interaction_center_fock!`.

Adapted for HFDMRG's `RHFOuterBlock`, two-group moving center, banded-H1
owner, and represented-energy result contract.  This is not byte-identical:
PPP's separated left/center/right panels are evaluated directly from the
equivalent block feature records already owned by HFDMRG.  In particular no
total-center packed-pair matrix is formed.
"""
mutable struct FastRHFPreparedCenter
    maximum_rank::Int
    maximum_rhs::Int
    active_rank::Int
    constant_energy::Float64
    h1::Matrix{Float64}
    density::Matrix{Float64}
    coulomb::Matrix{Float64}
    exchange::Matrix{Float64}
    fock::Matrix{Float64}
    action::Matrix{Float64}
    propagated::Matrix{Float64}
    left_center_exchange::Array{Float64,3}
    right_center_exchange::Array{Float64,3}
    left_right_exchange::Matrix{Float64}
    transfer_zero::Matrix{Float64}
    transfer_one::Matrix{Float64}
    transfer_two::Matrix{Float64}
    transfer::Matrix{Float64}
    transfer_scratch::Matrix{Float64}
    channel_first::Vector{Float64}
    channel_second::Vector{Float64}
    interaction_scratch_first::Vector{Float64}
    interaction_scratch_second::Vector{Float64}
    left::Union{Nothing,RHFOuterBlock}
    right::Union{Nothing,RHFOuterBlock}
    first_rows::UnitRange{Int}
    second_rows::UnitRange{Int}
    operator::UnitCellInteraction
    preparation_calls::Int
    fock_calls::Int
end

function FastRHFPreparedCenter(maximum_rank::Int,
        operator::UnitCellInteraction;
        maximum_rhs::Int=1)
    maximum_rank >= 1 || throw(ArgumentError("maximum rank must be positive"))
    maximum_rhs >= 1 || throw(ArgumentError("RHS capacity must be positive"))
    channels = max(_maximum_channel_rank(operator), 1)
    residual = max(size(operator.residual_transfer, 1), 1)
    maximum_side = max((maximum_rank - 2UNIT_CELL_WIDTH) ÷ 2, 1)
    FastRHFPreparedCenter(maximum_rank, maximum_rhs, 0, 0.0,
        zeros(maximum_rank, maximum_rank),
        zeros(maximum_rank, maximum_rank),
        zeros(maximum_rank, maximum_rank),
        zeros(maximum_rank, maximum_rank),
        zeros(maximum_rank, maximum_rank),
        zeros(maximum_rank, maximum_rhs),
        zeros(pair_dimension(maximum_rank), channels),
        zeros(maximum_side, maximum_side, 2UNIT_CELL_WIDTH),
        zeros(maximum_side, maximum_side, 2UNIT_CELL_WIDTH),
        zeros(maximum_side^2, maximum_side^2),
        Matrix{Float64}(I, channels, channels), zeros(channels, channels),
        zeros(channels, channels), zeros(channels, channels),
        zeros(channels, channels),
        zeros(channels), zeros(channels), zeros(residual), zeros(residual),
        nothing, nothing, 1:0, 1:0, operator, 0, 0)
end

function _prepare_fast_exchange_maps!(prepared::FastRHFPreparedCenter,
        left::RHFOuterBlock, right::RHFOuterBlock, center_rows, operator)
    channels = _maximum_channel_rank(operator)
    first_width = length(prepared.first_rows)
    @inbounds for (block, on_left) in ((left, true), (right, false))
        block.rank == 0 && continue
        feature = on_left ? block.right_feature : block.left_feature
        target = on_left ? prepared.left_center_exchange :
            prepared.right_center_exchange
        for (local_index, site) in enumerate(center_rows)
            gap = on_left ? (local_index <= first_width ? 0 : 1) :
                (local_index <= first_width ? 1 : 0)
            transfer = gap == 0 ? prepared.transfer_zero : prepared.transfer_one
            propagated = @view prepared.propagated[1:pair_dimension(block.rank),
                1:channels]
            mul!(propagated, feature,
                on_left ? transfer : transpose(transfer))
            phase = _phase(operator, site)
            for c = 1:block.rank, a = 1:block.rank
                value = 0.0
                for channel = 1:channels
                    value = muladd(propagated[pair_index(a, c), channel],
                        on_left ? operator.phase_sink[channel, phase] :
                            operator.phase_source[phase, channel], value)
                end
                target[a, c, local_index] = value / pair_scale(a, c)
            end
        end
    end
    if left.rank > 0 && right.rank > 0
        propagated = @view prepared.propagated[1:pair_dimension(left.rank),
            1:channels]
        mul!(propagated, left.right_feature, prepared.transfer_two)
        @inbounds for d = 1:right.rank, b = 1:right.rank,
                c = 1:left.rank, a = 1:left.rank
            value = 0.0
            for channel = 1:channels
                value = muladd(propagated[pair_index(a, c), channel],
                    right.left_feature[pair_index(b, d), channel], value)
            end
            prepared.left_right_exchange[(a-1)*left.rank+c,
                (b-1)*right.rank+d] = value /
                (pair_scale(a, c) * pair_scale(b, d))
        end
    end
    nothing
end

function _fast_h1_boundary!(h1, left::RHFOuterBlock, right::RHFOuterBlock,
        center_rows::UnitRange{Int}, rank_left::Int, offset_right::Int,
        one_body::BandedOneBody)
    left_boundary = min(h1_bandwidth(one_body), size(left.right_edge, 1))
    @inbounds for a = 1:left.rank, b in eachindex(center_rows)
        value = 0.0
        for edge = 1:left_boundary
            site = last(left.rows) - left_boundary + edge
            value = muladd(left.right_edge[edge, a],
                h1_entry(one_body, site, center_rows[b]), value)
        end
        h1[a, rank_left+b] = value
        h1[rank_left+b, a] = value
    end
    right_boundary = min(h1_bandwidth(one_body), size(right.left_edge, 1))
    @inbounds for a in eachindex(center_rows), b = 1:right.rank
        value = 0.0
        for edge = 1:right_boundary
            site = first(right.rows) + edge - 1
            value = muladd(h1_entry(one_body, center_rows[a], site),
                right.left_edge[edge, b], value)
        end
        h1[rank_left+a, offset_right+b] = value
        h1[offset_right+b, rank_left+a] = value
    end
    nothing
end

function _prepare_lifecycle_center!(prepared::FastRHFPreparedCenter,
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
    last(first_rows) + 1 == first(second_rows) || throw(ArgumentError(
        "lifecycle physical intervals must be adjacent"))
    rank_left, rank_right = left.rank, right.rank
    center_rows = first(first_rows):last(second_rows)
    rank = rank_left + length(center_rows) + rank_right
    rank <= prepared.maximum_rank || throw(DimensionMismatch(
        "fast RHF center capacity is too small"))
    h1 = @view prepared.h1[1:rank, 1:rank]
    fill!(h1, 0.0)
    rank_left > 0 && copyto!(@view(h1[1:rank_left, 1:rank_left]), left.h1)
    offset_right = rank_left + length(center_rows)
    rank_right > 0 && copyto!(@view(h1[offset_right+1:rank,
        offset_right+1:rank]), right.h1)
    @inbounds for column in eachindex(center_rows), row in eachindex(center_rows)
        h1[rank_left+row, rank_left+column] = h1_entry(one_body,
            center_rows[row], center_rows[column])
    end
    _fast_h1_boundary!(h1, left, right, center_rows, rank_left, offset_right,
        one_body)

    channels = _maximum_channel_rank(operator)
    copyto!(prepared.transfer_two, _channel_transfer!(prepared.transfer,
        prepared.transfer_scratch, operator, 2))
    transfer2 = prepared.transfer_two
    propagated_core = @view prepared.channel_first[1:channels]
    mul!(propagated_core, transpose(transfer2), left.right_core)
    constant = left.constant_energy + right.constant_energy +
        _core_cross(propagated_core, right.left_core)
    _add_core_hartree!(h1, right, offset_right, propagated_core,
        right.left_feature)
    mul!(@view(prepared.channel_second[1:channels]), transfer2,
        right.left_core)
    _add_core_hartree!(h1, left, 0,
        @view(prepared.channel_second[1:channels]), left.right_feature)

    copyto!(prepared.transfer_one, _channel_transfer!(prepared.transfer,
        prepared.transfer_scratch, operator, 1))
    transfer1 = prepared.transfer_one
    @inbounds for (local_index, site) in enumerate(center_rows)
        phase = _phase(operator, site)
        from_left = 0.0
        to_right = 0.0
        if local_index <= length(first_rows)
            for channel = 1:channels
                from_left = muladd(left.right_core[channel],
                    operator.phase_sink[channel, phase], from_left)
            end
            for incoming = 1:channels, outgoing = 1:channels
                to_right = muladd(operator.phase_source[phase, incoming] *
                    transfer1[incoming, outgoing],
                    right.left_core[outgoing], to_right)
            end
        else
            for incoming = 1:channels, outgoing = 1:channels
                from_left = muladd(left.right_core[incoming] *
                    transfer1[incoming, outgoing],
                    operator.phase_sink[outgoing, phase], from_left)
            end
            for channel = 1:channels
                to_right = muladd(operator.phase_source[phase, channel],
                    right.left_core[channel], to_right)
            end
        end
        h1[rank_left+local_index, rank_left+local_index] +=
            2.0(from_left + to_right)
    end
    prepared.active_rank = rank
    prepared.constant_energy = constant
    prepared.left = left; prepared.right = right
    prepared.first_rows = first_rows; prepared.second_rows = second_rows
    _prepare_fast_exchange_maps!(prepared, left, right, center_rows, operator)
    prepared.preparation_calls += 1
    prepared
end

function _prepare_lifecycle_center!(prepared::FastRHFPreparedCenter,
        left::RHFOuterBlock, right::RHFOuterBlock, cell::Int,
        one_body::BandedOneBody, operator::UnitCellInteraction)
    first_rows = (cell-1)*UNIT_CELL_WIDTH+1:cell*UNIT_CELL_WIDTH
    second_rows = cell*UNIT_CELL_WIDTH+1:(cell+1)*UNIT_CELL_WIDTH
    _prepare_lifecycle_center!(prepared, left, right, cell, first_rows,
        second_rows, one_body, operator)
end

function _fast_block_density_channels!(output, density, block::RHFOuterBlock,
        feature)
    fill!(output, 0.0)
    @inbounds for channel in eachindex(output), j = 1:block.rank, i = 1:j
        output[channel] = muladd(pair_scale(i, j) * density[i, j],
            feature[pair_index(i, j), channel], output[channel])
    end
    output
end

function _fast_block_physical_hartree!(J, block::RHFOuterBlock,
        block_offset::Int, physical_offset::Int, rows::UnitRange{Int},
        density, physical_density, operator, transfer, prepared,
        block_on_left::Bool)
    channels = _maximum_channel_rank(operator)
    block_feature = block_on_left ? block.right_feature : block.left_feature
    physical_channel = @view prepared.channel_first[1:channels]
    fill!(physical_channel, 0.0)
    @inbounds for (local_index, site) in enumerate(rows), channel = 1:channels
        value = physical_density[local_index, local_index]
        physical_channel[channel] = muladd(value,
            block_on_left ? operator.phase_sink[channel, _phase(operator, site)] :
                operator.phase_source[_phase(operator, site), channel],
            physical_channel[channel])
    end
    transformed = @view prepared.channel_second[1:channels]
    if block_on_left
        mul!(transformed, transfer, physical_channel)
    else
        mul!(transformed, transpose(transfer), physical_channel)
    end
    @inbounds for j = 1:block.rank, i = 1:j
        value = 0.0
        for channel = 1:channels
            value = muladd(block_feature[pair_index(i, j), channel],
                transformed[channel], value)
        end
        value /= pair_scale(i, j)
        J[block_offset+i, block_offset+j] += value
        i != j && (J[block_offset+j, block_offset+i] += value)
    end

    block_channel = @view prepared.channel_first[1:channels]
    _fast_block_density_channels!(block_channel, density, block, block_feature)
    if block_on_left
        mul!(transformed, transpose(transfer), block_channel)
    else
        mul!(transformed, transfer, block_channel)
    end
    @inbounds for (local_index, site) in enumerate(rows)
        value = 0.0; phase = _phase(operator, site)
        for channel = 1:channels
            value = muladd(transformed[channel],
                block_on_left ? operator.phase_sink[channel, phase] :
                    operator.phase_source[phase, channel], value)
        end
        J[physical_offset+local_index, physical_offset+local_index] += value
    end
    nothing
end

function _fast_block_block_hartree!(J, left::RHFOuterBlock,
        right::RHFOuterBlock, right_offset::Int, left_density, right_density,
        transfer, prepared)
    channels = size(left.right_feature, 2)
    left_channel = @view prepared.channel_first[1:channels]
    right_channel = @view prepared.channel_second[1:channels]
    _fast_block_density_channels!(left_channel, left_density, left,
        left.right_feature)
    _fast_block_density_channels!(right_channel, right_density, right,
        right.left_feature)
    propagated = @view prepared.propagated[1:pair_dimension(left.rank),
        1:channels]
    mul!(propagated, left.right_feature, transfer)
    @inbounds for j = 1:left.rank, i = 1:j
        value = 0.0
        for channel = 1:channels
            value = muladd(propagated[pair_index(i, j), channel],
                right_channel[channel], value)
        end
        value /= pair_scale(i, j)
        J[i, j] += value; i != j && (J[j, i] += value)
    end
    mul!(right_channel, transpose(transfer), left_channel)
    @inbounds for j = 1:right.rank, i = 1:j
        value = 0.0
        for channel = 1:channels
            value = muladd(right.left_feature[pair_index(i, j), channel],
                right_channel[channel], value)
        end
        value /= pair_scale(i, j)
        J[right_offset+i, right_offset+j] += value
        i != j && (J[right_offset+j, right_offset+i] += value)
    end
    nothing
end

function _fast_block_physical_exchange!(K, density, block::RHFOuterBlock,
        block_offset::Int, physical_offset::Int, center_offset::Int,
        rows::UnitRange{Int}, exchange)
    @inbounds for physical in eachindex(rows), a = 1:block.rank
        value = 0.0
        for c = 1:block.rank
            value = muladd(density[block_offset+c,
                physical_offset+physical], exchange[a, c,
                center_offset+physical], value)
        end
        K[block_offset+a, physical_offset+physical] += value
        K[physical_offset+physical, block_offset+a] += value
    end
    nothing
end

function _fast_block_block_exchange!(K, density, left::RHFOuterBlock,
        right::RHFOuterBlock, right_offset::Int, exchange)
    @inbounds for b = 1:right.rank, a = 1:left.rank
        value = 0.0
        for d = 1:right.rank, c = 1:left.rank
            value = muladd(density[c, right_offset+d],
                exchange[(a-1)*left.rank+c, (b-1)*right.rank+d], value)
        end
        K[a, right_offset+b] += value
        K[right_offset+b, a] += value
    end
    nothing
end

function _fast_rhf_center_energy!(prepared::FastRHFPreparedCenter, covariance,
        occupied::Int; construct_fock::Bool)
    rank = prepared.active_rank
    size(covariance) == (rank, rank) || throw(DimensionMismatch(
        "root and fast RHF center ranks disagree"))
    left = something(prepared.left); right = something(prepared.right)
    operator = something(prepared.operator)
    rows = first(prepared.first_rows):last(prepared.second_rows)
    ml, width, mr = left.rank, length(rows), right.rank
    density = @view prepared.density[1:rank, 1:rank]
    copyto!(density, covariance)
    J = @view prepared.coulomb[1:rank, 1:rank]
    K = @view prepared.exchange[1:rank, 1:rank]
    fill!(J, 0.0); fill!(K, 0.0)
    ml > 0 && _self_spin_fields!(@view(J[1:ml, 1:ml]),
        @view(K[1:ml, 1:ml]), @view(density[1:ml, 1:ml]),
        left.pair_field, ml)
    ro = ml + width
    mr > 0 && _self_spin_fields!(@view(J[ro+1:rank, ro+1:rank]),
        @view(K[ro+1:rank, ro+1:rank]),
        @view(density[ro+1:rank, ro+1:rank]), right.pair_field, mr)
    residual = size(operator.residual_transfer, 1)
    scratch1 = @view prepared.interaction_scratch_first[1:residual]
    scratch2 = @view prepared.interaction_scratch_second[1:residual]
    @inbounds for b = 1:width, a = 1:width
        value = _interaction_entry(operator, rows[a], rows[b], scratch1,
            scratch2)
        K[ml+a, ml+b] = value * density[ml+a, ml+b]
        J[ml+a, ml+a] = muladd(value, density[ml+b, ml+b],
            J[ml+a, ml+a])
    end
    identity_transfer = @view prepared.transfer_zero[1:_maximum_channel_rank(operator),
        1:_maximum_channel_rank(operator)]
    transfer1 = @view prepared.transfer_one[1:_maximum_channel_rank(operator),
        1:_maximum_channel_rank(operator)]
    transfer2 = @view prepared.transfer_two[1:_maximum_channel_rank(operator),
        1:_maximum_channel_rank(operator)]
    first_width = length(prepared.first_rows)
    ml > 0 && _fast_block_physical_hartree!(J, left, 0, ml,
        prepared.first_rows, @view(density[1:ml, 1:ml]),
        @view(density[ml+1:ml+first_width, ml+1:ml+first_width]), operator,
        identity_transfer, prepared, true)
    ml > 0 && _fast_block_physical_exchange!(K, density, left, 0, ml,
        0, prepared.first_rows, prepared.left_center_exchange)
    ml > 0 && _fast_block_physical_hartree!(J, left, 0, ml+first_width,
        prepared.second_rows, @view(density[1:ml, 1:ml]),
        @view(density[ml+first_width+1:ml+width,
        ml+first_width+1:ml+width]), operator, transfer1, prepared, true)
    ml > 0 && _fast_block_physical_exchange!(K, density, left, 0,
        ml+first_width, first_width, prepared.second_rows,
        prepared.left_center_exchange)
    mr > 0 && _fast_block_physical_hartree!(J, right, ro, ml+first_width,
        prepared.second_rows, @view(density[ro+1:rank, ro+1:rank]),
        @view(density[ml+first_width+1:ml+width,
            ml+first_width+1:ml+width]), operator, identity_transfer,
        prepared, false)
    mr > 0 && _fast_block_physical_exchange!(K, density, right, ro,
        ml+first_width, first_width, prepared.second_rows,
        prepared.right_center_exchange)
    mr > 0 && _fast_block_physical_hartree!(J, right, ro, ml,
        prepared.first_rows, @view(density[ro+1:rank, ro+1:rank]),
        @view(density[ml+1:ml+first_width, ml+1:ml+first_width]), operator,
        transfer1, prepared, false)
    mr > 0 && _fast_block_physical_exchange!(K, density, right, ro, ml,
        0, prepared.first_rows, prepared.right_center_exchange)
    ml > 0 && mr > 0 && _fast_block_block_hartree!(J, left, right, ro,
        @view(density[1:ml, 1:ml]), @view(density[ro+1:rank, ro+1:rank]),
        transfer2, prepared)
    ml > 0 && mr > 0 && _fast_block_block_exchange!(K, density, left, right,
        ro, prepared.left_right_exchange)
    h1 = @view prepared.h1[1:rank, 1:rank]
    if construct_fock
        fock = @view prepared.fock[1:rank, 1:rank]
        @. fock = h1 + 2.0J - K
        prepared.fock_calls += 1
    end
    one = 2.0dot(h1, density)
    hartree = 2.0dot(J, density)
    exchange = -dot(K, density)
    (total=prepared.constant_energy + one + hartree + exchange,
        core=prepared.constant_energy, one_body=one, hartree, exchange)
end

function _center_energy_fock!(prepared::FastRHFPreparedCenter, covariance,
        occupied::Int)
    _fast_rhf_center_energy!(prepared, covariance, occupied;
        construct_fock=true)
end

function _center_true_energy!(prepared::FastRHFPreparedCenter, covariance,
        occupied::Int)
    _fast_rhf_center_energy!(prepared, covariance, occupied;
        construct_fock=false)
end

function apply_center_fock!(output::AbstractVecOrMat{Float64},
        prepared::FastRHFPreparedCenter, input::AbstractVecOrMat{Float64})
    rank = prepared.active_rank
    size(input, 1) == rank && size(output) == size(input) ||
        throw(DimensionMismatch("fast RHF center action dimensions disagree"))
    mul!(output, @view(prepared.fock[1:rank, 1:rank]), input)
end

function fast_rhf_center_storage_bytes(prepared::FastRHFPreparedCenter)
    sum(sizeof, (prepared.h1, prepared.density, prepared.coulomb,
        prepared.exchange, prepared.fock, prepared.action,
        prepared.propagated, prepared.left_center_exchange,
        prepared.right_center_exchange, prepared.left_right_exchange,
        prepared.transfer_zero, prepared.transfer_one,
        prepared.transfer_two, prepared.transfer, prepared.transfer_scratch,
        prepared.channel_first, prepared.channel_second,
        prepared.interaction_scratch_first,
        prepared.interaction_scratch_second))
end
