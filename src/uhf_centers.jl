mutable struct PreparedUHFCenter
    alpha::PreparedRHFCenter
    beta::PreparedRHFCenter
    cross_alpha::Matrix{Float64}
    cross_beta::Matrix{Float64}
    cross_weight::Vector{Float64}
    alpha_factor::Vector{Float64}
    beta_factor::Vector{Float64}
    transfer::Matrix{Float64}
    transfer_scratch::Matrix{Float64}
    cross_rank::Int
    constant_energy::Float64
end

function PreparedUHFCenter(maximum_alpha::Int, maximum_beta::Int,
        channel_rank::Int; maximum_rhs::Int=3)
    min(maximum_alpha, maximum_beta) >= 1 || throw(ArgumentError(
        "UHF center capacities must be positive"))
    Ka, Kb = pair_dimension(maximum_alpha), pair_dimension(maximum_beta)
    columns = 2min(Ka, Kb) + 2UNIT_CELL_WIDTH^2 + 12channel_rank
    PreparedUHFCenter(PreparedRHFCenter(maximum_alpha, channel_rank;
            maximum_rhs), PreparedRHFCenter(maximum_beta, channel_rank;
            maximum_rhs), zeros(Ka, columns), zeros(Kb, columns),
        zeros(columns), zeros(columns), zeros(columns),
        zeros(channel_rank, channel_rank), zeros(channel_rank, channel_rank),
        0, 0.0)
end

@inline function _component_rank(block::UHFOuterBlock, spin::Symbol)
    spin === :alpha ? block.alpha.rank : block.beta.rank
end

@inline function _component_feature(block::UHFOuterBlock, spin::Symbol,
        local_pair::Int, channel::Int, left::Bool,
        operator::UnitCellInteraction)
    record = spin === :alpha ? block.alpha : block.beta
    left ? record.left_feature[local_pair, channel] :
        record.right_feature[local_pair, channel]
end

@inline function _physical_feature_value(rows::UnitRange{Int}, local_pair::Int,
        channel::Int, left::Bool, operator::UnitCellInteraction)
    local_index = 0
    @inbounds for index in eachindex(rows)
        pair_index(index, index) == local_pair && (local_index = index; break)
    end
    local_index == 0 && return 0.0
    phase = _phase(operator, rows[local_index])
    left ? operator.phase_sink[channel, phase] :
        operator.phase_source[phase, channel]
end

@inline function _center_component_rank(component::Int, spin::Symbol,
        left::UHFOuterBlock, first_rows, second_rows,
        right::UHFOuterBlock)
    component == 1 && return _component_rank(left, spin)
    component == 2 && return length(first_rows)
    component == 3 && return length(second_rows)
    _component_rank(right, spin)
end

@inline function _center_component_feature(component::Int, spin::Symbol,
        local_pair::Int, channel::Int, use_left::Bool,
        left::UHFOuterBlock, first_rows, second_rows,
        right::UHFOuterBlock, operator::UnitCellInteraction)
    component == 1 && return _component_feature(left, spin, local_pair,
        channel, use_left, operator)
    component == 2 && return _physical_feature_value(first_rows, local_pair,
        channel, use_left, operator)
    component == 3 && return _physical_feature_value(second_rows, local_pair,
        channel, use_left, operator)
    _component_feature(right, spin, local_pair, channel, use_left, operator)
end

function _append_internal_cross!(prepared::PreparedUHFCenter, column::Int,
        block::UHFOuterBlock, alpha_offset::Int, beta_offset::Int)
    record = block.cross_hartree
    @inbounds for factor = 1:cross_hartree_rank(record)
        column += 1
        for j = 1:block.alpha.rank, i = 1:j
            prepared.cross_alpha[pair_index(alpha_offset+i, alpha_offset+j),
                column] = record.alpha[pair_index(i, j), factor]
        end
        for j = 1:block.beta.rank, i = 1:j
            prepared.cross_beta[pair_index(beta_offset+i, beta_offset+j),
                column] = record.beta[pair_index(i, j), factor]
        end
        prepared.cross_weight[column] = record.weight[factor]
    end
    column
end

function _append_physical_cross!(prepared::PreparedUHFCenter, column::Int,
        rows::UnitRange{Int}, alpha_offset::Int, beta_offset::Int,
        operator::UnitCellInteraction, scratch1, scratch2)
    @inbounds for beta_local in eachindex(rows), alpha_local in eachindex(rows)
        column += 1
        prepared.cross_alpha[pair_index(alpha_offset+alpha_local,
            alpha_offset+alpha_local), column] = 1.0
        prepared.cross_beta[pair_index(beta_offset+beta_local,
            beta_offset+beta_local), column] = 1.0
        prepared.cross_weight[column] = _interaction_entry(operator,
            rows[alpha_local], rows[beta_local], scratch1, scratch2)
    end
    column
end

function _append_component_cross!(prepared::PreparedUHFCenter, column::Int,
        first_component::Int, second_component::Int, alpha_offsets,
        beta_offsets, left::UHFOuterBlock, first_rows, second_rows,
        right::UHFOuterBlock, operator::UnitCellInteraction)
    channels = _maximum_channel_rank(operator)
    power = second_component - first_component - 1
    transfer = _channel_transfer!(prepared.transfer,
        prepared.transfer_scratch, operator, power)
    first_alpha = _center_component_rank(first_component, :alpha, left,
        first_rows, second_rows, right)
    second_alpha = _center_component_rank(second_component, :alpha, left,
        first_rows, second_rows, right)
    first_beta = _center_component_rank(first_component, :beta, left,
        first_rows, second_rows, right)
    second_beta = _center_component_rank(second_component, :beta, left,
        first_rows, second_rows, right)
    @inbounds for outgoing = 1:channels
        column += 1
        for j = 1:first_alpha, i = 1:j
            pair = pair_index(i, j); value = 0.0
            for incoming = 1:channels
                value = muladd(_center_component_feature(first_component,
                    :alpha, pair, incoming, false, left, first_rows,
                    second_rows, right, operator), transfer[incoming, outgoing],
                    value)
            end
            prepared.cross_alpha[pair_index(alpha_offsets[first_component]+i,
                alpha_offsets[first_component]+j), column] = value
        end
        for j = 1:second_beta, i = 1:j
            pair = pair_index(i, j)
            prepared.cross_beta[pair_index(beta_offsets[second_component]+i,
                beta_offsets[second_component]+j), column] =
                _center_component_feature(second_component, :beta, pair,
                    outgoing, true, left, first_rows, second_rows, right,
                    operator)
        end
        prepared.cross_weight[column] = 1.0

        column += 1
        for j = 1:second_alpha, i = 1:j
            pair = pair_index(i, j)
            prepared.cross_alpha[pair_index(alpha_offsets[second_component]+i,
                alpha_offsets[second_component]+j), column] =
                _center_component_feature(second_component, :alpha, pair,
                    outgoing, true, left, first_rows, second_rows, right,
                    operator)
        end
        for j = 1:first_beta, i = 1:j
            pair = pair_index(i, j); value = 0.0
            for incoming = 1:channels
                value = muladd(_center_component_feature(first_component,
                    :beta, pair, incoming, false, left, first_rows,
                    second_rows, right, operator), transfer[incoming, outgoing],
                    value)
            end
            prepared.cross_beta[pair_index(beta_offsets[first_component]+i,
                beta_offsets[first_component]+j), column] = value
        end
        prepared.cross_weight[column] = 1.0
    end
    column
end

function _prepare_center_cross!(prepared::PreparedUHFCenter,
        left::UHFOuterBlock, first_rows::UnitRange{Int},
        second_rows::UnitRange{Int}, right::UHFOuterBlock,
        operator::UnitCellInteraction)
    ranks_alpha = (left.alpha.rank, length(first_rows), length(second_rows),
        right.alpha.rank)
    ranks_beta = (left.beta.rank, length(first_rows), length(second_rows),
        right.beta.rank)
    offsets_alpha = (0, ranks_alpha[1], ranks_alpha[1]+ranks_alpha[2],
        ranks_alpha[1]+ranks_alpha[2]+ranks_alpha[3])
    offsets_beta = (0, ranks_beta[1], ranks_beta[1]+ranks_beta[2],
        ranks_beta[1]+ranks_beta[2]+ranks_beta[3])
    Ka, Kb = pair_dimension(sum(ranks_alpha)), pair_dimension(sum(ranks_beta))
    fill!(@view(prepared.cross_alpha[1:Ka, :]), 0.0)
    fill!(@view(prepared.cross_beta[1:Kb, :]), 0.0)
    column = _append_internal_cross!(prepared, 0, left, 0, 0)
    residual = size(operator.residual_transfer, 1)
    scratch1 = @view prepared.alpha.core_scratch[1:residual]
    scratch2 = @view prepared.beta.core_scratch[1:residual]
    column = _append_physical_cross!(prepared, column, first_rows,
        offsets_alpha[2], offsets_beta[2], operator, scratch1, scratch2)
    column = _append_physical_cross!(prepared, column, second_rows,
        offsets_alpha[3], offsets_beta[3], operator, scratch1, scratch2)
    column = _append_internal_cross!(prepared, column, right,
        offsets_alpha[4], offsets_beta[4])
    for second_component = 2:4, first_component = 1:second_component-1
        column = _append_component_cross!(prepared, column, first_component,
            second_component, offsets_alpha, offsets_beta, left, first_rows,
            second_rows, right, operator)
    end
    column <= size(prepared.cross_alpha, 2) || throw(DimensionMismatch(
        "UHF center cross-Hartree capacity is too small"))
    prepared.cross_rank = column
    nothing
end

function prepare_uhf_lifecycle_center!(prepared::PreparedUHFCenter,
        left::UHFOuterBlock, right::UHFOuterBlock, cell::Int,
        first_rows::UnitRange{Int}, second_rows::UnitRange{Int},
        one_body::BandedOneBody, operator::UnitCellInteraction)
    _validate_uhf_block(left, operator); _validate_uhf_block(right, operator)
    left.provenance == right.provenance || throw(ArgumentError(
        "UHF center provenance is inconsistent"))
    isempty(left.alpha.cells) || last(left.alpha.cells) + 1 == cell ||
        throw(ArgumentError("left UHF center interval is inconsistent"))
    isempty(right.alpha.cells) || cell + 2 == first(right.alpha.cells) ||
        throw(ArgumentError("right UHF center interval is inconsistent"))
    !isempty(first_rows) && !isempty(second_rows) &&
        last(first_rows) + 1 == first(second_rows) || throw(ArgumentError(
        "UHF center physical intervals must be adjacent"))
    ma = left.alpha.rank + length(first_rows) + length(second_rows) +
        right.alpha.rank
    mb = left.beta.rank + length(first_rows) + length(second_rows) +
        right.beta.rank
    ma <= prepared.alpha.maximum_rank && mb <= prepared.beta.maximum_rank ||
        throw(DimensionMismatch("UHF prepared-center capacity is too small"))
    required_columns = cross_hartree_rank(left.cross_hartree) +
        cross_hartree_rank(right.cross_hartree) + length(first_rows)^2 +
        length(second_rows)^2 + 12 * _maximum_channel_rank(operator)
    required_columns <= size(prepared.cross_alpha, 2) &&
        required_columns <= size(prepared.cross_beta, 2) ||
        throw(DimensionMismatch("UHF cross-Hartree capacity is too small"))
    _prepare_lifecycle_center!(prepared.alpha, left.alpha, right.alpha, cell,
        first_rows, second_rows, one_body, operator)
    _prepare_lifecycle_center!(prepared.beta, left.beta, right.beta, cell,
        first_rows, second_rows, one_body, operator)
    abs(prepared.alpha.constant_energy - prepared.beta.constant_energy) <=
        2048eps(Float64) * max(abs(prepared.alpha.constant_energy), 1.0) ||
        throw(ArgumentError("spin center core scalars disagree"))
    _prepare_center_cross!(prepared, left, first_rows, second_rows, right,
        operator)
    prepared.constant_energy = prepared.alpha.constant_energy
    prepared
end

function _uhf_center_preflight(prepared::PreparedUHFCenter, alpha_covariance,
        beta_covariance, alpha_occupied::Int, beta_occupied::Int)
    ma, mb = prepared.alpha.active_rank, prepared.beta.active_rank
    size(alpha_covariance) == (ma, ma) &&
        size(beta_covariance) == (mb, mb) || throw(DimensionMismatch(
        "UHF center covariance dimensions disagree"))
    0 <= alpha_occupied <= ma && 0 <= beta_occupied <= mb ||
        throw(DimensionMismatch("UHF center occupation exceeds a spin rank"))
    nothing
end

function uhf_center_energy_fock!(prepared::PreparedUHFCenter,
        alpha_covariance, beta_covariance, alpha_occupied::Int,
        beta_occupied::Int)
    _uhf_center_preflight(prepared, alpha_covariance, beta_covariance,
        alpha_occupied, beta_occupied)
    ma, mb = prepared.alpha.active_rank, prepared.beta.active_rank
    Da = @view prepared.alpha.density[1:ma, 1:ma]
    Db = @view prepared.beta.density[1:mb, 1:mb]
    copyto!(Da, alpha_covariance); copyto!(Db, beta_covariance)
    Ga = @view prepared.alpha.pair_field[1:pair_dimension(ma),
        1:pair_dimension(ma)]
    Gb = @view prepared.beta.pair_field[1:pair_dimension(mb),
        1:pair_dimension(mb)]
    _self_spin_fields!(prepared.alpha.coulomb, prepared.alpha.exchange, Da,
        Ga, ma)
    _self_spin_fields!(prepared.beta.coulomb, prepared.beta.exchange, Db,
        Gb, mb)
    pa = _pack_density!(prepared.alpha.pair_density, Da, ma)
    pb = _pack_density!(prepared.beta.pair_density, Db, mb)
    record = UHFCrossHartree(@view(prepared.cross_alpha[1:pair_dimension(ma),
            1:prepared.cross_rank]),
        @view(prepared.cross_beta[1:pair_dimension(mb), 1:prepared.cross_rank]),
        @view(prepared.cross_weight[1:prepared.cross_rank]), ma, mb)
    cross_energy = _cross_packed_fields!(
        @view(prepared.alpha.feature_scratch[1:pair_dimension(ma), 1]),
        @view(prepared.beta.feature_scratch[1:pair_dimension(mb), 1]),
        pa, pb, record, prepared.alpha_factor, prepared.beta_factor)
    h1a = @view prepared.alpha.h1[1:ma, 1:ma]
    h1b = @view prepared.beta.h1[1:mb, 1:mb]
    Fa = @view prepared.alpha.fock[1:ma, 1:ma]
    Fb = @view prepared.beta.fock[1:mb, 1:mb]
    @inbounds for j = 1:ma, i = 1:ma
        Fa[i, j] = h1a[i, j] + prepared.alpha.coulomb[i, j] -
            prepared.alpha.exchange[i, j] +
            prepared.alpha.feature_scratch[pair_index(i, j), 1] /
                pair_scale(i, j)
    end
    @inbounds for j = 1:mb, i = 1:mb
        Fb[i, j] = h1b[i, j] + prepared.beta.coulomb[i, j] -
            prepared.beta.exchange[i, j] +
            prepared.beta.feature_scratch[pair_index(i, j), 1] /
                pair_scale(i, j)
    end
    one = dot(h1a, Da) + dot(h1b, Db)
    hartree_alpha = 0.0; exchange_alpha = 0.0
    hartree_beta = 0.0; exchange_beta = 0.0
    @inbounds for j = 1:ma, i = 1:ma
        hartree_alpha = muladd(0.5Da[i, j], prepared.alpha.coulomb[i, j],
            hartree_alpha)
        exchange_alpha = muladd(-0.5Da[i, j],
            prepared.alpha.exchange[i, j], exchange_alpha)
    end
    @inbounds for j = 1:mb, i = 1:mb
        hartree_beta = muladd(0.5Db[i, j], prepared.beta.coulomb[i, j],
            hartree_beta)
        exchange_beta = muladd(-0.5Db[i, j],
            prepared.beta.exchange[i, j], exchange_beta)
    end
    hartree = hartree_alpha + hartree_beta + cross_energy
    exchange = exchange_alpha + exchange_beta
    (total=prepared.constant_energy + one + hartree + exchange,
        core=prepared.constant_energy, one_body=one,
        hartree=hartree, exchange=exchange,
        cross_hartree=cross_energy)
end

function apply_uhf_center_fock!(alpha_output::AbstractVecOrMat{Float64},
        beta_output::AbstractVecOrMat{Float64},
        prepared::PreparedUHFCenter, alpha_input::AbstractVecOrMat{Float64},
        beta_input::AbstractVecOrMat{Float64})
    ma, mb = prepared.alpha.active_rank, prepared.beta.active_rank
    size(alpha_input, 1) == ma && size(alpha_output) == size(alpha_input) &&
        size(beta_input, 1) == mb && size(beta_output) == size(beta_input) ||
        throw(DimensionMismatch("UHF center action dimensions disagree"))
    (Base.mightalias(alpha_output, alpha_input) ||
        Base.mightalias(beta_output, beta_input) ||
        Base.mightalias(alpha_output, beta_output)) && throw(ArgumentError(
        "UHF center action outputs must not alias inputs or each other"))
    mul!(alpha_output, @view(prepared.alpha.fock[1:ma, 1:ma]), alpha_input)
    mul!(beta_output, @view(prepared.beta.fock[1:mb, 1:mb]), beta_input)
    nothing
end
