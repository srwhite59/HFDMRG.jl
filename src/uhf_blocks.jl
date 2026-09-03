struct UHFCrossHartree{A<:AbstractMatrix{Float64},
        B<:AbstractMatrix{Float64},W<:AbstractVector{Float64}}
    alpha::A
    beta::B
    weight::W
    alpha_rank::Int
    beta_rank::Int
end

@inline cross_hartree_rank(record::UHFCrossHartree) = length(record.weight)

function _check_cross_hartree(record::UHFCrossHartree)
    pairs_alpha, pairs_beta = pair_dimension(record.alpha_rank),
        pair_dimension(record.beta_rank)
    size(record.alpha, 1) == pairs_alpha &&
        size(record.beta, 1) == pairs_beta &&
        size(record.alpha, 2) == size(record.beta, 2) ==
            length(record.weight) || throw(DimensionMismatch(
            "cross-Hartree factor dimensions disagree"))
    all(isfinite, record.alpha) && all(isfinite, record.beta) &&
        all(isfinite, record.weight) || throw(ArgumentError(
        "cross-Hartree factors must be finite"))
    nothing
end

mutable struct UHFOuterBlock
    alpha::RHFOuterBlock
    beta::RHFOuterBlock
    cross_hartree::UHFCrossHartree{Matrix{Float64},Matrix{Float64},
        Vector{Float64}}
    provenance::UInt
end

function _validate_uhf_block(block::UHFOuterBlock,
        operator::UnitCellInteraction)
    _validate_block(block.alpha, operator); _validate_block(block.beta, operator)
    block.alpha.cells == block.beta.cells && block.alpha.rows == block.beta.rows ||
        throw(ArgumentError("spin blocks describe different physical intervals"))
    block.alpha.reflected == block.beta.reflected || throw(ArgumentError(
        "spin block orientations disagree"))
    block.alpha.provenance == block.beta.provenance == block.provenance ||
        throw(ArgumentError("spin block provenance disagrees"))
    block.alpha.left_core === block.beta.left_core &&
        block.alpha.right_core === block.beta.right_core || throw(ArgumentError(
        "UHF blocks must share one total-density Hartree connector"))
    block.alpha.constant_energy == block.beta.constant_energy ||
        throw(ArgumentError("spin blocks disagree on the core scalar"))
    record = block.cross_hartree
    _check_cross_hartree(record)
    record.alpha_rank == block.alpha.rank &&
        record.beta_rank == block.beta.rank || throw(DimensionMismatch(
        "cross-Hartree and spin ranks disagree"))
    cross_hartree_rank(record) <= min(pair_dimension(block.alpha.rank),
        pair_dimension(block.beta.rank)) || throw(DimensionMismatch(
        "cross-Hartree rank exceeds a spin response space"))
    nothing
end

mutable struct UHFBlockBuildWorkspace
    alpha::BlockBuildWorkspace
    beta::BlockBuildWorkspace
    maximum_alpha::Int
    maximum_beta::Int
    maximum_columns::Int
    alpha_density::Matrix{Float64}
    beta_density::Matrix{Float64}
    alpha_j::Matrix{Float64}
    alpha_k::Matrix{Float64}
    beta_j::Matrix{Float64}
    beta_k::Matrix{Float64}
    parent_alpha::Matrix{Float64}
    parent_beta::Matrix{Float64}
    transformed_alpha::Matrix{Float64}
    transformed_beta::Matrix{Float64}
    weights::Vector{Float64}
    factor_alpha::Vector{Float64}
    factor_beta::Vector{Float64}
    factor_work_alpha::Vector{Float64}
    factor_work_beta::Vector{Float64}
end

function UHFBlockBuildWorkspace(maximum_alpha::Int, maximum_beta::Int,
        channel_rank::Int)
    min(maximum_alpha, maximum_beta) >= 1 || throw(ArgumentError(
        "spin build capacities must be positive"))
    Ka, Kb = pair_dimension(maximum_alpha), pair_dimension(maximum_beta)
    columns = max(2min(Ka, Kb) + 2channel_rank, UNIT_CELL_WIDTH^2)
    UHFBlockBuildWorkspace(BlockBuildWorkspace(maximum_alpha, channel_rank),
        BlockBuildWorkspace(maximum_beta, channel_rank), maximum_alpha,
        maximum_beta, columns, zeros(maximum_alpha, maximum_alpha),
        zeros(maximum_beta, maximum_beta), zeros(maximum_alpha, maximum_alpha),
        zeros(maximum_alpha, maximum_alpha), zeros(maximum_beta, maximum_beta),
        zeros(maximum_beta, maximum_beta), zeros(Ka, columns),
        zeros(Kb, columns), zeros(Ka, columns), zeros(Kb, columns),
        zeros(columns), zeros(columns), zeros(columns), zeros(columns),
        zeros(columns))
end

function _compress_cross_hartree(alpha::Matrix{Float64},
        beta::Matrix{Float64}, weight::Vector{Float64}, ma::Int, mb::Int)
    Ka, columns = size(alpha); Kb = size(beta, 1)
    size(beta, 2) == columns == length(weight) || throw(DimensionMismatch(
        "cross-Hartree candidate columns disagree"))
    Ka == pair_dimension(ma) && Kb == pair_dimension(mb) ||
        throw(DimensionMismatch("cross-Hartree candidate ranks disagree"))
    if columns == 0 || Ka == 0 || Kb == 0
        return UHFCrossHartree(zeros(Ka, 0), zeros(Kb, 0), Float64[], ma, mb)
    end
    @inbounds for column = 1:columns, row = 1:Ka
        alpha[row, column] *= weight[column]
    end
    first = svd!(alpha; full=false)
    # An intermediate singular value is not itself a physical cross-Hartree
    # weight: its partner in `beta` may carry the reciprocal scale.  Retain
    # the complete numerical range until both oriented factors have been
    # combined, or a harmless factor gauge can become a physical truncation.
    first_rank = length(first.S)
    first_rank == 0 && return UHFCrossHartree(zeros(Ka, 0), zeros(Kb, 0),
        Float64[], ma, mb)
    carried = beta * @view(first.V[:, 1:first_rank])
    @inbounds for column = 1:first_rank, row = 1:Kb
        carried[row, column] *= first.S[column]
    end
    second = svd!(carried; full=false)
    rank = length(second.S)
    output_alpha = Matrix{Float64}(undef, Ka, rank)
    rank > 0 && mul!(output_alpha, @view(first.U[:, 1:first_rank]),
        @view(second.V[:, 1:rank]))
    UHFCrossHartree(output_alpha, Matrix(@view(second.U[:, 1:rank])),
        Vector(@view(second.S[1:rank])), ma, mb)
end

function _seed_cross_hartree(operator::UnitCellInteraction,
        rows::UnitRange{Int})
    width = length(rows); pairs = pair_dimension(width); columns = width^2
    alpha = zeros(pairs, columns); beta = zeros(pairs, columns)
    weight = Vector{Float64}(undef, columns)
    scratch1 = zeros(size(operator.residual_transfer, 1))
    scratch2 = similar(scratch1)
    column = 0
    @inbounds for second in eachindex(rows), first in eachindex(rows)
        column += 1
        alpha[pair_index(first, first), column] = 1.0
        beta[pair_index(second, second), column] = 1.0
        weight[column] = _interaction_entry(operator, rows[first],
            rows[second], scratch1, scratch2)
    end
    _compress_cross_hartree(alpha, beta, weight, width, width)
end

function seed_uhf_interval_block(one_body::BandedOneBody,
        operator::UnitCellInteraction, cell::Int, rows::UnitRange{Int};
        provenance::UInt=UInt(0), reflected::Bool=false)
    alpha = seed_interval_block(one_body, operator, cell, rows; provenance,
        reflected)
    beta = seed_interval_block(one_body, operator, cell, rows; provenance,
        reflected)
    channels = _maximum_channel_rank(operator)
    left_core = zeros(channels); right_core = zeros(channels)
    alpha.left_core = left_core; beta.left_core = left_core
    alpha.right_core = right_core; beta.right_core = right_core
    block = UHFOuterBlock(alpha, beta, _seed_cross_hartree(operator, rows),
        provenance)
    _validate_uhf_block(block, operator)
    block
end

function _empty_uhf_block(operator::UnitCellInteraction, provenance::UInt)
    alpha = _empty_outer_block(operator, provenance)
    beta = _empty_outer_block(operator, provenance)
    left_core = alpha.left_core; right_core = alpha.right_core
    beta.left_core = left_core; beta.right_core = right_core
    block = UHFOuterBlock(alpha, beta,
        UHFCrossHartree(zeros(0, 0), zeros(0, 0), Float64[], 0, 0),
        provenance)
    _validate_uhf_block(block, operator)
    block
end

function _pack_density!(output, density, rank::Int)
    @inbounds for j = 1:rank, i = 1:j
        output[pair_index(i, j)] = pair_scale(i, j) * density[i, j]
    end
    @view output[1:pair_dimension(rank)]
end

function _self_spin_fields!(jfield, kfield, density, pair_field, rank::Int)
    fill!(@view(jfield[1:rank, 1:rank]), 0.0)
    fill!(@view(kfield[1:rank, 1:rank]), 0.0)
    @inbounds for b = 1:rank, a = 1:rank, d = 1:rank, c = 1:rank
        value = density[c, d]
        jfield[a, b] = muladd(value,
            _raw_integral(pair_field, a, b, c, d), jfield[a, b])
        kfield[a, b] = muladd(value,
            _raw_integral(pair_field, a, c, b, d), kfield[a, b])
    end
    nothing
end

function _cross_packed_fields!(alpha_output, beta_output,
        alpha_density, beta_density, record::UHFCrossHartree,
        alpha_factor, beta_factor)
    rank = cross_hartree_rank(record)
    fill!(@view(alpha_output[1:pair_dimension(record.alpha_rank)]), 0.0)
    fill!(@view(beta_output[1:pair_dimension(record.beta_rank)]), 0.0)
    @inbounds for factor = 1:rank
        avalue = 0.0; bvalue = 0.0
        for pair in eachindex(alpha_density)
            avalue = muladd(record.alpha[pair, factor], alpha_density[pair],
                avalue)
        end
        for pair in eachindex(beta_density)
            bvalue = muladd(record.beta[pair, factor], beta_density[pair],
                bvalue)
        end
        alpha_factor[factor] = avalue
        beta_factor[factor] = bvalue
        for pair in eachindex(alpha_output)
            alpha_output[pair] = muladd(record.alpha[pair, factor],
                record.weight[factor] * bvalue, alpha_output[pair])
        end
        for pair in eachindex(beta_output)
            beta_output[pair] = muladd(record.beta[pair, factor],
                record.weight[factor] * avalue, beta_output[pair])
        end
    end
    energy = 0.0
    @inbounds for factor = 1:rank
        energy = muladd(alpha_factor[factor] * record.weight[factor],
            beta_factor[factor], energy)
    end
    energy
end

function _append_cross_parent!(workspace::UHFBlockBuildWorkspace,
        left::UHFOuterBlock, right::UHFOuterBlock, old_alpha::Int,
        old_beta::Int, channels::Int)
    fill!(@view(workspace.parent_alpha[1:pair_dimension(old_alpha), :]), 0.0)
    fill!(@view(workspace.parent_beta[1:pair_dimension(old_beta), :]), 0.0)
    column = 0
    for (child, alpha_offset, beta_offset) in
            ((left, 0, 0), (right, left.alpha.rank, left.beta.rank))
        record = child.cross_hartree
        for factor = 1:cross_hartree_rank(record)
            column += 1
            @inbounds for j = 1:child.alpha.rank, i = 1:j
                workspace.parent_alpha[pair_index(alpha_offset+i,
                    alpha_offset+j), column] = record.alpha[
                    pair_index(i, j), factor]
            end
            @inbounds for j = 1:child.beta.rank, i = 1:j
                workspace.parent_beta[pair_index(beta_offset+i,
                    beta_offset+j), column] = record.beta[
                    pair_index(i, j), factor]
            end
            workspace.weights[column] = record.weight[factor]
        end
    end
    for channel = 1:channels
        column += 1
        @inbounds for j = 1:left.alpha.rank, i = 1:j
            workspace.parent_alpha[pair_index(i, j), column] =
                left.alpha.right_feature[pair_index(i, j), channel]
        end
        @inbounds for j = 1:right.beta.rank, i = 1:j
            workspace.parent_beta[pair_index(left.beta.rank+i,
                left.beta.rank+j), column] =
                right.beta.left_feature[pair_index(i, j), channel]
        end
        workspace.weights[column] = 1.0
        column += 1
        @inbounds for j = 1:right.alpha.rank, i = 1:j
            workspace.parent_alpha[pair_index(left.alpha.rank+i,
                left.alpha.rank+j), column] =
                right.alpha.left_feature[pair_index(i, j), channel]
        end
        @inbounds for j = 1:left.beta.rank, i = 1:j
            workspace.parent_beta[pair_index(i, j), column] =
                left.beta.right_feature[pair_index(i, j), channel]
        end
        workspace.weights[column] = 1.0
    end
    column <= workspace.maximum_columns || throw(DimensionMismatch(
        "cross-Hartree build capacity is too small"))
    column
end

function _spin_publication(workspace::BlockBuildWorkspace,
        left::RHFOuterBlock, right::RHFOuterBlock, link::RHFStateLink,
        one_body::BandedOneBody, constant::Float64, left_core, right_core)
    old_rank = left.rank + right.rank; new_rank = link.active_rank
    old_pairs = pair_dimension(old_rank); new_pairs = pair_dimension(new_rank)
    transform = _pair_transform!(workspace.pair_map, link.close_map,
        old_rank, new_rank)
    h1 = Matrix{Float64}(undef, new_rank, new_rank)
    scratch_h1 = @view workspace.h1_scratch[1:old_rank, 1:new_rank]
    mul!(scratch_h1, @view(workspace.h1[1:old_rank, 1:old_rank]),
        link.close_map)
    mul!(h1, transpose(link.close_map), scratch_h1)
    scratch_pair = @view workspace.pair_scratch[1:old_pairs, 1:new_pairs]
    mul!(scratch_pair, @view(workspace.pair_field[1:old_pairs, 1:old_pairs]),
        transform)
    pair_field = Matrix{Float64}(undef, new_pairs, new_pairs)
    mul!(pair_field, transpose(transform), scratch_pair)
    channels = size(left.left_feature, 2)
    left_feature = Matrix{Float64}(undef, new_pairs, channels)
    right_feature = similar(left_feature)
    mul!(left_feature, transpose(transform),
        @view(workspace.feature_left[1:old_pairs, 1:channels]))
    mul!(right_feature, transpose(transform),
        @view(workspace.feature_right[1:old_pairs, 1:channels]))
    rows = isempty(left.rows) ? right.rows : isempty(right.rows) ? left.rows :
        first(left.rows):last(right.rows)
    edge_width = min(length(rows), max(h1_bandwidth(one_body), 1))
    left_edge = zeros(edge_width, new_rank)
    right_edge = zeros(edge_width, new_rank)
    for (child, offset) in ((left, 0), (right, left.rank))
        @inbounds for edge = 1:edge_width, column = 1:new_rank,
                old = 1:child.rank
            left_edge[edge, column] = muladd(_edge_value(child,
                first(rows)+edge-1, old, true),
                link.close_map[offset+old, column], left_edge[edge, column])
            right_edge[edge, column] = muladd(_edge_value(child,
                last(rows)-edge_width+edge, old, false),
                link.close_map[offset+old, column], right_edge[edge, column])
        end
    end
    cells = isempty(left.cells) ? right.cells : isempty(right.cells) ?
        left.cells : first(left.cells):last(right.cells)
    RHFOuterBlock(cells, rows, new_rank,
        left.completed_count + right.completed_count + link.selection.completed,
        constant, h1, pair_field, left_feature, right_feature, left_core,
        right_core, left_edge, right_edge, link, false, left.provenance)
end

function build_uhf_outer_block(left::UHFOuterBlock, right::UHFOuterBlock,
        alpha_link::RHFStateLink, beta_link::RHFStateLink,
        one_body::BandedOneBody, operator::UnitCellInteraction,
        workspace::UHFBlockBuildWorkspace)
    _validate_uhf_block(left, operator); _validate_uhf_block(right, operator)
    left.provenance == right.provenance || throw(ArgumentError(
        "UHF block provenance is inconsistent"))
    olda, pairsa, consta, _ = _fill_merged!(workspace.alpha, left.alpha,
        right.alpha, one_body, operator)
    oldb, pairsb, constb, _ = _fill_merged!(workspace.beta, left.beta,
        right.beta, one_body, operator)
    abs(consta - constb) <= 1024eps(Float64) * max(abs(consta), 1.0) ||
        throw(ArgumentError("shared-core scalar differs between spin bases"))
    olda == alpha_link.old_rank + alpha_link.physical_width &&
        oldb == beta_link.old_rank + beta_link.physical_width ||
        throw(DimensionMismatch("UHF links do not describe the parent blocks"))
    channels = _maximum_channel_rank(operator)
    columns = _append_cross_parent!(workspace, left, right, olda, oldb,
        channels)
    Ca = @view workspace.alpha_density[1:olda, 1:olda]
    Cb = @view workspace.beta_density[1:oldb, 1:oldb]
    mul!(Ca, alpha_link.completed_map, transpose(alpha_link.completed_map))
    mul!(Cb, beta_link.completed_map, transpose(beta_link.completed_map))
    Ga = @view workspace.alpha.pair_field[1:pairsa, 1:pairsa]
    Gb = @view workspace.beta.pair_field[1:pairsb, 1:pairsb]
    _self_spin_fields!(workspace.alpha_j, workspace.alpha_k, Ca, Ga, olda)
    _self_spin_fields!(workspace.beta_j, workspace.beta_k, Cb, Gb, oldb)
    pa = _pack_density!(workspace.factor_alpha, Ca, olda)
    pb = _pack_density!(workspace.factor_beta, Cb, oldb)
    cross = UHFCrossHartree(
        @view(workspace.parent_alpha[1:pairsa, 1:columns]),
        @view(workspace.parent_beta[1:pairsb, 1:columns]),
        @view(workspace.weights[1:columns]), olda, oldb)
    cross_energy = _cross_packed_fields!(
        @view(workspace.transformed_alpha[1:pairsa, 1]),
        @view(workspace.transformed_beta[1:pairsb, 1]), pa, pb, cross,
        workspace.factor_work_alpha, workspace.factor_work_beta)
    h1a = @view workspace.alpha.h1[1:olda, 1:olda]
    h1b = @view workspace.beta.h1[1:oldb, 1:oldb]
    constant = consta + dot(h1a, Ca) + dot(h1b, Cb)
    constant += 0.5dot(Ca, @view(workspace.alpha_j[1:olda, 1:olda]) -
        @view(workspace.alpha_k[1:olda, 1:olda]))
    constant += 0.5dot(Cb, @view(workspace.beta_j[1:oldb, 1:oldb]) -
        @view(workspace.beta_k[1:oldb, 1:oldb])) + cross_energy
    @inbounds for j = 1:olda, i = 1:olda
        pair = pair_index(i, j)
        h1a[i, j] += workspace.alpha_j[i, j] - workspace.alpha_k[i, j] +
            workspace.transformed_alpha[pair, 1] / pair_scale(i, j)
    end
    @inbounds for j = 1:oldb, i = 1:oldb
        pair = pair_index(i, j)
        h1b[i, j] += workspace.beta_j[i, j] - workspace.beta_k[i, j] +
            workspace.transformed_beta[pair, 1] / pair_scale(i, j)
    end
    left_core = copy(@view workspace.alpha.left_core[1:channels])
    right_core = copy(@view workspace.alpha.right_core[1:channels])
    @inbounds for channel = 1:channels
        for pair in eachindex(pa)
            left_core[channel] = muladd(0.5pa[pair],
                workspace.alpha.feature_left[pair, channel], left_core[channel])
            right_core[channel] = muladd(0.5pa[pair],
                workspace.alpha.feature_right[pair, channel], right_core[channel])
        end
        for pair in eachindex(pb)
            left_core[channel] = muladd(0.5pb[pair],
                workspace.beta.feature_left[pair, channel], left_core[channel])
            right_core[channel] = muladd(0.5pb[pair],
                workspace.beta.feature_right[pair, channel], right_core[channel])
        end
    end
    alpha = _spin_publication(workspace.alpha, left.alpha, right.alpha,
        alpha_link, one_body, constant, left_core, right_core)
    beta = _spin_publication(workspace.beta, left.beta, right.beta, beta_link,
        one_body, constant, left_core, right_core)
    Ta = _pair_transform!(workspace.alpha.pair_map, alpha_link.close_map,
        olda, alpha.rank)
    Tb = _pair_transform!(workspace.beta.pair_map, beta_link.close_map,
        oldb, beta.rank)
    transformed_a = Matrix{Float64}(undef, pair_dimension(alpha.rank), columns)
    transformed_b = Matrix{Float64}(undef, pair_dimension(beta.rank), columns)
    mul!(transformed_a, transpose(Ta),
        @view(workspace.parent_alpha[1:pairsa, 1:columns]))
    mul!(transformed_b, transpose(Tb),
        @view(workspace.parent_beta[1:pairsb, 1:columns]))
    record = _compress_cross_hartree(transformed_a, transformed_b,
        Vector(@view(workspace.weights[1:columns])), alpha.rank, beta.rank)
    block = UHFOuterBlock(alpha, beta, record, left.provenance)
    _validate_uhf_block(block, operator)
    block
end

function reflect_uhf_block(block::UHFOuterBlock, total_cells::Int,
        total_rows::Int, operator::UnitCellInteraction)
    alpha = reflect_block(block.alpha, total_cells, total_rows)
    beta = reflect_block(block.beta, total_cells, total_rows)
    left_core = alpha.left_core; right_core = alpha.right_core
    beta.left_core = left_core; beta.right_core = right_core
    record = UHFCrossHartree(copy(block.cross_hartree.alpha),
        copy(block.cross_hartree.beta), copy(block.cross_hartree.weight),
        block.alpha.rank, block.beta.rank)
    reflected = UHFOuterBlock(alpha, beta, record, block.provenance)
    _validate_uhf_block(reflected, operator)
    reflected
end

function uhf_block_storage_bytes(block::UHFOuterBlock)
    block_storage_bytes(block.alpha) + block_storage_bytes(block.beta) -
        sizeof(block.alpha.left_core) - sizeof(block.alpha.right_core) +
        sizeof(block.cross_hartree.alpha) + sizeof(block.cross_hartree.beta) +
        sizeof(block.cross_hartree.weight)
end
