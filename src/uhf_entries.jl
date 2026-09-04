"""Separated compact-UHF entry records and reusable donor-derived builder."""
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
    Ka, Kb = pair_dimension(record.alpha_rank), pair_dimension(record.beta_rank)
    size(record.alpha) == (Ka, length(record.weight)) &&
        size(record.beta) == (Kb, length(record.weight)) ||
        throw(DimensionMismatch("cross-Hartree factor dimensions disagree"))
    all(isfinite, record.alpha) && all(isfinite, record.beta) &&
        all(isfinite, record.weight) || throw(ArgumentError(
        "cross-Hartree factors must be finite"))
    nothing
end

function _pack_density!(output, density, rank::Int)
    @inbounds for j = 1:rank, i = 1:j
        output[pair_index(i, j)] = pair_scale(i, j) * density[i, j]
    end
    @view output[1:pair_dimension(rank)]
end

function _self_spin_fields!(jfield, kfield, density, internal, rank::Int)
    fill!(jfield, 0.0); fill!(kfield, 0.0)
    @inbounds for b = 1:rank, a = 1:rank, d = 1:rank, c = 1:rank
        value = density[c, d]
        jfield[a, b] = muladd(value,
            _raw_integral(internal, a, b, c, d), jfield[a, b])
        kfield[a, b] = muladd(value,
            _raw_integral(internal, a, c, b, d), kfield[a, b])
    end
    nothing
end

function _cross_packed_fields!(alpha_output, beta_output, alpha_density,
        beta_density, record::UHFCrossHartree, alpha_factor, beta_factor)
    rank = cross_hartree_rank(record)
    fill!(alpha_output, 0.0); fill!(beta_output, 0.0)
    energy = 0.0
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
        alpha_factor[factor] = avalue; beta_factor[factor] = bvalue
        energy = muladd(avalue * record.weight[factor], bvalue, energy)
        for pair in eachindex(alpha_output)
            alpha_output[pair] = muladd(record.alpha[pair, factor],
                record.weight[factor] * bvalue, alpha_output[pair])
        end
        for pair in eachindex(beta_output)
            beta_output[pair] = muladd(record.beta[pair, factor],
                record.weight[factor] * avalue, beta_output[pair])
        end
    end
    energy
end

function _compress_cross_hartree(alpha::Matrix{Float64},
        beta::Matrix{Float64}, weight::Vector{Float64}, ma::Int, mb::Int)
    Ka, columns = size(alpha); Kb = size(beta, 1)
    size(beta, 2) == columns == length(weight) &&
        Ka == pair_dimension(ma) && Kb == pair_dimension(mb) ||
        throw(DimensionMismatch("cross-Hartree candidate dimensions disagree"))
    columns == 0 && return UHFCrossHartree(zeros(Ka, 0), zeros(Kb, 0),
        Float64[], ma, mb)
    @inbounds for column = 1:columns, row = 1:Ka
        alpha[row, column] *= weight[column]
    end
    first = svd!(alpha; full=false); first_rank = length(first.S)
    carried = beta * @view(first.V[:, 1:first_rank])
    @inbounds for column = 1:first_rank, row = 1:Kb
        carried[row, column] *= first.S[column]
    end
    second = svd!(carried; full=false); rank = length(second.S)
    output_alpha = Matrix{Float64}(undef, Ka, rank)
    rank > 0 && mul!(output_alpha, @view(first.U[:, 1:first_rank]),
        @view(second.V[:, 1:rank]))
    UHFCrossHartree(output_alpha, Matrix(@view(second.U[:, 1:rank])),
        Vector(@view(second.S[1:rank])), ma, mb)
end
mutable struct UHFEntrySpinBlock
    cells::UnitRange{Int}
    rows::UnitRange{Int}
    rank::Int
    completed_count::Int
    empty_count::Int
    one_body_energy::Float64
    h1::Matrix{Float64}
    internal::Matrix{Float64}
    connector::Matrix{Float64}
    core_field::Matrix{Float64}
    left_edge::Matrix{Float64}
    right_edge::Matrix{Float64}
    link::RHFStateLink
    reflected::Bool
    provenance::UInt
end

mutable struct UHFEntryBlock
    alpha::UHFEntrySpinBlock
    beta::UHFEntrySpinBlock
    cross_hartree::UHFCrossHartree{Matrix{Float64},Matrix{Float64},
        Vector{Float64}}
    hartree_channel::Vector{Float64}
    core_energy::Float64
    provenance::UInt
end

function _validate_uhf_entry(block::UHFEntryBlock,
        operator::UnitCellInteraction)
    channels = _maximum_channel_rank(operator)
    for spin in (block.alpha, block.beta)
        K = pair_dimension(spin.rank)
        size(spin.h1) == size(spin.core_field) == (spin.rank, spin.rank) &&
            size(spin.internal) == (K, K) &&
            size(spin.connector) == (K, channels) || throw(DimensionMismatch(
            "separated UHF entry dimensions disagree"))
        size(spin.left_edge, 2) == size(spin.right_edge, 2) == spin.rank &&
            spin.link.active_rank == spin.rank || throw(DimensionMismatch(
            "separated UHF entry link or edge dimensions disagree"))
        min(spin.completed_count, spin.empty_count) >= 0 &&
            isfinite(spin.one_body_energy) || throw(ArgumentError(
            "separated UHF entry sector data is invalid"))
        all(isfinite, spin.h1) && all(isfinite, spin.internal) &&
            all(isfinite, spin.connector) && all(isfinite, spin.core_field) ||
            throw(ArgumentError("separated UHF entry must be finite"))
    end
    block.alpha.cells == block.beta.cells && block.alpha.rows == block.beta.rows &&
        block.alpha.reflected == block.beta.reflected || throw(ArgumentError(
        "separated spin entries disagree on geometry"))
    block.alpha.provenance == block.beta.provenance == block.provenance ||
        throw(ArgumentError("separated spin provenance disagrees"))
    length(block.hartree_channel) == channels &&
        all(isfinite, block.hartree_channel) && isfinite(block.core_energy) ||
        throw(ArgumentError("separated UHF shared core is invalid"))
    _check_cross_hartree(block.cross_hartree)
    block.cross_hartree.alpha_rank == block.alpha.rank &&
        block.cross_hartree.beta_rank == block.beta.rank ||
        throw(DimensionMismatch("cross-Hartree and spin ranks disagree"))
    nothing
end

function _empty_uhf_entry(operator::UnitCellInteraction, provenance::UInt)
    channels = _maximum_channel_rank(operator)
    selection = StateSelection(0, 0, 0, 0.0, :structural_boundary)
    link = RHFStateLink(0, 0, zeros(0, 0), zeros(0, 0), selection)
    spin() = UHFEntrySpinBlock(1:0, 1:0, 0, 0, 0, 0.0, zeros(0, 0),
        zeros(0, 0), zeros(0, channels), zeros(0, 0), zeros(0, 0),
        zeros(0, 0), link, false, provenance)
    block = UHFEntryBlock(spin(), spin(),
        UHFCrossHartree(zeros(0, 0), zeros(0, 0), Float64[], 0, 0),
        zeros(channels), 0.0, provenance)
    _validate_uhf_entry(block, operator)
    block
end

@inline function _uhf_entry_edge(block::UHFEntrySpinBlock, physical::Int,
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

mutable struct _IndependentSpinParentBuildCapacity
    maximum_parent::Int
    internal::Matrix{Float64}
    connector::Matrix{Float64}
    transformed_connector::Matrix{Float64}
    h1::Matrix{Float64}
    inherited::Matrix{Float64}
    completed::Matrix{Float64}
    coulomb::Matrix{Float64}
    exchange::Matrix{Float64}
    pair_map::Matrix{Float64}
    pair_scratch::Matrix{Float64}
    dense_scratch::Matrix{Float64}
    packed::Vector{Float64}
    cross_field::Vector{Float64}
end

function _IndependentSpinParentBuildCapacity(maximum_old::Int,
        maximum_width::Int, channels::Int)
    parent = maximum_old + maximum_width; K = pair_dimension(parent)
    _IndependentSpinParentBuildCapacity(parent, zeros(K, K),
        zeros(K, channels), zeros(K, channels), zeros(parent, parent),
        zeros(parent, parent), zeros(parent, parent), zeros(parent, parent),
        zeros(parent, parent), zeros(K, K), zeros(K, K),
        zeros(parent, parent), zeros(K), zeros(K))
end

mutable struct _IndependentSpinEntryBuildWorkspace
    alpha::_IndependentSpinParentBuildCapacity
    beta::_IndependentSpinParentBuildCapacity
    maximum_width::Int
    maximum_columns::Int
    local_interaction::Matrix{Float64}
    transfer::Matrix{Float64}
    transfer_scratch::Matrix{Float64}
    transferred_hartree::Vector{Float64}
    parent_alpha::Matrix{Float64}
    parent_beta::Matrix{Float64}
    transformed_alpha::Matrix{Float64}
    transformed_beta::Matrix{Float64}
    weights::Vector{Float64}
    factor_alpha::Vector{Float64}
    factor_beta::Vector{Float64}
    interaction_scratch_first::Vector{Float64}
    interaction_scratch_second::Vector{Float64}
end

function _IndependentSpinEntryBuildWorkspace(maximum_alpha::Int,
        maximum_beta::Int, maximum_width::Int, operator::UnitCellInteraction)
    channels = _maximum_channel_rank(operator)
    pa, pb = maximum_alpha + maximum_width, maximum_beta + maximum_width
    Ka, Kb = pair_dimension(pa), pair_dimension(pb)
    columns = min(pair_dimension(maximum_alpha), pair_dimension(maximum_beta)) +
        3maximum_width
    residual = size(operator.residual_transfer, 1)
    _IndependentSpinEntryBuildWorkspace(
        _IndependentSpinParentBuildCapacity(maximum_alpha, maximum_width,
            channels),
        _IndependentSpinParentBuildCapacity(maximum_beta, maximum_width,
            channels), maximum_width, columns, zeros(maximum_width,
            maximum_width), zeros(channels, channels), zeros(channels,
            channels), zeros(channels), zeros(Ka, columns),
        zeros(Kb, columns), zeros(Ka, columns), zeros(Kb, columns),
        zeros(columns), zeros(columns), zeros(columns), zeros(residual),
        zeros(residual))
end

@inline _entry_old_position(index::Int, width::Int, left::Bool) =
    left ? index : width + index
@inline _entry_site_position(index::Int, old_rank::Int, left::Bool) =
    left ? old_rank + index : index

function _entry_local_interaction!(workspace, operator, rows)
    values = @view workspace.local_interaction[1:length(rows), 1:length(rows)]
    @inbounds for j in eachindex(rows), i in eachindex(rows)
        values[i, j] = _interaction_entry(operator, rows[i], rows[j],
            workspace.interaction_scratch_first,
            workspace.interaction_scratch_second)
    end
    values
end

function _entry_transfer_connector!(output, input, workspace, operator,
        left::Bool)
    T = _channel_transfer!(workspace.transfer, workspace.transfer_scratch,
        operator, 1)
    left ? mul!(output, input, T) : mul!(output, input, transpose(T))
    output
end

function _bare_spin_parent!(capacity::_IndependentSpinParentBuildCapacity,
        old::Union{Nothing,UHFEntrySpinBlock}, rows, operator, left::Bool,
        workspace)
    old_rank = old === nothing ? 0 : old.rank
    width = length(rows); parent = old_rank + width
    parent <= capacity.maximum_parent || throw(DimensionMismatch(
        "independent-spin parent capacity"))
    K = pair_dimension(parent); channels = _maximum_channel_rank(operator)
    internal = @view capacity.internal[1:K, 1:K]
    connector = @view capacity.connector[1:K, 1:channels]
    transformed = @view capacity.transformed_connector[1:K, 1:channels]
    fill!(internal, 0.0); fill!(connector, 0.0); fill!(transformed, 0.0)
    if old !== nothing && old_rank > 0
        @inbounds for j = 1:old_rank, i = 1:j
            source = pair_index(i, j)
            target = pair_index(_entry_old_position(i, width, left),
                _entry_old_position(j, width, left))
            copyto!(@view(transformed[target, :]), @view(old.connector[source, :]))
            for l = 1:old_rank, k = 1:l
                source2 = pair_index(k, l)
                target2 = pair_index(_entry_old_position(k, width, left),
                    _entry_old_position(l, width, left))
                internal[target, target2] = old.internal[source, source2]
            end
        end
        _entry_transfer_connector!(connector, transformed, workspace,
            operator, left)
    end
    @inbounds for s in eachindex(rows)
        ps = _entry_site_position(s, old_rank, left); qs = pair_index(ps, ps)
        phase = _phase(operator, rows[s])
        for t in eachindex(rows)
            pt = _entry_site_position(t, old_rank, left)
            internal[qs, pair_index(pt, pt)] = workspace.local_interaction[s, t]
        end
        for channel = 1:channels
            connector[qs, channel] += left ?
                operator.phase_source[phase, channel] :
                operator.phase_sink[channel, phase]
        end
        old === nothing && continue
        for j = 1:old_rank, i = 1:j
            source = pair_index(i, j)
            target = pair_index(_entry_old_position(i, width, left),
                _entry_old_position(j, width, left))
            value = 0.0
            for channel = 1:channels
                value = muladd(old.connector[source, channel], left ?
                    operator.phase_sink[channel, phase] :
                    operator.phase_source[phase, channel], value)
            end
            internal[target, qs] = value; internal[qs, target] = value
        end
    end
    (internal=internal, connector=connector, transformed=transformed,
        old_rank=old_rank, parent=parent)
end

function _entry_parent_h1!(capacity, old, rows, one_body, left::Bool)
    old_rank = old === nothing ? 0 : old.rank
    width = length(rows); parent = old_rank + width
    h1 = @view capacity.h1[1:parent, 1:parent]; fill!(h1, 0.0)
    offset = left ? 0 : width
    old_rank > 0 && copyto!(@view(h1[offset+1:offset+old_rank,
        offset+1:offset+old_rank]), old.h1)
    @inbounds for j in eachindex(rows), i in eachindex(rows)
        h1[_entry_site_position(i, old_rank, left),
            _entry_site_position(j, old_rank, left)] =
            h1_entry(one_body, rows[i], rows[j])
    end
    if old !== nothing && old_rank > 0
        boundary = min(h1_bandwidth(one_body), left ?
            size(old.right_edge, 1) : size(old.left_edge, 1))
        @inbounds for a = 1:old_rank, b in eachindex(rows)
            value = 0.0
            for edge = 1:boundary
                site = left ? last(old.rows)-boundary+edge :
                    first(old.rows)+edge-1
                coefficient = left ? old.right_edge[edge, a] :
                    old.left_edge[edge, a]
                value = muladd(coefficient,
                    h1_entry(one_body, site, rows[b]), value)
            end
            po = _entry_old_position(a, width, left)
            pp = _entry_site_position(b, old_rank, left)
            h1[po, pp] = value; h1[pp, po] = value
        end
    end
    h1
end

function _fill_parent_inherited!(capacity, old, rows, operator, left::Bool,
        alpha::Bool)
    old_rank = old === nothing ? 0 : (alpha ? old.alpha.rank : old.beta.rank)
    width = length(rows); parent = old_rank + width
    output = @view capacity.inherited[1:parent, 1:parent]; fill!(output, 0.0)
    old === nothing && return output
    spin = alpha ? old.alpha : old.beta; offset = left ? 0 : width
    old_rank > 0 && copyto!(@view(output[offset+1:offset+old_rank,
        offset+1:offset+old_rank]), spin.core_field)
    @inbounds for site in eachindex(rows)
        value = 0.0; phase = _phase(operator, rows[site])
        for channel in eachindex(old.hartree_channel)
            value = muladd(left ? operator.phase_sink[channel, phase] :
                operator.phase_source[phase, channel],
                old.hartree_channel[channel], value)
        end
        index = _entry_site_position(site, old_rank, left)
        output[index, index] += value
    end
    output
end

function _transferred_hartree!(output, old, operator, left::Bool, workspace)
    fill!(output, 0.0); old === nothing && return output
    T = _channel_transfer!(workspace.transfer, workspace.transfer_scratch,
        operator, 1)
    left ? mul!(output, transpose(T), old.hartree_channel) :
        mul!(output, T, old.hartree_channel)
    output
end

function _entry_parent_cross!(workspace, old, alpha_parent, beta_parent,
        rows, operator, left::Bool)
    ma, mb = alpha_parent.parent, beta_parent.parent
    A = @view workspace.parent_alpha[1:pair_dimension(ma), :]
    B = @view workspace.parent_beta[1:pair_dimension(mb), :]
    fill!(A, 0.0); fill!(B, 0.0); column = 0
    olda, oldb, width = alpha_parent.old_rank, beta_parent.old_rank,
        length(rows)
    if old !== nothing
        record = old.cross_hartree
        @inbounds for factor = 1:cross_hartree_rank(record)
            column += 1
            for j = 1:olda, i = 1:j
                A[pair_index(_entry_old_position(i, width, left),
                    _entry_old_position(j, width, left)), column] =
                    record.alpha[pair_index(i, j), factor]
            end
            for j = 1:oldb, i = 1:j
                B[pair_index(_entry_old_position(i, width, left),
                    _entry_old_position(j, width, left)), column] =
                    record.beta[pair_index(i, j), factor]
            end
            workspace.weights[column] = record.weight[factor]
        end
    end
    channels = _maximum_channel_rank(operator)
    @inbounds for site = 1:width
        phase = _phase(operator, rows[site])
        qa = pair_index(_entry_site_position(site, olda, left),
            _entry_site_position(site, olda, left))
        qb = pair_index(_entry_site_position(site, oldb, left),
            _entry_site_position(site, oldb, left))
        column += 1
        for pair in axes(A, 1), channel = 1:channels
            A[pair, column] = muladd(alpha_parent.transformed[pair, channel],
                left ? operator.phase_sink[channel, phase] :
                operator.phase_source[phase, channel], A[pair, column])
        end
        B[qb, column] = 1.0; workspace.weights[column] = 1.0
        column += 1; A[qa, column] = 1.0
        for pair in axes(B, 1), channel = 1:channels
            B[pair, column] = muladd(beta_parent.transformed[pair, channel],
                left ? operator.phase_sink[channel, phase] :
                operator.phase_source[phase, channel], B[pair, column])
        end
        workspace.weights[column] = 1.0
    end
    @inbounds for other = 1:width
        column += 1
        qb = pair_index(_entry_site_position(other, oldb, left),
            _entry_site_position(other, oldb, left))
        B[qb, column] = 1.0
        for site = 1:width
            qa = pair_index(_entry_site_position(site, olda, left),
                _entry_site_position(site, olda, left))
            A[qa, column] = workspace.local_interaction[site, other]
        end
        workspace.weights[column] = 1.0
    end
    column <= workspace.maximum_columns || throw(DimensionMismatch(
        "cross-Hartree entry capacity"))
    UHFCrossHartree(@view(A[:, 1:column]), @view(B[:, 1:column]),
        @view(workspace.weights[1:column]), ma, mb)
end

function _entry_project_spin!(capacity, parent, link, old, rows, one_body,
        left::Bool, effective_core, provenance::UInt)
    m, active = parent.parent, link.active_rank
    K, Ka = pair_dimension(m), pair_dimension(active)
    T = _pair_transform!(capacity.pair_map, link.close_map, m, active)
    pair_scratch = @view capacity.pair_scratch[1:K, 1:Ka]
    mul!(pair_scratch, parent.internal, T)
    internal = Matrix{Float64}(undef, Ka, Ka)
    mul!(internal, transpose(T), pair_scratch)
    connector = Matrix{Float64}(undef, Ka, size(parent.connector, 2))
    mul!(connector, transpose(T), parent.connector)
    parent_h1 = _entry_parent_h1!(capacity, old, rows, one_body, left)
    dense = @view capacity.dense_scratch[1:m, 1:active]
    mul!(dense, parent_h1, link.close_map)
    h1 = Matrix{Float64}(undef, active, active)
    mul!(h1, transpose(link.close_map), dense)
    mul!(dense, effective_core, link.close_map)
    core_field = Matrix{Float64}(undef, active, active)
    mul!(core_field, transpose(link.close_map), dense)
    width = length(rows); old_rank = parent.old_rank
    rows_out = old === nothing || isempty(old.rows) ? rows :
        left ? (first(old.rows):last(rows)) : (first(rows):last(old.rows))
    cells_out = (((first(rows_out)-1) ÷ UNIT_CELL_WIDTH + 1):
        ((last(rows_out)-1) ÷ UNIT_CELL_WIDTH + 1))
    edge_width = min(length(rows_out), max(h1_bandwidth(one_body), 1))
    left_edge = zeros(edge_width, active); right_edge = zeros(edge_width, active)
    @inbounds for edge = 1:edge_width, column = 1:active, index = 1:m
        lsite = first(rows_out) + edge - 1
        rsite = last(rows_out) - edge_width + edge
        lvalue = if lsite in rows
            index == _entry_site_position(lsite-first(rows)+1, old_rank, left)
        elseif old !== nothing
            oi = index - (left ? 0 : width)
            1 <= oi <= old_rank ? _uhf_entry_edge(old, lsite, oi, true) : 0.0
        else
            0.0
        end
        rvalue = if rsite in rows
            index == _entry_site_position(rsite-first(rows)+1, old_rank, left)
        elseif old !== nothing
            oi = index - (left ? 0 : width)
            1 <= oi <= old_rank ? _uhf_entry_edge(old, rsite, oi, false) : 0.0
        else
            0.0
        end
        left_edge[edge, column] = muladd(lvalue, link.close_map[index, column],
            left_edge[edge, column])
        right_edge[edge, column] = muladd(rvalue, link.close_map[index, column],
            right_edge[edge, column])
    end
    completed_energy = dot(parent_h1,
        @view(capacity.completed[1:m, 1:m]))
    UHFEntrySpinBlock(cells_out, rows_out, active,
        (old === nothing ? 0 : old.completed_count) + link.selection.completed,
        (old === nothing ? 0 : old.empty_count) + link.selection.empty,
        (old === nothing ? 0.0 : old.one_body_energy) + completed_energy,
        h1, internal, connector, core_field, left_edge, right_edge, link,
        !left, provenance)
end

function _build_independent_spin_entry!(workspace,
        old::Union{Nothing,UHFEntryBlock}, alpha_link::RHFStateLink,
        beta_link::RHFStateLink, rows, one_body::BandedOneBody,
        operator::UnitCellInteraction, left::Bool, provenance::UInt)
    width = length(rows)
    width <= workspace.maximum_width || throw(DimensionMismatch(
        "independent-spin entry width capacity"))
    _entry_local_interaction!(workspace, operator, rows)
    alpha_parent = _bare_spin_parent!(workspace.alpha,
        old === nothing ? nothing : old.alpha, rows, operator, left, workspace)
    beta_parent = _bare_spin_parent!(workspace.beta,
        old === nothing ? nothing : old.beta, rows, operator, left, workspace)
    parent_cross = _entry_parent_cross!(workspace, old, alpha_parent,
        beta_parent, rows, operator, left)
    ma, mb = alpha_parent.parent, beta_parent.parent
    Ca = @view workspace.alpha.completed[1:ma, 1:ma]
    Cb = @view workspace.beta.completed[1:mb, 1:mb]
    mul!(Ca, alpha_link.completed_map, transpose(alpha_link.completed_map))
    mul!(Cb, beta_link.completed_map, transpose(beta_link.completed_map))
    Ja = @view workspace.alpha.coulomb[1:ma, 1:ma]
    Ka = @view workspace.alpha.exchange[1:ma, 1:ma]
    Jb = @view workspace.beta.coulomb[1:mb, 1:mb]
    Kb = @view workspace.beta.exchange[1:mb, 1:mb]
    _self_spin_fields!(Ja, Ka, Ca, alpha_parent.internal, ma)
    _self_spin_fields!(Jb, Kb, Cb, beta_parent.internal, mb)
    pa = _pack_density!(workspace.alpha.packed, Ca, ma)
    pb = _pack_density!(workspace.beta.packed, Cb, mb)
    fa = @view workspace.alpha.cross_field[1:pair_dimension(ma)]
    fb = @view workspace.beta.cross_field[1:pair_dimension(mb)]
    cross_energy = _cross_packed_fields!(fa, fb, pa, pb, parent_cross,
        workspace.factor_alpha, workspace.factor_beta)
    inherited_a = _fill_parent_inherited!(workspace.alpha, old, rows,
        operator, left, true)
    inherited_b = _fill_parent_inherited!(workspace.beta, old, rows,
        operator, left, false)
    self_energy = 0.0
    @inbounds for j = 1:ma, i = 1:ma
        self_energy += 0.5(Ja[i, j] - Ka[i, j]) * Ca[i, j]
    end
    @inbounds for j = 1:mb, i = 1:mb
        self_energy += 0.5(Jb[i, j] - Kb[i, j]) * Cb[i, j]
    end
    core_energy = (old === nothing ? 0.0 : old.core_energy) +
        dot(inherited_a, Ca) + dot(inherited_b, Cb) + self_energy +
        cross_energy
    @inbounds for j = 1:ma, i = 1:ma
        inherited_a[i, j] += Ja[i, j] - Ka[i, j] +
            fa[pair_index(i, j)] / pair_scale(i, j)
    end
    @inbounds for j = 1:mb, i = 1:mb
        inherited_b[i, j] += Jb[i, j] - Kb[i, j] +
            fb[pair_index(i, j)] / pair_scale(i, j)
    end
    transferred = _transferred_hartree!(workspace.transferred_hartree, old,
        operator, left, workspace)
    hartree = Vector{Float64}(undef, length(transferred))
    mul!(hartree, transpose(alpha_parent.connector), pa)
    mul!(hartree, transpose(beta_parent.connector), pb, 1.0, 1.0)
    hartree .+= transferred
    alpha = _entry_project_spin!(workspace.alpha, alpha_parent, alpha_link,
        old === nothing ? nothing : old.alpha, rows, one_body, left,
        inherited_a, provenance)
    beta = _entry_project_spin!(workspace.beta, beta_parent, beta_link,
        old === nothing ? nothing : old.beta, rows, one_body, left,
        inherited_b, provenance)
    Ta = _pair_transform!(workspace.alpha.pair_map, alpha_link.close_map,
        ma, alpha.rank)
    Tb = _pair_transform!(workspace.beta.pair_map, beta_link.close_map,
        mb, beta.rank)
    columns = cross_hartree_rank(parent_cross)
    Aa = @view workspace.transformed_alpha[1:pair_dimension(alpha.rank),
        1:columns]
    Bb = @view workspace.transformed_beta[1:pair_dimension(beta.rank),
        1:columns]
    mul!(Aa, transpose(Ta), parent_cross.alpha)
    mul!(Bb, transpose(Tb), parent_cross.beta)
    cross = _compress_cross_hartree(Matrix(Aa), Matrix(Bb),
        Vector(parent_cross.weight), alpha.rank, beta.rank)
    block = UHFEntryBlock(alpha, beta, cross, hartree, core_energy, provenance)
    _validate_uhf_entry(block, operator)
    block
end

function reflect_uhf_entry(block::UHFEntryBlock, total_cells::Int,
        total_rows::Int, operator::UnitCellInteraction)
    function reflect_spin(spin)
        rows = isempty(spin.rows) ? (1:0) :
            (total_rows-last(spin.rows)+1):(total_rows-first(spin.rows)+1)
        cells = isempty(spin.cells) ? (1:0) :
            (total_cells-last(spin.cells)+1):(total_cells-first(spin.cells)+1)
        UHFEntrySpinBlock(cells, rows, spin.rank, spin.completed_count,
            spin.empty_count, spin.one_body_energy, copy(spin.h1),
            copy(spin.internal), copy(spin.connector), copy(spin.core_field),
            copy(spin.right_edge), copy(spin.left_edge), spin.link,
            !spin.reflected, spin.provenance)
    end
    output = UHFEntryBlock(reflect_spin(block.alpha), reflect_spin(block.beta),
        UHFCrossHartree(copy(block.cross_hartree.alpha),
            copy(block.cross_hartree.beta), copy(block.cross_hartree.weight),
            block.alpha.rank, block.beta.rank), copy(block.hartree_channel),
        block.core_energy, block.provenance)
    _validate_uhf_entry(output, operator); output
end

function uhf_entry_storage_bytes(block::UHFEntryBlock)
    spin_bytes(spin) = sum(sizeof, (spin.h1, spin.internal, spin.connector,
        spin.core_field, spin.left_edge, spin.right_edge, spin.link.close_map,
        spin.link.completed_map))
    spin_bytes(block.alpha) + spin_bytes(block.beta) + sum(sizeof,
        (block.cross_hartree.alpha, block.cross_hartree.beta,
        block.cross_hartree.weight, block.hartree_channel))
end
