"""
Reusable direct two-sided center for independent alpha/beta spaces.

Algorithmic donor: PPP 0bcd4f11831a98ae4bf4fb3d09eac5968df9d96c,
`_IndependentSpinTwoSidedCenters.jl`: `_IndependentSpinSideView`,
`_PreparedIndependentOneSpinTwoSided`, `_IndependentSpinTwoSidedWorkspace`,
`_prepare_one_spin_two_sided!`, `_self_spin_fock!`,
`_internal_cross_field!`, `_independent_cross_fields!`,
`_prepare_independent_spin_two_sided!`, and
`_independent_spin_two_sided_fock_core!`.

Adapted to HFDMRG's `UHFOuterBlock`, two-group moving center, banded-H1
owner, and result contract.  Opposite-spin Hartree fields are contracted
directly from block-local factors, MPO channels, and physical charges.  No
cap-sized packed-pair cross panel is materialized.
"""
mutable struct FastUHFPreparedCenter
    alpha::FastRHFPreparedCenter
    beta::FastRHFPreparedCenter
    cross_alpha::Matrix{Float64}
    cross_beta::Matrix{Float64}
    alpha_left_density::Vector{Float64}
    alpha_right_density::Vector{Float64}
    beta_left_density::Vector{Float64}
    beta_right_density::Vector{Float64}
    alpha_left_field::Vector{Float64}
    alpha_right_field::Vector{Float64}
    beta_left_field::Vector{Float64}
    beta_right_field::Vector{Float64}
    alpha_factor::Vector{Float64}
    beta_factor::Vector{Float64}
    channel::Vector{Float64}
    left::Union{Nothing,UHFOuterBlock}
    right::Union{Nothing,UHFOuterBlock}
    constant_energy::Float64
    preparation_calls::Int
    fock_calls::Int
end

function FastUHFPreparedCenter(maximum_alpha::Int, maximum_beta::Int,
        operator::UnitCellInteraction; maximum_rhs::Int=3,
        maximum_alpha_side::Int=max(maximum_alpha - 2UNIT_CELL_WIDTH, 1),
        maximum_beta_side::Int=max(maximum_beta - 2UNIT_CELL_WIDTH, 1))
    min(maximum_alpha, maximum_beta) >= 1 || throw(ArgumentError(
        "UHF center capacities must be positive"))
    maximum_side = max(maximum_alpha_side, maximum_beta_side)
    pairs = pair_dimension(maximum_side)
    vector() = zeros(pairs)
    FastUHFPreparedCenter(
        FastRHFPreparedCenter(maximum_alpha, operator; maximum_rhs,
            maximum_side_rank=maximum_alpha_side),
        FastRHFPreparedCenter(maximum_beta, operator; maximum_rhs,
            maximum_side_rank=maximum_beta_side),
        zeros(maximum_alpha, maximum_alpha),
        zeros(maximum_beta, maximum_beta),
        vector(), vector(), vector(), vector(), vector(), vector(), vector(),
        vector(), vector(), vector(), zeros(_maximum_channel_rank(operator)),
        nothing, nothing, 0.0, 0, 0)
end

function prepare_uhf_lifecycle_center!(prepared::FastUHFPreparedCenter,
        left::UHFOuterBlock, right::UHFOuterBlock, cell::Int,
        first_rows::UnitRange{Int}, second_rows::UnitRange{Int},
        one_body::BandedOneBody, operator::UnitCellInteraction)
    _validate_uhf_block(left, operator); _validate_uhf_block(right, operator)
    left.provenance == right.provenance || throw(ArgumentError(
        "UHF center provenance is inconsistent"))
    _prepare_lifecycle_center!(prepared.alpha, left.alpha, right.alpha, cell,
        first_rows, second_rows, one_body, operator)
    _prepare_lifecycle_center!(prepared.beta, left.beta, right.beta, cell,
        first_rows, second_rows, one_body, operator)
    envelope = 2048eps(Float64) * max(abs(prepared.alpha.constant_energy), 1.0)
    abs(prepared.alpha.constant_energy - prepared.beta.constant_energy) <=
        envelope || throw(ArgumentError("spin center core scalars disagree"))
    prepared.left = left; prepared.right = right
    prepared.constant_energy = prepared.alpha.constant_energy
    prepared.preparation_calls += 1
    prepared
end

function _fast_uhf_preflight(prepared::FastUHFPreparedCenter,
        alpha_covariance, beta_covariance, alpha_occupied::Int,
        beta_occupied::Int)
    ma, mb = prepared.alpha.active_rank, prepared.beta.active_rank
    size(alpha_covariance) == (ma, ma) && size(beta_covariance) == (mb, mb) ||
        throw(DimensionMismatch("UHF center covariance dimensions disagree"))
    0 <= alpha_occupied <= ma && 0 <= beta_occupied <= mb ||
        throw(DimensionMismatch("UHF center occupation exceeds a spin rank"))
    nothing
end

function _fast_one_spin_fields!(spin::FastRHFPreparedCenter, covariance)
    rank = spin.active_rank
    left = something(spin.left); right = something(spin.right)
    operator = spin.operator
    rows = first(spin.first_rows):last(spin.second_rows)
    ml, width, mr = left.rank, length(rows), right.rank
    density = @view spin.density[1:rank, 1:rank]
    copyto!(density, covariance)
    J = @view spin.coulomb[1:rank, 1:rank]
    K = @view spin.exchange[1:rank, 1:rank]
    fill!(J, 0.0); fill!(K, 0.0)
    ml > 0 && _self_spin_fields!(@view(J[1:ml, 1:ml]),
        @view(K[1:ml, 1:ml]), @view(density[1:ml, 1:ml]),
        left.pair_field, ml)
    offset_right = ml + width
    mr > 0 && _self_spin_fields!(@view(J[offset_right+1:rank,
        offset_right+1:rank]), @view(K[offset_right+1:rank,
        offset_right+1:rank]), @view(density[offset_right+1:rank,
        offset_right+1:rank]), right.pair_field, mr)
    residual = size(operator.residual_transfer, 1)
    scratch1 = @view spin.interaction_scratch_first[1:residual]
    scratch2 = @view spin.interaction_scratch_second[1:residual]
    @inbounds for b = 1:width, a = 1:width
        value = _interaction_entry(operator, rows[a], rows[b], scratch1,
            scratch2)
        K[ml+a, ml+b] = value * density[ml+a, ml+b]
        J[ml+a, ml+a] = muladd(value, density[ml+b, ml+b], J[ml+a, ml+a])
    end
    channels = _maximum_channel_rank(operator)
    transfer0 = @view spin.transfer_zero[1:channels, 1:channels]
    transfer1 = @view spin.transfer_one[1:channels, 1:channels]
    transfer2 = @view spin.transfer_two[1:channels, 1:channels]
    first_width = length(spin.first_rows)
    ml > 0 && _fast_block_physical_hartree!(J, left, 0, ml,
        spin.first_rows, @view(density[1:ml, 1:ml]),
        @view(density[ml+1:ml+first_width, ml+1:ml+first_width]), operator,
        transfer0, spin, true)
    ml > 0 && _fast_block_physical_exchange!(K, density, left, 0, ml, 0,
        spin.first_rows, spin.left_center_exchange)
    ml > 0 && _fast_block_physical_hartree!(J, left, 0, ml+first_width,
        spin.second_rows, @view(density[1:ml, 1:ml]),
        @view(density[ml+first_width+1:ml+width,
            ml+first_width+1:ml+width]), operator, transfer1, spin, true)
    ml > 0 && _fast_block_physical_exchange!(K, density, left, 0,
        ml+first_width, first_width, spin.second_rows,
        spin.left_center_exchange)
    mr > 0 && _fast_block_physical_hartree!(J, right, offset_right,
        ml+first_width, spin.second_rows,
        @view(density[offset_right+1:rank, offset_right+1:rank]),
        @view(density[ml+first_width+1:ml+width,
            ml+first_width+1:ml+width]), operator, transfer0, spin, false)
    mr > 0 && _fast_block_physical_exchange!(K, density, right,
        offset_right, ml+first_width, first_width, spin.second_rows,
        spin.right_center_exchange)
    mr > 0 && _fast_block_physical_hartree!(J, right, offset_right, ml,
        spin.first_rows,
        @view(density[offset_right+1:rank, offset_right+1:rank]),
        @view(density[ml+1:ml+first_width, ml+1:ml+first_width]), operator,
        transfer1, spin, false)
    mr > 0 && _fast_block_physical_exchange!(K, density, right,
        offset_right, ml, 0, spin.first_rows, spin.right_center_exchange)
    ml > 0 && mr > 0 && _fast_block_block_hartree!(J, left, right,
        offset_right, @view(density[1:ml, 1:ml]),
        @view(density[offset_right+1:rank, offset_right+1:rank]), transfer2,
        spin)
    ml > 0 && mr > 0 && _fast_block_block_exchange!(K, density, left, right,
        offset_right, spin.left_right_exchange)
    density, J, K
end

function _add_packed_cross_block!(output, packed, offset::Int, rank::Int)
    @inbounds for j = 1:rank, i = 1:j
        value = packed[pair_index(i, j)] / pair_scale(i, j)
        output[offset+i, offset+j] += value
        i != j && (output[offset+j, offset+i] += value)
    end
    nothing
end

function _block_cross_fields!(cross_alpha, cross_beta, alpha_density,
        beta_density, block::UHFOuterBlock, alpha_offset::Int,
        beta_offset::Int, prepared::FastUHFPreparedCenter, side::Symbol)
    ma, mb = block.alpha.rank, block.beta.rank
    ma == 0 && mb == 0 && return nothing
    pa_store = side === :left ? prepared.alpha_left_density :
        prepared.alpha_right_density
    pb_store = side === :left ? prepared.beta_left_density :
        prepared.beta_right_density
    fa_store = side === :left ? prepared.alpha_left_field :
        prepared.alpha_right_field
    fb_store = side === :left ? prepared.beta_left_field :
        prepared.beta_right_field
    pa = _pack_density!(pa_store,
        @view(alpha_density[alpha_offset+1:alpha_offset+ma,
            alpha_offset+1:alpha_offset+ma]), ma)
    pb = _pack_density!(pb_store,
        @view(beta_density[beta_offset+1:beta_offset+mb,
            beta_offset+1:beta_offset+mb]), mb)
    fa = @view fa_store[1:pair_dimension(ma)]
    fb = @view fb_store[1:pair_dimension(mb)]
    _cross_packed_fields!(fa, fb, pa, pb, block.cross_hartree,
        prepared.alpha_factor, prepared.beta_factor)
    nothing
end

function _fast_uhf_cross_fields!(prepared::FastUHFPreparedCenter,
        left::UHFOuterBlock, right::UHFOuterBlock)
    alpha, beta = prepared.alpha, prepared.beta
    Da = @view alpha.density[1:alpha.active_rank, 1:alpha.active_rank]
    Db = @view beta.density[1:beta.active_rank, 1:beta.active_rank]
    Ca = @view prepared.cross_alpha[1:alpha.active_rank, 1:alpha.active_rank]
    Cb = @view prepared.cross_beta[1:beta.active_rank, 1:beta.active_rank]
    fill!(Ca, 0.0); fill!(Cb, 0.0)
    la, lb = left.alpha.rank, left.beta.rank
    ra, rb = right.alpha.rank, right.beta.rank
    width = length(alpha.first_rows) + length(alpha.second_rows)
    oa, ob = la + width, lb + width
    _block_cross_fields!(Ca, Cb, Da, Db, left, 0, 0, prepared, :left)
    _block_cross_fields!(Ca, Cb, Da, Db, right, oa, ob, prepared, :right)

    pa_l = @view prepared.alpha_left_density[1:pair_dimension(la)]
    pa_r = @view prepared.alpha_right_density[1:pair_dimension(ra)]
    pb_l = @view prepared.beta_left_density[1:pair_dimension(lb)]
    pb_r = @view prepared.beta_right_density[1:pair_dimension(rb)]
    fa_l = @view prepared.alpha_left_field[1:pair_dimension(la)]
    fa_r = @view prepared.alpha_right_field[1:pair_dimension(ra)]
    fb_l = @view prepared.beta_left_field[1:pair_dimension(lb)]
    fb_r = @view prepared.beta_right_field[1:pair_dimension(rb)]
    channels = length(prepared.channel)
    channel = prepared.channel
    if la > 0 && rb > 0
        propagated = @view alpha.propagated[1:pair_dimension(la), 1:channels]
        mul!(channel, transpose(right.beta.left_feature), pb_r)
        mul!(fa_l, propagated, channel, 1.0, 1.0)
        mul!(channel, transpose(propagated), pa_l)
        mul!(fb_r, right.beta.left_feature, channel, 1.0, 1.0)
    end
    if lb > 0 && ra > 0
        propagated = @view beta.propagated[1:pair_dimension(lb), 1:channels]
        mul!(channel, transpose(propagated), pb_l)
        mul!(fa_r, right.alpha.left_feature, channel, 1.0, 1.0)
        mul!(channel, transpose(right.alpha.left_feature), pa_r)
        mul!(fb_l, propagated, channel, 1.0, 1.0)
    end
    la > 0 && _add_packed_cross_block!(Ca, fa_l, 0, la)
    lb > 0 && _add_packed_cross_block!(Cb, fb_l, 0, lb)
    ra > 0 && _add_packed_cross_block!(Ca, fa_r, oa, ra)
    rb > 0 && _add_packed_cross_block!(Cb, fb_r, ob, rb)

    @inbounds for site = 1:width
        alpha_charge = Da[la+site, la+site]
        beta_charge = Db[lb+site, lb+site]
        for c = 1:la, a = 1:la
            Ca[a, c] = muladd(alpha.left_center_exchange[a, c, site],
                beta_charge, Ca[a, c])
            Cb[lb+site, lb+site] = muladd(Da[a, c],
                alpha.left_center_exchange[a, c, site],
                Cb[lb+site, lb+site])
        end
        for c = 1:lb, a = 1:lb
            Cb[a, c] = muladd(beta.left_center_exchange[a, c, site],
                alpha_charge, Cb[a, c])
            Ca[la+site, la+site] = muladd(Db[a, c],
                beta.left_center_exchange[a, c, site],
                Ca[la+site, la+site])
        end
        for c = 1:ra, a = 1:ra
            Ca[oa+a, oa+c] = muladd(alpha.right_center_exchange[a, c, site],
                beta_charge, Ca[oa+a, oa+c])
            Cb[lb+site, lb+site] = muladd(Da[oa+a, oa+c],
                alpha.right_center_exchange[a, c, site],
                Cb[lb+site, lb+site])
        end
        for c = 1:rb, a = 1:rb
            Cb[ob+a, ob+c] = muladd(beta.right_center_exchange[a, c, site],
                alpha_charge, Cb[ob+a, ob+c])
            Ca[la+site, la+site] = muladd(Db[ob+a, ob+c],
                beta.right_center_exchange[a, c, site],
                Ca[la+site, la+site])
        end
    end
    operator = alpha.operator
    rows = first(alpha.first_rows):last(alpha.second_rows)
    residual = size(operator.residual_transfer, 1)
    scratch1 = @view alpha.interaction_scratch_first[1:residual]
    scratch2 = @view alpha.interaction_scratch_second[1:residual]
    @inbounds for site = 1:width, other = 1:width
        value = _interaction_entry(operator, rows[site], rows[other],
            scratch1, scratch2)
        Ca[la+site, la+site] = muladd(value, Db[lb+other, lb+other],
            Ca[la+site, la+site])
        Cb[lb+site, lb+site] = muladd(value, Da[la+other, la+other],
            Cb[lb+site, lb+site])
    end
    nothing
end

function _fast_uhf_energy!(prepared::FastUHFPreparedCenter,
        alpha_covariance, beta_covariance, alpha_occupied::Int,
        beta_occupied::Int;
        construct_fock::Bool)
    _fast_uhf_preflight(prepared, alpha_covariance, beta_covariance,
        alpha_occupied, beta_occupied)
    Da, Ja, Ka = _fast_one_spin_fields!(prepared.alpha, alpha_covariance)
    Db, Jb, Kb = _fast_one_spin_fields!(prepared.beta, beta_covariance)
    _fast_uhf_cross_fields!(prepared, something(prepared.left),
        something(prepared.right))
    ma, mb = prepared.alpha.active_rank, prepared.beta.active_rank
    Ca = @view prepared.cross_alpha[1:ma, 1:ma]
    Cb = @view prepared.cross_beta[1:mb, 1:mb]
    h1a = @view prepared.alpha.h1[1:ma, 1:ma]
    h1b = @view prepared.beta.h1[1:mb, 1:mb]
    if construct_fock
        Fa = @view prepared.alpha.fock[1:ma, 1:ma]
        Fb = @view prepared.beta.fock[1:mb, 1:mb]
        @. Fa = h1a + Ja - Ka + Ca
        @. Fb = h1b + Jb - Kb + Cb
        prepared.fock_calls += 1
    end
    one = dot(h1a, Da) + dot(h1b, Db)
    self_hartree = 0.5dot(Ja, Da) + 0.5dot(Jb, Db)
    exchange = -0.5dot(Ka, Da) - 0.5dot(Kb, Db)
    cross_a = dot(Ca, Da); cross_b = dot(Cb, Db)
    envelope = 4096eps(Float64) * max(abs(cross_a), abs(cross_b), 1.0)
    abs(cross_a - cross_b) <= envelope || throw(ArgumentError(
        "opposite-spin Hartree derivatives are inconsistent"))
    hartree = self_hartree + cross_a
    (total=prepared.constant_energy + one + hartree + exchange,
        core=prepared.constant_energy, one_body=one, hartree, exchange,
        cross_hartree=cross_a)
end

function uhf_center_energy_fock!(prepared::FastUHFPreparedCenter,
        alpha_covariance, beta_covariance, alpha_occupied::Int,
        beta_occupied::Int)
    _fast_uhf_energy!(prepared, alpha_covariance, beta_covariance,
        alpha_occupied, beta_occupied; construct_fock=true)
end

function uhf_center_true_energy!(prepared::FastUHFPreparedCenter,
        alpha_covariance, beta_covariance, alpha_occupied::Int,
        beta_occupied::Int)
    _fast_uhf_energy!(prepared, alpha_covariance, beta_covariance,
        alpha_occupied, beta_occupied; construct_fock=false)
end

function apply_uhf_center_fock!(alpha_output::AbstractVecOrMat{Float64},
        beta_output::AbstractVecOrMat{Float64},
        prepared::FastUHFPreparedCenter, alpha_input::AbstractVecOrMat{Float64},
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

function fast_uhf_center_storage_bytes(prepared::FastUHFPreparedCenter)
    fast_rhf_center_storage_bytes(prepared.alpha) +
        fast_rhf_center_storage_bytes(prepared.beta) + sum(sizeof, (
        prepared.cross_alpha, prepared.cross_beta,
        prepared.alpha_left_density, prepared.alpha_right_density,
        prepared.beta_left_density, prepared.beta_right_density,
        prepared.alpha_left_field, prepared.alpha_right_field,
        prepared.beta_left_field, prepared.beta_right_field,
        prepared.alpha_factor, prepared.beta_factor, prepared.channel))
end
