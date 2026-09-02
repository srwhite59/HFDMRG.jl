const _LOCAL_TRIANGLE_COUNT = UNIT_CELL_WIDTH * (UNIT_CELL_WIDTH + 1) ÷ 2

struct UnitCellInteraction
    sites::Int
    phase_at_first::Int
    tail_transfer::Vector{Float64}
    residual_transfer::Matrix{Float64}
    phase_source::Matrix{Float64}
    phase_sink::Matrix{Float64}
    local_triangle::Vector{Float64}
    provenance::UInt
end

function UnitCellInteraction(sites::Int, phase_at_first::Int,
        tail_transfer::AbstractVector{Float64},
        residual_transfer::AbstractMatrix{Float64},
        phase_source::AbstractMatrix{Float64},
        phase_sink::AbstractMatrix{Float64},
        local_triangle::AbstractVector{Float64})
    sites >= 2 || throw(ArgumentError("interaction length must be at least two"))
    1 <= phase_at_first <= UNIT_CELL_WIDTH || throw(ArgumentError(
        "first phase must lie in 1:18"))
    size(residual_transfer, 1) == size(residual_transfer, 2) ||
        throw(DimensionMismatch("residual transfer must be square"))
    tail_rank = length(tail_transfer)
    residual_rank = size(residual_transfer, 1)
    rank = Base.Checked.checked_add(tail_rank, residual_rank)
    size(phase_source) == (UNIT_CELL_WIDTH, rank) ||
        throw(DimensionMismatch("phase source dimensions are wrong"))
    size(phase_sink) == (rank, UNIT_CELL_WIDTH) ||
        throw(DimensionMismatch("phase sink dimensions are wrong"))
    length(local_triangle) == _LOCAL_TRIANGLE_COUNT ||
        throw(DimensionMismatch("local triangle must contain 171 values"))
    arrays = (tail_transfer, residual_transfer, phase_source, phase_sink,
        local_triangle)
    all(array -> all(isfinite, array), arrays) || throw(ArgumentError(
        "interaction coefficients must be finite"))
    owned = (Vector{Float64}(tail_transfer),
        Matrix{Float64}(residual_transfer), Matrix{Float64}(phase_source),
        Matrix{Float64}(phase_sink), Vector{Float64}(local_triangle))
    provenance = hash((sites, phase_at_first, owned...))
    UnitCellInteraction(sites, phase_at_first, owned..., provenance)
end

@inline _phase(operator::UnitCellInteraction, site::Int) =
    1 + mod(operator.phase_at_first + site - 2, UNIT_CELL_WIDTH)

@inline function _local_index(first_phase::Int, second_phase::Int)
    left, right = minmax(first_phase, second_phase)
    right * (right - 1) ÷ 2 + left
end

@inline function _first_group_count(operator::UnitCellInteraction)
    min(operator.sites, UNIT_CELL_WIDTH - operator.phase_at_first + 1)
end

@inline function _site_group_phase(operator::UnitCellInteraction, site::Int)
    1 <= site <= operator.sites || throw(BoundsError(1:operator.sites, site))
    first_count = _first_group_count(operator)
    site <= first_count && return 1, operator.phase_at_first + site - 1
    relative = site - first_count - 1
    2 + relative ÷ UNIT_CELL_WIDTH, 1 + relative % UNIT_CELL_WIDTH
end

function _interaction_entry(operator::UnitCellInteraction, first_site::Int,
        second_site::Int, first_scratch::AbstractVector{Float64},
        second_scratch::AbstractVector{Float64})
    left, right = minmax(first_site, second_site)
    left_group, left_phase = _site_group_phase(operator, left)
    right_group, right_phase = _site_group_phase(operator, right)
    left_group == right_group && return operator.local_triangle[
        _local_index(left_phase, right_phase)]
    residual_rank = size(operator.residual_transfer, 1)
    length(first_scratch) >= residual_rank &&
        length(second_scratch) >= residual_rank || throw(DimensionMismatch(
            "interaction scratch capacity is too small"))
    gap = right_group - left_group
    value = 0.0
    @inbounds for channel in eachindex(operator.tail_transfer)
        value = muladd(operator.phase_source[left_phase, channel] *
            operator.tail_transfer[channel]^(gap - 1),
            operator.phase_sink[channel, right_phase], value)
    end
    residual_rank == 0 && return value
    offset = length(operator.tail_transfer)
    old, new = first_scratch, second_scratch
    @inbounds for channel = 1:residual_rank
        old[channel] = operator.phase_source[left_phase, offset + channel]
    end
    @inbounds for _ = 1:gap-1
        for new_channel = 1:residual_rank
            accumulator = 0.0
            for old_channel = 1:residual_rank
                accumulator = muladd(operator.residual_transfer[
                    old_channel, new_channel], old[old_channel], accumulator)
            end
            new[new_channel] = accumulator
        end
        old, new = new, old
    end
    @inbounds for channel = 1:residual_rank
        value = muladd(old[channel],
            operator.phase_sink[offset + channel, right_phase], value)
    end
    value
end

function _maximum_channel_rank(operator::UnitCellInteraction)
    length(operator.tail_transfer) + size(operator.residual_transfer, 1)
end
