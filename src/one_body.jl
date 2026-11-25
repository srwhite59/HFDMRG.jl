struct BandedOneBody
    bands::Matrix{Float64}

    function BandedOneBody(bands::AbstractMatrix{Float64})
        rows, sites = size(bands)
        1 <= rows <= sites || throw(DimensionMismatch(
            "banded H1 dimensions are invalid"))
        owned = Matrix{Float64}(bands)
        all(isfinite, owned) || throw(ArgumentError("H1 bands must be finite"))
        @inbounds for offset = 1:rows-1, site = sites-offset+1:sites
            iszero(owned[offset + 1, site]) || throw(ArgumentError(
                "inactive H1 band storage must be zero"))
        end
        new(owned)
    end
end

@inline h1_sites(one_body::BandedOneBody) = size(one_body.bands, 2)
@inline h1_bandwidth(one_body::BandedOneBody) = size(one_body.bands, 1) - 1
@inline function h1_entry(one_body::BandedOneBody, first::Int, second::Int)
    left, right = minmax(first, second)
    offset = right - left
    offset > h1_bandwidth(one_body) && return 0.0
    @inbounds one_body.bands[offset + 1, left]
end

function _fill_h1_outer_coupling!(h1, left, right, right_offset::Int,
        one_body::BandedOneBody)
    (left.rank == 0 || right.rank == 0) && return nothing
    left_last, right_first = last(left.rows), first(right.rows)
    reach = h1_bandwidth(one_body) - (right_first - left_last) + 1
    reach <= 0 && return nothing

    # Project the band between the two stored edge bases. Only pairs with
    # physical separation at most the bandwidth contribute.
    left_width = size(left.right_edge, 1)
    right_width = size(right.left_edge, 1)
    left_offset = left_last - left_width
    first_edge = max(1, left_width - reach + 1)
    @inbounds for a = 1:left.rank, b = 1:right.rank
        coupling = 0.0
        for i = first_edge:left_width, j = 1:min(right_width, reach - left_width + i)
            coupling += left.right_edge[i, a] *
                h1_entry(one_body, left_offset + i, right_first + j - 1) *
                right.left_edge[j, b]
        end
        h1[a, right_offset + b] = coupling
        h1[right_offset + b, a] = coupling
    end
    nothing
end
