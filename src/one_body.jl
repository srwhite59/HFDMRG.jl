struct BandedOneBody
    bands::Matrix{Float64}
end

function BandedOneBody(bands::AbstractMatrix{Float64})
    rows, sites = size(bands)
    1 <= rows <= sites || throw(DimensionMismatch(
        "banded H1 dimensions are invalid"))
    all(isfinite, bands) || throw(ArgumentError("H1 bands must be finite"))
    owned = Matrix{Float64}(bands)
    @inbounds for offset = 1:rows-1, site = sites-offset+1:sites
        iszero(owned[offset + 1, site]) || throw(ArgumentError(
            "inactive H1 band storage must be zero"))
    end
    BandedOneBody(owned)
end

@inline h1_sites(one_body::BandedOneBody) = size(one_body.bands, 2)
@inline h1_bandwidth(one_body::BandedOneBody) = size(one_body.bands, 1) - 1
@inline function h1_entry(one_body::BandedOneBody, first::Int, second::Int)
    left, right = minmax(first, second)
    offset = right - left
    offset > h1_bandwidth(one_body) && return 0.0
    @inbounds one_body.bands[offset + 1, left]
end
