"""
Slice layout for grouping orbitals into contiguous slices.

dims lists the number of orbitals per slice. offs is the prefix sum with
offs[1] == 0 and offs[i+1] == sum(dims[1:i]).
"""
struct SliceLayout
    dims::Vector{Int}
    offs::Vector{Int}
    function SliceLayout(dims::Vector{Int})
        isempty(dims) && error("dims must be non-empty")
        any(d -> d <= 0, dims) && error("dims must be positive")
        offs = Vector{Int}(undef, length(dims) + 1)
        offs[1] = 0
        for i = 1:length(dims)
            offs[i + 1] = offs[i] + dims[i]
        end
        new(dims, offs)
    end
end

"""
Number of slices in the layout.
"""
nslices(layout::SliceLayout) = length(layout.dims)

"""
Size of slice n in orbitals.
"""
function slice_dims(layout::SliceLayout, n::Integer)
    1 <= n <= length(layout.dims) || error("slice index out of bounds")
    layout.dims[n]
end

"""
Orbital range for a single slice n.
"""
function orb_range(layout::SliceLayout, n::Int)
    1 <= n <= length(layout.dims) || error("slice index out of bounds")
    (layout.offs[n] + 1):layout.offs[n + 1]
end

"""
Orbital range for a contiguous slice range n1:n2.
"""
function orb_range(layout::SliceLayout, ns::UnitRange{Int})
    n1, n2 = first(ns), last(ns)
    1 <= n1 <= n2 <= length(layout.dims) || error("slice range out of bounds")
    (layout.offs[n1] + 1):layout.offs[n2 + 1]
end

"""
Map a flattened orbital index p to (slice, local index) in the layout.
"""
function slice_local(layout::SliceLayout, p::Int)
    1 <= p <= layout.offs[end] || error("orbital index out of bounds")
    n = searchsortedlast(layout.offs, p - 1)
    a = p - layout.offs[n]
    n, a
end
