"""
Cached sliced-basis backend.

Uses retained per-slice and environment tensors for aligned local Fock builds;
windows that split a slice delegate to the projection backend.
"""
struct SlicedBasisBackendCached{T}
    layout::SliceLayout
    V::T
end

struct SlicedCachedBlockState
    ra::UnitRange{Int}
    phi::Matrix{Float64}
    G::Matrix{Float64}
    exterior_slices::UnitRange{Int}
    pair_offs::Vector{Int}
    W::Matrix{Float64}
    # Packed Wswap[:,pair(s,a,b)] uses the swapped orientation
    # V6[a,b,c,d,s,n] without assuming additional physical symmetries.
    Wswap::Matrix{Float64}
end

struct SlicedCachedCrossScratch
    density::Matrix{Float64}
    physical_work::Matrix{Float64}
    mixed::Matrix{Float64}
    exchange::Matrix{Float64}
    pair_field::Vector{Float64}
end

struct SlicedBasisCachedWindow{T}
    V::T
    layout::SliceLayout
    ml::Int
    lc::Int
    mr::Int
    center_slices::UnitRange{Int}
    L::SlicedCachedBlockState
    R::SlicedCachedBlockState
    scratch::SlicedCachedCrossScratch
end

const _bench_timing = Dict{Symbol, Float64}(
    :cache_init_calls => 0.0,
    :cache_incremental_calls => 0.0,
    :cache_fallback_calls => 0.0,
    :window_local_calls => 0.0,
    :window_projection_calls => 0.0,
)

_bench_timing_enabled() = !isempty(get(ENV, "HFDMRG_BENCH_TIMING", ""))

function _bench_timing_reset!()
    for key in keys(_bench_timing)
        _bench_timing[key] = 0.0
    end
    nothing
end

_bench_timing_snapshot() = Dict(key => value for (key, value) in _bench_timing)

_pair_offs(layout, slices) = [0; cumsum(abs2.(layout.dims[slices]))]
_pair_range(offs, s) = offs[s] + 1:offs[s + 1]
_state_pair_range(state, s) = _pair_range(state.pair_offs, s - first(state.exterior_slices) + 1)
_packed_tensor(W, state, s, m, d) = reshape(view(W, :, _state_pair_range(state, s)), m, m, d, d)

function SlicedBasisBackendCached(layout::SliceLayout, V6::Array{Float64,6})
    ns = nslices(layout)
    nj = layout.dims[1]
    any(d -> d != nj, layout.dims) && error("layout must have fixed slice size")
    size(V6) == (nj, nj, nj, nj, ns, ns) || error("V6 must have size (nj,nj,nj,nj,ns,ns)")
    SlicedBasisBackendCached{Array{Float64,6}}(layout, V6)
end

function SlicedBasisBackendCached(layout::SliceLayout, vee::SlicedVeeRagged)
    layout.dims == vee.layout.dims || error("layout does not match ragged interaction")
    SlicedBasisBackendCached{SlicedVeeRagged}(layout, vee)
end

SlicedBasisBackendCached(vee::SlicedVeeRagged) = SlicedBasisBackendCached(vee.layout, vee)

SlicedBasisBackendCached(layout::SliceLayout, Vblocks::Vector{Vector{Array{Float64,4}}}) =
    SlicedBasisBackendCached(SlicedVeeRagged(layout, Vblocks))

_slice_vee(V::Array{Float64,6}, n::Int, m::Int) = @view V[:, :, :, :, n, m]
_slice_vee(vee::SlicedVeeRagged, n::Int, m::Int) = vee.V[n][m]

_slice_local(layout::SliceLayout, p::Int) = slice_local(layout, p)

function _build_slice_basis(layout::SliceLayout, ra::UnitRange{Int}, phi::Matrix{Float64})
    ns = nslices(layout)
    m = size(phi, 2)
    P = [zeros(Float64, layout.dims[s], m) for s in 1:ns]
    for (row, p) in enumerate(ra)
        n, a = _slice_local(layout, p)
        P[n][a, :] = phi[row, :]
    end
    P
end

function _pair_map_mat!(U::AbstractMatrix{Float64}, P::AbstractMatrix{Float64})
    nj, m = size(P)
    size(U, 1) == nj * nj || error("pair map has wrong row count")
    size(U, 2) == m * m || error("pair map has wrong column count")
    for a = 1:nj, b = 1:nj
        row = a + (b - 1) * nj
        for k = 1:m, l = 1:m
            col = k + (l - 1) * m
            U[row, col] = P[a, k] * P[b, l]
        end
    end
    U
end

function _pair_map_mat(P::AbstractMatrix{Float64})
    nj, m = size(P)
    U = zeros(Float64, nj * nj, m * m)
    _pair_map_mat!(U, P)
end

function _active_slices(layout::SliceLayout, ra::UnitRange{Int})
    first(_slice_local(layout, first(ra))):first(_slice_local(layout, last(ra)))
end

_add_grouped!(G, L, M, R, tmp) = (work = view(tmp, :, 1:size(M, 2));
    mul!(work, L', M); mul!(G, work, R, 1.0, 1.0))
function _direct_grouped(P, active, m, backend)
    G, tmp = zeros(Float64, m^2, m^2), zeros(Float64, m^2, maximum(abs2.(backend.layout.dims[active])))
    for n in active
        Rn = _pair_map_mat(P[n])
        for s in active
            ds, dn = backend.layout.dims[s], backend.layout.dims[n]
            M = reshape(_slice_vee(backend.V, n, s), dn^2, ds^2)
            _add_grouped!(G, Rn, M, _pair_map_mat(P[s]), tmp)
        end
    end
    G
end

function _build_cached_state(ra::UnitRange{Int}, raV::UnitRange{Int},
    phi_in::AbstractMatrix, backend::SlicedBasisBackendCached)
    phi = Float64.(phi_in); layout, dims = backend.layout, backend.layout.dims
    m = size(phi, 2); P = _build_slice_basis(layout, ra, phi)
    active, exterior = _active_slices(layout, ra), _active_slices(layout, raV)
    pair_offs = _pair_offs(layout, exterior)
    W, Wswap = zeros(Float64, m^2, pair_offs[end]), zeros(Float64, m^2, pair_offs[end])
    for n in active
        dn = dims[n]; Rn = _pair_map_mat(P[n])
        for (j, s) in enumerate(exterior)
            ds = dims[s]
            mul!(view(W, :, _pair_range(pair_offs, j)), Rn',
                reshape(_slice_vee(backend.V, n, s), dn * dn, ds * ds), 1.0, 1.0)
            mul!(view(Wswap, :, _pair_range(pair_offs, j)), Rn',
                transpose(reshape(_slice_vee(backend.V, s, n), ds * ds, dn * dn)),
                1.0, 1.0)
        end
    end
    SlicedCachedBlockState(ra, phi, _direct_grouped(P, active, m, backend), exterior,
        pair_offs, W, Wswap)
end

_whole_slice_range(layout, ra) = first(ra) - 1 in layout.offs && last(ra) in layout.offs

function _check_cached_exterior(side, ra, raV, N)
    expected = side === :left ? (last(ra) + 1:N) : (1:first(ra) - 1)
    raV == expected || error("raV must be the complete one-sided exterior range")
end

function _aligned_absorption(side, oldra, cra, raV, layout)
    _whole_slice_range(layout, oldra) && _whole_slice_range(layout, cra) || return false
    N = layout.offs[end]
    side === :left && first(oldra) != 1 && return false
    side === :right && last(oldra) != N && return false
    side === :left && return last(oldra) + 1 == first(cra) && raV == last(cra) + 1:N
    last(cra) + 1 == first(oldra) && raV == 1:first(cra) - 1
end

function _old_basis_map(side, old::SlicedCachedBlockState, cra, phi_new)
    rows = side === :left ? (1:length(old.ra)) : (length(cra) + 1:size(phi_new, 1))
    old_rows = @view phi_new[rows, :]
    A = old.phi' * old_rows
    err = norm(old_rows - old.phi * A)
    scale = max(1.0, norm(old_rows)) * max(size(old_rows)...)
    err <= 32 * eps(Float64) * scale ? A : nothing
end

_slice_phi(layout, ra, phi, s) = view(phi, (layout.offs[s] - first(ra) + 2):(layout.offs[s + 1] - first(ra) + 1), :)

function _two_sided!(Y, X, A, work, beta = 0.0)
    tmp = view(work, 1:size(A, 1), 1:size(A, 2))
    mul!(tmp, X, A)
    mul!(Y, transpose(A), tmp, 1.0, beta)
end

function _transform_channels!(dest, source, A, source_cols, work)
    mold, mnew = size(A)
    for (j, oldj) in enumerate(source_cols)
        _two_sided!(reshape(view(dest, :, j), mnew, mnew),
            reshape(view(source, :, oldj), mold, mold), A, work)
    end
end

function _add_tensor_transform!(G, X, L, R, stage, work, input, output;
        pairtranspose = false)
    p, r = size(L); q = size(R, 1)
    size(R, 2) == r || error("tensor transforms must have the same output rank")
    S = view(stage, 1:r^2, 1:q^2)
    for col = 1:q^2
        _two_sided!(reshape(view(S, :, col), r, r),
            reshape(view(X, :, col), p, p), L, work)
    end
    Xin, Yout = view(input, 1:q, 1:q), view(output, 1:r, 1:r)
    for row = 1:r^2
        for b = 1:q, a = 1:q
            Xin[a, b] = S[row, a + (b - 1) * q]
        end
        _two_sided!(Yout, Xin, R, work)
        if pairtranspose
            for col = 1:r^2
                G[col, row] += Yout[col]
            end
        else
            for col = 1:r^2
                G[row, col] += Yout[col]
            end
        end
    end
    nothing
end

function _add_absorbed_channels!(W, Wswap, pair_offs, absorbed, exterior, newra,
        phi, backend, work, input)
    dims = backend.layout.dims
    for n in absorbed
        Cn = _slice_phi(backend.layout, newra, phi, n)
        dn = dims[n]; X = view(input, 1:dn, 1:dn)
        for (j, s) in enumerate(exterior)
            ds = dims[s]; Vns = _slice_vee(backend.V, n, s)
            Vsn = _slice_vee(backend.V, s, n)
            for b = 1:ds, a = 1:ds
                col = first(_pair_range(pair_offs, j)) + a + (b - 1) * ds - 1
                for d = 1:dn, c = 1:dn
                    X[c, d] = Vns[c, d, a, b]
                end
                _two_sided!(reshape(view(W, :, col), size(phi, 2), size(phi, 2)),
                    X, Cn, work, 1.0)
                for d = 1:dn, c = 1:dn
                    X[c, d] = Vsn[a, b, c, d]
                end
                _two_sided!(reshape(view(Wswap, :, col), size(phi, 2), size(phi, 2)),
                    X, Cn, work, 1.0)
            end
        end
    end
end

function _incremental_cached_state(side, old::SlicedCachedBlockState, cra, phi_in, A,
    raV, backend::SlicedBasisBackendCached)
    phi = Float64.(phi_in); layout, dims = backend.layout, backend.layout.dims
    newra = side === :left ? (first(old.ra):last(cra)) : (first(cra):last(old.ra))
    absorbed, exterior = _active_slices(layout, cra), _active_slices(layout, raV)
    expected = side === :left ? (first(absorbed):last(exterior)) : (first(exterior):last(absorbed))
    old.exterior_slices == expected || error("cached exterior span is inconsistent")
    mnew, mold = size(phi, 2), size(old.phi, 2); pair_offs = _pair_offs(layout, exterior)
    W, Wswap = zeros(Float64, mnew^2, pair_offs[end]), zeros(Float64, mnew^2, pair_offs[end])
    oldcols = first(_state_pair_range(old, first(exterior))):last(_state_pair_range(old, last(exterior)))
    maxdim = max(mold, maximum(dims[absorbed]))
    work = zeros(Float64, maxdim, mnew)
    input = zeros(Float64, maxdim, maxdim)
    _transform_channels!(W, old.W, A, oldcols, work)
    _transform_channels!(Wswap, old.Wswap, A, oldcols, work)
    _add_absorbed_channels!(W, Wswap, pair_offs, absorbed, exterior, newra, phi,
        backend, work, input)

    G = zeros(Float64, mnew^2, mnew^2)
    stage = zeros(Float64, mnew^2, max(mold^2, maximum(abs2.(dims[absorbed]))))
    output = zeros(Float64, mnew, mnew)
    _add_tensor_transform!(G, old.G, A, A, stage, work, input, output)
    for n in absorbed
        Cn = _slice_phi(layout, newra, phi, n)
        _add_tensor_transform!(G, view(old.W, :, _state_pair_range(old, n)),
            A, Cn, stage, work, input, output)
        _add_tensor_transform!(G, view(old.Wswap, :, _state_pair_range(old, n)),
            A, Cn, stage, work, input, output; pairtranspose = true)
        for s in absorbed
            Cs = _slice_phi(layout, newra, phi, s)
            M = reshape(_slice_vee(backend.V, n, s), dims[n]^2, dims[s]^2)
            _add_tensor_transform!(G, M, Cn, Cs, stage, work, input, output)
        end
    end
    SlicedCachedBlockState(newra, phi, G, exterior, pair_offs, W, Wswap)
end

function vee_init_block(side, ra, raV, phi, backend::SlicedBasisBackendCached)
    _check_side(side)
    _check_range(ra, "ra")
    _check_range(raV, "raV")
    _check_cached_exterior(side, ra, raV, backend.layout.offs[end])
    size(phi, 1) == length(ra) || error("phi has wrong row count for ra")
    _bench_timing_enabled() && (_bench_timing[:cache_init_calls] += 1)
    _build_cached_state(ra, raV, phi, backend)
end

function vee_absorb_block(side, vee_old::SlicedCachedBlockState, cra, Phi_old, Phi_C,
    phi_new, raV_new, backend::SlicedBasisBackendCached)
    _check_side(side)
    _check_range(cra, "cra")
    _check_range(raV_new, "raV_new")
    oldra = vee_old.ra
    newra = side == :left ? (oldra[1]:cra[end]) : (cra[1]:oldra[end])
    _check_cached_exterior(side, newra, raV_new, backend.layout.offs[end])
    size(phi_new, 1) == length(newra) || error("phi_new has wrong row count for new ra")
    aligned = _aligned_absorption(side, oldra, cra, raV_new, backend.layout)
    A = aligned ? _old_basis_map(side, vee_old, cra, phi_new) : nothing
    if A === nothing
        _bench_timing_enabled() && (_bench_timing[:cache_fallback_calls] += 1)
        return _build_cached_state(newra, raV_new, phi_new, backend)
    end
    _bench_timing_enabled() && (_bench_timing[:cache_incremental_calls] += 1)
    _incremental_cached_state(side, vee_old, cra, phi_new, A, raV_new, backend)
end

function vee_window(Lvee::SlicedCachedBlockState, Rvee::SlicedCachedBlockState,
    Cra, backend::SlicedBasisBackendCached)
    _check_range(Cra, "Cra")
    Lra = Lvee.ra
    Rra = Rvee.ra
    Lra[end] < Cra[1] || error("Lra must end before Cra")
    Cra[end] < Rra[1] || error("Cra must end before Rra")
    N = backend.layout.offs[end]
    Rra[end] <= N || error("ranges exceed full basis size")

    aligned = first(Lra) == 1 && last(Lra) + 1 == first(Cra) &&
        last(Cra) + 1 == first(Rra) && last(Rra) == N &&
        _whole_slice_range(backend.layout, Lra) &&
        _whole_slice_range(backend.layout, Cra) &&
        _whole_slice_range(backend.layout, Rra)
    if !aligned
        _bench_timing_enabled() && (_bench_timing[:window_projection_calls] += 1)
        projection = SlicedBasisBackend(backend.layout, backend.V)
        return vee_window(SlicedBlockState(Lra, Lvee.phi),
            SlicedBlockState(Rra, Rvee.phi), Cra, projection)
    end
    _bench_timing_enabled() && (_bench_timing[:window_local_calls] += 1)

    ml = size(Lvee.phi, 2)
    mr = size(Rvee.phi, 2)
    lc = length(Cra)
    dims = backend.layout.dims
    center_slices = (first(_slice_local(backend.layout, first(Cra))):
        first(_slice_local(backend.layout, last(Cra))))

    maxd = maximum(dims)
    maxk = max(ml, mr)
    scratch = SlicedCachedCrossScratch(zeros(maxd, maxd), zeros(maxd, maxk),
        zeros(maxk, maxd), zeros(maxk, maxd), zeros(maxk^2))

    SlicedBasisCachedWindow(backend.V, backend.layout, ml, lc, mr, center_slices,
        Lvee, Rvee, scratch)
end

# Pair-grouped tensors store T[x,z;y,w]; Wswap reverses the grouped region pairs.
@inline _sector_value(V, x, z, y, w, nx, ny, ::Val{:pair_grouped}) = V[x, z, y, w]
@inline _sector_value(V, x, z, y, w, nx, ny, ::Val{:pair_matrix}) = V[x + (z - 1) * nx, y + (w - 1) * ny]
@inline _sector_value(V, x, z, y, w, nx, ny, ::Val{:pair_swapped}) = V[y, w, x, z]

function _add_cached_sector!(Fup, Fdn, rhoup, rhodn, V, xoff, yoff, nx, ny,
    ordering, ::Val{restricted}) where {restricted}
    @inbounds for w = 1:ny, y = 1:ny, z = 1:nx, x = 1:nx
        px, pz = xoff + x, xoff + z
        py, pw = yoff + y, yoff + w
        v = _sector_value(V, x, z, y, w, nx, ny, ordering)
        if restricted
            Fup[px, pz] += 2.0 * rhoup[py, pw] * v
            Fup[px, pw] -= rhoup[pz, py] * v
        else
            direct = (rhoup[py, pw] + rhodn[py, pw]) * v
            Fup[px, pz] += direct
            Fdn[px, pz] += direct
            Fup[px, pw] -= rhoup[pz, py] * v
            Fdn[px, pw] -= rhodn[pz, py] * v
        end
    end
    nothing
end

function _add_cross_direct!(Fup, Fdn, rhoup, rhodn, X, Y, xoff, yoff,
        layout, scratch, ::Val{restricted}) where {restricted}
    kx, ky = size(X.phi, 2), size(Y.phi, 2)
    yrange = yoff + 1:yoff + ky
    for s in _active_slices(layout, Y.ra)
        P = _slice_phi(layout, Y.ra, Y.phi, s); d = size(P, 1)
        density = view(scratch.density, 1:d, 1:d)
        work = view(scratch.physical_work, 1:d, 1:ky)
        mul!(work, P, view(rhoup, yrange, yrange))
        mul!(density, work, transpose(P))
        if !restricted
            mul!(work, P, view(rhodn, yrange, yrange))
            mul!(density, work, transpose(P), 1.0, 1.0)
        end
        W = _packed_tensor(X.W, X, s, kx, d)
        pair = view(scratch.pair_field, 1:kx^2)
        for z = 1:kx, x = 1:kx
            value = 0.0
            for b = 1:d, a = 1:d
                value = muladd(W[x, z, a, b], density[a, b], value)
            end
            pair[x + (z - 1) * kx] = value
        end
        factor = restricted ? 2.0 : 1.0
        for z = 1:kx, x = 1:kx
            value = factor * pair[x + (z - 1) * kx]
            Fup[xoff + x, xoff + z] += value
            restricted || (Fdn[xoff + x, xoff + z] += value)
        end
    end
    nothing
end

function _add_cross_exchange!(F, rho, X, Y, xoff, yoff, layout, scratch)
    kx, ky = size(X.phi, 2), size(Y.phi, 2)
    xrange, yrange = xoff + 1:xoff + kx, yoff + 1:yoff + ky
    for s in _active_slices(layout, Y.ra)
        P = _slice_phi(layout, Y.ra, Y.phi, s); d = size(P, 1)
        mixed = view(scratch.mixed, 1:kx, 1:d)
        field = view(scratch.exchange, 1:kx, 1:d)
        mul!(mixed, view(rho, xrange, yrange), transpose(P))
        W = _packed_tensor(X.W, X, s, kx, d)
        for b = 1:d, i = 1:kx
            value = 0.0
            for a = 1:d, k = 1:kx
                value = muladd(W[i, k, a, b], mixed[k, a], value)
            end
            field[i, b] = value
        end
        mul!(view(F, xrange, yrange), field, P, -1.0, 1.0)
    end
    nothing
end

function _add_cross!(Fup, Fdn, rhoup, rhodn, X, Y, xoff, yoff, layout,
        scratch, restricted::Val{r}) where r
    _add_cross_direct!(Fup, Fdn, rhoup, rhodn, X, Y, xoff, yoff, layout,
        scratch, restricted)
    _add_cross_exchange!(Fup, rhoup, X, Y, xoff, yoff, layout, scratch)
    r ||
        _add_cross_exchange!(Fdn, rhodn, X, Y, xoff, yoff, layout, scratch)
end

function _add_cached_fock!(Fup, Fdn, rhoup, rhodn, win, restricted)
    ml, lc, mr = win.ml, win.lc, win.mr
    roff = ml + lc
    dims = win.layout.dims
    grouped = Val(:pair_grouped)
    swapped = Val(:pair_swapped)
    matrix = Val(:pair_matrix)

    _add_cached_sector!(Fup, Fdn, rhoup, rhodn, win.L.G, 0, 0, ml, ml,
        matrix, restricted)
    coff = ml
    for s in win.center_slices
        ds = dims[s]
        _add_cached_sector!(Fup, Fdn, rhoup, rhodn,
            _packed_tensor(win.L.W, win.L, s, ml, ds), 0, coff, ml, ds,
            grouped, restricted)
        coff += ds
    end
    _add_cross!(Fup, Fdn, rhoup, rhodn, win.L, win.R, 0, roff, win.layout,
        win.scratch, restricted)

    coff = ml
    for s in win.center_slices
        ds = dims[s]
        _add_cached_sector!(Fup, Fdn, rhoup, rhodn,
            _packed_tensor(win.L.Wswap, win.L, s, ml, ds), coff, 0, ds,
            ml, swapped, restricted)
        coff += ds
    end
    xoff = ml
    for sx in win.center_slices
        yoff = ml
        for sy in win.center_slices
            _add_cached_sector!(Fup, Fdn, rhoup, rhodn, _slice_vee(win.V, sx, sy),
                xoff, yoff, dims[sx], dims[sy], grouped, restricted)
            yoff += dims[sy]
        end
        xoff += dims[sx]
    end
    coff = ml
    for s in win.center_slices
        ds = dims[s]
        _add_cached_sector!(Fup, Fdn, rhoup, rhodn,
            _packed_tensor(win.R.Wswap, win.R, s, mr, ds), coff, roff, ds,
            mr, swapped, restricted)
        coff += ds
    end

    _add_cross!(Fup, Fdn, rhoup, rhodn, win.R, win.L, roff, 0, win.layout,
        win.scratch, restricted)
    coff = ml
    for s in win.center_slices
        ds = dims[s]
        _add_cached_sector!(Fup, Fdn, rhoup, rhodn,
            _packed_tensor(win.R.W, win.R, s, mr, ds), roff, coff, mr, ds,
            grouped, restricted)
        coff += ds
    end
    _add_cached_sector!(Fup, Fdn, rhoup, rhodn, win.R.G, roff, roff, mr, mr,
        matrix, restricted)
    nothing
end

function vee_add_fock_r!(F, rho, win::SlicedBasisCachedWindow)
    n = size(F, 1)
    size(F) == (n, n) || error("F has wrong size")
    size(rho) == size(F) || error("rho has wrong size")
    win.ml + win.lc + win.mr == n || error("window dimensions do not match F")
    _add_cached_fock!(F, F, rho, rho, win, Val(true))
end

function vee_add_fock!(Fup, Fdn, rhoup, rhodn, win::SlicedBasisCachedWindow)
    n = size(Fup, 1)
    size(Fup) == (n, n) || error("Fup has wrong size")
    size(Fdn) == (n, n) || error("Fdn has wrong size")
    size(rhoup) == size(Fup) || error("rhoup has wrong size")
    size(rhodn) == size(Fup) || error("rhodn has wrong size")
    win.ml + win.lc + win.mr == n || error("window dimensions do not match F")
    _add_cached_fock!(Fup, Fdn, rhoup, rhodn, win, Val(false))
end

"""
solve_hfdmrg(H, backend::SlicedBasisBackendCached, psiup0; kwargs...) -> (psiup, psidn, energy)

Restricted HF (RHF) sweep using the cached sliced-basis backend.
"""
function solve_hfdmrg(H, backend::SlicedBasisBackendCached, psiup0; kwargs...)
    solve_hfdmrg_core(H, backend, psiup0, psiup0; restricted = true, kwargs...)
end

"""
solve_hfdmrg(H, backend::SlicedBasisBackendCached, psiup0, psidn0; kwargs...) -> (psiup, psidn, energy)

Unrestricted HF (UHF) sweep using the cached sliced-basis backend.
"""
function solve_hfdmrg(H, backend::SlicedBasisBackendCached, psiup0, psidn0; kwargs...)
    solve_hfdmrg_core(H, backend, psiup0, psidn0; restricted = false, kwargs...)
end

"""
solve_hfdmrg(Hup, Hdn, backend::SlicedBasisBackendCached, psiup0, psidn0; kwargs...) -> (psiup, psidn, energy)

Unrestricted HF (UHF) sweep with spin-dependent one-body Hamiltonians and the
cached sliced-basis backend. If Hup and Hdn are the same matrix object, this
routes to the common-H fast path.
"""
function solve_hfdmrg(Hup, Hdn, backend::SlicedBasisBackendCached, psiup0, psidn0;
    kwargs...)
    Hup === Hdn && return solve_hfdmrg(Hup, backend, psiup0, psidn0; kwargs...)
    solve_hfdmrg_core_split(Hup, Hdn, backend, psiup0, psidn0; kwargs...)
end
