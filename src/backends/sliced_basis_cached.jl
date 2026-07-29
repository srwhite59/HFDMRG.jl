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

struct SlicedBasisCachedWindow{T}
    V::T
    layout::SliceLayout
    ml::Int
    lc::Int
    mr::Int
    center_slices::UnitRange{Int}
    L::SlicedCachedBlockState
    R::SlicedCachedBlockState
    VLR::Array{Float64,4}
    VRL::Array{Float64,4}
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

function _incremental_cached_state(side, old::SlicedCachedBlockState, cra, phi_in, A,
    raV, backend::SlicedBasisBackendCached)
    phi = Float64.(phi_in); layout, dims = backend.layout, backend.layout.dims
    newra = side === :left ? (first(old.ra):last(cra)) : (first(cra):last(old.ra))
    absorbed, exterior = _active_slices(layout, cra), _active_slices(layout, raV)
    expected = side === :left ? (first(absorbed):last(exterior)) : (first(exterior):last(absorbed))
    old.exterior_slices == expected || error("cached exterior span is inconsistent")
    mnew, mold = size(phi, 2), size(old.phi, 2); pair_offs = _pair_offs(layout, exterior)
    W, Wswap = zeros(Float64, mnew^2, pair_offs[end]), zeros(Float64, mnew^2, pair_offs[end])
    RA = _pair_map_mat(A)
    oldcols = first(_state_pair_range(old, first(exterior))):last(_state_pair_range(old, last(exterior)))
    mul!(W, RA', view(old.W, :, oldcols))
    mul!(Wswap, RA', view(old.Wswap, :, oldcols))
    Rnew = [_pair_map_mat(_slice_phi(layout, newra, phi, n)) for n in absorbed]
    for (i, n) in enumerate(absorbed)
        dn = dims[n]
        for (j, s) in enumerate(exterior)
            ds = dims[s]
            mul!(view(W, :, _pair_range(pair_offs, j)), Rnew[i]',
                reshape(_slice_vee(backend.V, n, s), dn * dn, ds * ds), 1.0, 1.0)
            mul!(view(Wswap, :, _pair_range(pair_offs, j)), Rnew[i]',
                transpose(reshape(_slice_vee(backend.V, s, n), ds * ds, dn * dn)),
                1.0, 1.0)
        end
    end
    G = zeros(Float64, mnew^2, mnew^2)
    tmp = zeros(Float64, mnew^2, max(mold^2, maximum(abs2.(dims[absorbed]))))
    _add_grouped!(G, RA, old.G, RA, tmp)
    for (i, n) in enumerate(absorbed)
        Rn = Rnew[i]; _add_grouped!(
            G, RA, view(old.W, :, _state_pair_range(old, n)), Rn, tmp)
        _add_grouped!(G, Rn,
            transpose(view(old.Wswap, :, _state_pair_range(old, n))), RA, tmp)
        for (j, s) in enumerate(absorbed)
            M = reshape(_slice_vee(backend.V, n, s), dims[n]^2, dims[s]^2)
            _add_grouped!(G, Rn, M, Rnew[j], tmp)
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

    VLR = zeros(Float64, ml, ml, mr, mr)
    VRL = zeros(Float64, mr, mr, ml, ml)
    VLR_mat = reshape(VLR, ml * ml, mr * mr)
    VRL_mat = reshape(VRL, mr * mr, ml * ml)
    maxd = maximum(dims)
    Utmp_R = zeros(Float64, maxd * maxd, mr * mr)
    Utmp_L = zeros(Float64, maxd * maxd, ml * ml)
    for s in _active_slices(backend.layout, Rra)
        ds = dims[s]
        Uview = view(Utmp_R, 1:(ds * ds), :)
        _pair_map_mat!(Uview, _slice_phi(backend.layout, Rra, Rvee.phi, s))
        mul!(VLR_mat, view(Lvee.W, :, _state_pair_range(Lvee, s)),
            Uview, 1.0, 1.0)
    end
    for s in _active_slices(backend.layout, Lra)
        ds = dims[s]
        Uview = view(Utmp_L, 1:(ds * ds), :)
        _pair_map_mat!(Uview, _slice_phi(backend.layout, Lra, Lvee.phi, s))
        mul!(VRL_mat, view(Rvee.W, :, _state_pair_range(Rvee, s)),
            Uview, 1.0, 1.0)
    end

    SlicedBasisCachedWindow(backend.V, backend.layout, ml, lc, mr, center_slices,
        Lvee, Rvee, VLR, VRL)
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
    _add_cached_sector!(Fup, Fdn, rhoup, rhodn, win.VLR, 0, roff, ml, mr,
        grouped, restricted)

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

    _add_cached_sector!(Fup, Fdn, rhoup, rhodn, win.VRL, roff, 0, mr, ml,
        grouped, restricted)
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
