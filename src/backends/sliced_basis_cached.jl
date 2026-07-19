"""
Cached sliced-basis backend (correctness-first).

Uses per-slice basis coefficients to build superblock Fock contributions without
embedding into the full orbital basis.
"""
struct SlicedBasisBackendCached{T}
    layout::SliceLayout
    V::T
    ns::Int
end

struct SlicedCachedBlockState
    ra::UnitRange{Int}
    phi::Matrix{Float64}
    Vijkl::Array{Float64,4}
    W::Vector{Array{Float64,4}}
    # Wswap[s][i,j,a,b] uses the swapped orientation V6[a,b,c,d,s,n] to avoid
    # assuming additional symmetries when building block-center exchange terms.
    Wswap::Vector{Array{Float64,4}}
    P::Vector{Matrix{Float64}}
end

struct SlicedBasisCachedWindow{T}
    S::Vector{Matrix{Float64}}
    V::T
    layout::SliceLayout
    ns::Int
    ml::Int
    lc::Int
    mr::Int
    split_slices::Vector{Int}
    VLL::Array{Float64,4}
    VRR::Array{Float64,4}
    WL::Vector{Array{Float64,4}}
    WR::Vector{Array{Float64,4}}
    WLswap::Vector{Array{Float64,4}}
    WRswap::Vector{Array{Float64,4}}
    VLR::Array{Float64,4}
    VRL::Array{Float64,4}
    center_slice::Vector{Int}
    center_local::Vector{Int}
    center_cols::Vector{Vector{Int}}
    Srho::Vector{Matrix{Float64}}
    Srho_up::Vector{Matrix{Float64}}
    Srho_dn::Vector{Matrix{Float64}}
    rho_slice::Vector{Matrix{Float64}}
    rho_slice_up::Vector{Matrix{Float64}}
    rho_slice_dn::Vector{Matrix{Float64}}
    K::Vector{Matrix{Float64}}
    tmp_n::Vector{Matrix{Float64}}
end

const _bench_timing = Dict{Symbol, Float64}(
    :direct_LL_RR => 0.0,
    :direct_CL_CR => 0.0,
    :direct_LR => 0.0,
    :exchange => 0.0,
    :cache_init_calls => 0.0,
    :cache_incremental_calls => 0.0,
    :cache_fallback_calls => 0.0,
)

# Keep the historical timing keys for compatibility with scripts/bench_sliced.jl.

_bench_timing_enabled() = !isempty(get(ENV, "HFDMRG_BENCH_TIMING", ""))

function _bench_timing_reset!()
    for key in keys(_bench_timing)
        _bench_timing[key] = 0.0
    end
    nothing
end

_bench_timing_snapshot() = Dict(key => value for (key, value) in _bench_timing)

function SlicedBasisBackendCached(layout::SliceLayout, V6::Array{Float64,6})
    ns = nslices(layout)
    nj = layout.dims[1]
    any(d -> d != nj, layout.dims) && error("layout must have fixed slice size")
    size(V6) == (nj, nj, nj, nj, ns, ns) || error("V6 must have size (nj,nj,nj,nj,ns,ns)")
    SlicedBasisBackendCached{Array{Float64,6}}(layout, V6, ns)
end

function SlicedBasisBackendCached(layout::SliceLayout, vee::SlicedVeeRagged)
    layout.dims == vee.layout.dims || error("layout does not match ragged interaction")
    ns = nslices(layout)
    SlicedBasisBackendCached{SlicedVeeRagged}(layout, vee, ns)
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

function _pair_map_mat!(U::AbstractMatrix{Float64}, P::Matrix{Float64})
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

function _pair_map_mat(P::Matrix{Float64})
    nj, m = size(P)
    U = zeros(Float64, nj * nj, m * m)
    _pair_map_mat!(U, P)
end

function _active_slices(layout::SliceLayout, ra::UnitRange{Int})
    first(_slice_local(layout, first(ra))):first(_slice_local(layout, last(ra)))
end

function _cached_vijkl(W, P, active, m)
    grouped = zeros(Float64, m * m, m * m)
    for s in active
        ds = size(P[s], 1)
        mul!(grouped, reshape(W[s], m * m, ds * ds), _pair_map_mat(P[s]), 1.0, 1.0)
    end
    permutedims(reshape(grouped, m, m, m, m), (1, 3, 2, 4))
end

function _build_cached_state(ra::UnitRange{Int}, phi_in::AbstractMatrix,
    backend::SlicedBasisBackendCached)
    phi = Float64.(phi_in)
    layout = backend.layout
    dims = layout.dims
    ns = backend.ns
    m = size(phi, 2)
    P = _build_slice_basis(layout, ra, phi)
    active = _active_slices(layout, ra)
    W = [zeros(Float64, m, m, dims[s], dims[s]) for s = 1:ns]
    Wswap = [zeros(Float64, m, m, dims[s], dims[s]) for s = 1:ns]
    for n in active
        dn = dims[n]
        Rn = _pair_map_mat(P[n])
        for s = 1:ns
            ds = dims[s]
            mul!(reshape(W[s], m * m, ds * ds), Rn',
                reshape(_slice_vee(backend.V, n, s), dn * dn, ds * ds), 1.0, 1.0)
            mul!(reshape(Wswap[s], m * m, ds * ds), Rn',
                transpose(reshape(_slice_vee(backend.V, s, n), ds * ds, dn * dn)),
                1.0, 1.0)
        end
    end
    Vijkl = _cached_vijkl(W, P, active, m)
    SlicedCachedBlockState(ra, phi, Vijkl, W, Wswap, P)
end

_whole_slice_range(layout, ra) = first(ra) - 1 in layout.offs && last(ra) in layout.offs

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

function _incremental_cached_state(side, old::SlicedCachedBlockState, cra, phi_in, A,
    backend::SlicedBasisBackendCached)
    phi = Float64.(phi_in)
    layout, dims, ns = backend.layout, backend.layout.dims, backend.ns
    newra = side === :left ? (first(old.ra):last(cra)) : (first(cra):last(old.ra))
    mnew, mold = size(phi, 2), size(old.phi, 2)
    P = _build_slice_basis(layout, newra, phi)
    W = [zeros(Float64, mnew, mnew, dims[s], dims[s]) for s = 1:ns]
    Wswap = [zeros(Float64, mnew, mnew, dims[s], dims[s]) for s = 1:ns]
    RA = _pair_map_mat(A)
    for s = 1:ns
        ds = dims[s]
        mul!(reshape(W[s], mnew * mnew, ds * ds), RA',
            reshape(old.W[s], mold * mold, ds * ds))
        mul!(reshape(Wswap[s], mnew * mnew, ds * ds), RA',
            reshape(old.Wswap[s], mold * mold, ds * ds))
    end
    for n in _active_slices(layout, cra)
        dn = dims[n]
        Rn = _pair_map_mat(P[n])
        for s = 1:ns
            ds = dims[s]
            mul!(reshape(W[s], mnew * mnew, ds * ds), Rn',
                reshape(_slice_vee(backend.V, n, s), dn * dn, ds * ds), 1.0, 1.0)
            mul!(reshape(Wswap[s], mnew * mnew, ds * ds), Rn',
                transpose(reshape(_slice_vee(backend.V, s, n), ds * ds, dn * dn)),
                1.0, 1.0)
        end
    end
    active = _active_slices(layout, newra)
    Vijkl = _cached_vijkl(W, P, active, mnew)
    SlicedCachedBlockState(newra, phi, Vijkl, W, Wswap, P)
end

function vee_init_block(side, ra, raV, phi, backend::SlicedBasisBackendCached)
    _check_side(side)
    _check_range(ra, "ra")
    _check_range(raV, "raV")
    size(phi, 1) == length(ra) || error("phi has wrong row count for ra")
    _bench_timing_enabled() && (_bench_timing[:cache_init_calls] += 1)
    _build_cached_state(ra, phi, backend)
end

function vee_absorb_block(side, vee_old::SlicedCachedBlockState, cra, Phi_old, Phi_C,
    phi_new, raV_new, backend::SlicedBasisBackendCached)
    _check_side(side)
    _check_range(cra, "cra")
    _check_range(raV_new, "raV_new")
    oldra = vee_old.ra
    newra = side == :left ? (oldra[1]:cra[end]) : (cra[1]:oldra[end])
    size(phi_new, 1) == length(newra) || error("phi_new has wrong row count for new ra")
    aligned = _aligned_absorption(side, oldra, cra, raV_new, backend.layout)
    A = aligned ? _old_basis_map(side, vee_old, cra, phi_new) : nothing
    if A === nothing
        _bench_timing_enabled() && (_bench_timing[:cache_fallback_calls] += 1)
        return _build_cached_state(newra, phi_new, backend)
    end
    _bench_timing_enabled() && (_bench_timing[:cache_incremental_calls] += 1)
    _incremental_cached_state(side, vee_old, cra, phi_new, A, backend)
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

    ml = size(Lvee.phi, 2)
    mr = size(Rvee.phi, 2)
    lc = length(Cra)
    superdim = ml + lc + mr
    ns = backend.ns
    dims = backend.layout.dims
    S = [zeros(Float64, dims[s], superdim) for s = 1:ns]
    for s = 1:ns
        S[s][:, 1:ml] = Lvee.P[s]
        S[s][:, ml + lc + 1:end] = Rvee.P[s]
    end
    for (idx, p) in enumerate(Cra)
        n, a = _slice_local(backend.layout, p)
        S[n][a, ml + idx] = 1.0
    end
    center_slice = zeros(Int, lc)
    center_local = zeros(Int, lc)
    center_cols = [fill(0, dims[s]) for s in 1:ns]
    for (idx, p) in enumerate(Cra)
        n, a = _slice_local(backend.layout, p)
        center_slice[idx] = n
        center_local[idx] = a
        center_cols[n][a] = ml + idx
    end

    split_slices = Int[]
    for s = 1:ns
        has_center = any(c -> c != 0, center_cols[s])
        has_block = any(!iszero, Lvee.P[s]) || any(!iszero, Rvee.P[s])
        if has_center && has_block
            push!(split_slices, s)
        end
    end

    Srho = [zeros(Float64, dims[s], superdim) for s = 1:ns]
    Srho_up = [zeros(Float64, dims[s], superdim) for s = 1:ns]
    Srho_dn = [zeros(Float64, dims[s], superdim) for s = 1:ns]
    rho_slice = [zeros(Float64, dims[s], dims[s]) for s = 1:ns]
    rho_slice_up = [zeros(Float64, dims[s], dims[s]) for s = 1:ns]
    rho_slice_dn = [zeros(Float64, dims[s], dims[s]) for s = 1:ns]
    K = [zeros(Float64, dims[s], dims[s]) for s = 1:ns]
    tmp_n = [zeros(Float64, dims[s], superdim) for s = 1:ns]

    VLR = zeros(Float64, ml, ml, mr, mr)
    VRL = zeros(Float64, mr, mr, ml, ml)
    VLR_mat = reshape(VLR, ml * ml, mr * mr)
    VRL_mat = reshape(VRL, mr * mr, ml * ml)
    maxd = maximum(dims)
    Utmp_R = zeros(Float64, maxd * maxd, mr * mr)
    Utmp_L = zeros(Float64, maxd * maxd, ml * ml)
    for s = 1:ns
        ds = dims[s]
        if any(!iszero, Rvee.P[s])
            Uview = view(Utmp_R, 1:(ds * ds), :)
            _pair_map_mat!(Uview, Rvee.P[s])
            WL_mat = reshape(Lvee.W[s], ml * ml, ds * ds)
            mul!(VLR_mat, WL_mat, Uview, 1.0, 1.0)
        end
        if any(!iszero, Lvee.P[s])
            Uview = view(Utmp_L, 1:(ds * ds), :)
            _pair_map_mat!(Uview, Lvee.P[s])
            WR_mat = reshape(Rvee.W[s], mr * mr, ds * ds)
            mul!(VRL_mat, WR_mat, Uview, 1.0, 1.0)
        end
    end

    SlicedBasisCachedWindow(S, backend.V, backend.layout, ns, ml, lc, mr, split_slices,
        Lvee.Vijkl, Rvee.Vijkl, Lvee.W, Rvee.W, Lvee.Wswap, Rvee.Wswap, VLR, VRL,
        center_slice, center_local, center_cols, Srho, Srho_up, Srho_dn, rho_slice,
        rho_slice_up, rho_slice_dn, K, tmp_n)
end

function vee_add_fock_r!(F, rho, win::SlicedBasisCachedWindow)
    S = win.S
    V = win.V
    dims = win.layout.dims
    ns = win.ns
    ml, lc, mr = win.ml, win.lc, win.mr
    n = size(F, 1)
    size(F) == (n, n) || error("F has wrong size")
    size(rho) == size(F) || error("rho has wrong size")
    ml + lc + mr == n || error("window dimensions do not match F")

    center_cols = win.center_cols
    Srho = win.Srho
    @views for n1 = 1:ns
        mul!(Srho[n1], S[n1], rho, 1.0, 0.0)
    end
    split_slices = win.split_slices
    timing = _bench_timing_enabled()
    if isempty(split_slices)
        # Exchange (K) contributions via cached tensors.
        VLL = win.VLL
        if timing
            t0 = time_ns()
        end
        for i = 1:ml, l = 1:ml
            acc = 0.0
            for j = 1:ml, k = 1:ml
                acc += rho[k, j] * VLL[i, j, k, l]
            end
            F[i, l] -= acc
        end

        VRR = win.VRR
        for i = 1:mr, l = 1:mr
            acc = 0.0
            for j = 1:mr, k = 1:mr
                acc += rho[ml + lc + k, ml + lc + j] * VRR[i, j, k, l]
            end
            F[ml + lc + i, ml + lc + l] -= acc
        end
        if timing
            _bench_timing[:direct_LL_RR] += (time_ns() - t0) * 1e-9
        end

        if timing
            t0 = time_ns()
        end
        for np = 1:ns
            cols_p = center_cols[np]
            any(c -> c != 0, cols_p) || continue
            for nq = 1:ns
                cols_q = center_cols[nq]
                any(c -> c != 0, cols_q) || continue
                rho_nm = Matrix{Float64}(undef, length(cols_p), length(cols_q))
                mul!(rho_nm, Srho[np], S[nq]', 1.0, 0.0)
                Vpq = _slice_vee(V, np, nq)
                for a = 1:length(cols_p)
                    col_a = cols_p[a]
                    col_a == 0 && continue
                    for d = 1:length(cols_q)
                        col_d = cols_q[d]
                        col_d == 0 && continue
                        acc = 0.0
                        for b = 1:length(cols_p), c = 1:length(cols_q)
                            acc += rho_nm[b, c] * Vpq[a, b, c, d]
                        end
                        F[col_a, col_d] -= acc
                    end
                end
            end
        end

        WL = win.WL
        WLswap = win.WLswap
        for s = 1:ns
            cols = center_cols[s]
            for b = 1:length(cols)
                col_b = cols[b]
                col_b == 0 && continue
                for i = 1:ml
                    acc = 0.0
                    for a = 1:length(cols)
                        col_a = cols[a]
                        col_a == 0 && continue
                        for j = 1:ml
                            acc += rho[j, col_a] * WL[s][i, j, a, b]
                        end
                    end
                    F[i, col_b] -= acc
                end
            end
            for a = 1:length(cols)
                col_a = cols[a]
                col_a == 0 && continue
                for j = 1:ml
                    acc = 0.0
                    for b = 1:length(cols)
                        col_b = cols[b]
                        col_b == 0 && continue
                        for i = 1:ml
                            acc += rho[col_b, i] * WLswap[s][i, j, a, b]
                        end
                    end
                    F[col_a, j] -= acc
                end
            end
        end

        WR = win.WR
        WRswap = win.WRswap
        for s = 1:ns
            cols = center_cols[s]
            for b = 1:length(cols)
                col_b = cols[b]
                col_b == 0 && continue
                for i = 1:mr
                    acc = 0.0
                    for a = 1:length(cols)
                        col_a = cols[a]
                        col_a == 0 && continue
                        for j = 1:mr
                            acc += rho[ml + lc + j, col_a] * WR[s][i, j, a, b]
                        end
                    end
                    F[ml + lc + i, col_b] -= acc
                end
            end
            for a = 1:length(cols)
                col_a = cols[a]
                col_a == 0 && continue
                for j = 1:mr
                    acc = 0.0
                    for b = 1:length(cols)
                        col_b = cols[b]
                        col_b == 0 && continue
                        for i = 1:mr
                            acc += rho[col_b, ml + lc + i] * WRswap[s][i, j, a, b]
                        end
                    end
                    F[col_a, ml + lc + j] -= acc
                end
            end
        end
        if timing
            _bench_timing[:direct_CL_CR] += (time_ns() - t0) * 1e-9
        end

        VLR = win.VLR
        VRL = win.VRL

        if timing
            t0 = time_ns()
        end
        for i = 1:ml, l = 1:mr
            acc = 0.0
            for j = 1:ml, k = 1:mr
                acc += rho[j, ml + lc + k] * VLR[i, j, k, l]
            end
            F[i, ml + lc + l] -= acc
        end
        for k = 1:mr, j = 1:ml
            acc = 0.0
            for l = 1:mr, i = 1:ml
                acc += rho[ml + lc + l, i] * VRL[k, l, i, j]
            end
            F[ml + lc + k, j] -= acc
        end
        if timing
            _bench_timing[:direct_LR] += (time_ns() - t0) * 1e-9
        end
    else
        # Correctness-first path: when any slice is split, exchange falls back to
        # full slice-projection to maintain parity with the projection backend.
        # This is intentionally non-local.
        if timing
            t0 = time_ns()
        end
        superdim = size(F, 1)
        for n1 = 1:ns
            dn = dims[n1]
            for m = 1:ns
                dm = dims[m]
                rho_nm = Matrix{Float64}(undef, dn, dm)
                mul!(rho_nm, Srho[n1], S[m]', 1.0, 0.0)
                Vnm = _slice_vee(V, n1, m)
                Knm = Matrix{Float64}(undef, dn, dm)
                for a = 1:dn, d = 1:dm
                    acc = 0.0
                    for b = 1:dn, c = 1:dm
                        acc += rho_nm[b, c] * Vnm[a, b, c, d]
                    end
                    Knm[a, d] = acc
                end
                tmp = Matrix{Float64}(undef, dn, superdim)
                mul!(tmp, Knm, S[m], 1.0, 0.0)
                mul!(F, S[n1]', tmp, -1.0, 1.0)
            end
        end
        if timing
            _bench_timing[:direct_CL_CR] += (time_ns() - t0) * 1e-9
        end
    end

    # J[p,r] is block diagonal in the output slice pair (p,r).
    if timing
        t0 = time_ns()
    end
    rho_slice = win.rho_slice
    J = win.K
    tmp_n = win.tmp_n
    @views for n1 = 1:ns
        mul!(rho_slice[n1], Srho[n1], S[n1]', 1.0, 0.0)
    end
    @views for n1 = 1:ns
        J1 = J[n1]
        Jvec = reshape(J1, :)
        fill!(Jvec, 0.0)
        dn1 = dims[n1]
        dn1sq = dn1 * dn1
        for n2 = 1:ns
            dn2 = dims[n2]
            Vnm = _slice_vee(V, n1, n2)
            Vpair = reshape(Vnm, dn1sq, dn2 * dn2)
            rvec = reshape(rho_slice[n2], dn2 * dn2)
            mul!(Jvec, Vpair, rvec, 1.0, 1.0)
        end
        tmp = tmp_n[n1]
        mul!(tmp, J1, S[n1], 1.0, 0.0)
        mul!(F, S[n1]', tmp, 2.0, 1.0)
    end
    if timing
        _bench_timing[:exchange] += (time_ns() - t0) * 1e-9
    end

end

function vee_add_fock!(Fup, Fdn, rhoup, rhodn, win::SlicedBasisCachedWindow)
    S = win.S
    V = win.V
    dims = win.layout.dims
    ns = win.ns
    ml, lc, mr = win.ml, win.lc, win.mr
    n = size(Fup, 1)
    size(Fup) == (n, n) || error("Fup has wrong size")
    size(Fdn) == (n, n) || error("Fdn has wrong size")
    size(rhoup) == size(Fup) || error("rhoup has wrong size")
    size(rhodn) == size(Fup) || error("rhodn has wrong size")
    ml + lc + mr == n || error("window dimensions do not match F")

    center_cols = win.center_cols
    Srho = win.Srho
    Srho_up = win.Srho_up
    Srho_dn = win.Srho_dn
    VLR = win.VLR
    VRL = win.VRL

    @views for n1 = 1:ns
        S1 = S[n1]
        mul!(Srho_up[n1], S1, rhoup, 1.0, 0.0)
        mul!(Srho_dn[n1], S1, rhodn, 1.0, 0.0)
        Srho[n1] .= Srho_up[n1]
        Srho[n1] .+= Srho_dn[n1]
    end
    split_slices = win.split_slices
    timing = _bench_timing_enabled()
    if isempty(split_slices)
        # Exchange (K) contributions via cached tensors.
        VLL = win.VLL
        if timing
            t0 = time_ns()
        end
        for i = 1:ml, l = 1:ml
            acc_up = 0.0
            acc_dn = 0.0
            for j = 1:ml, k = 1:ml
                v = VLL[i, j, k, l]
                acc_up += rhoup[k, j] * v
                acc_dn += rhodn[k, j] * v
            end
            Fup[i, l] -= acc_up
            Fdn[i, l] -= acc_dn
        end

        VRR = win.VRR
        for i = 1:mr, l = 1:mr
            acc_up = 0.0
            acc_dn = 0.0
            for j = 1:mr, k = 1:mr
                v = VRR[i, j, k, l]
                acc_up += rhoup[ml + lc + k, ml + lc + j] * v
                acc_dn += rhodn[ml + lc + k, ml + lc + j] * v
            end
            Fup[ml + lc + i, ml + lc + l] -= acc_up
            Fdn[ml + lc + i, ml + lc + l] -= acc_dn
        end
        if timing
            _bench_timing[:direct_LL_RR] += (time_ns() - t0) * 1e-9
        end

        if timing
            t0 = time_ns()
        end
        for np = 1:ns
            cols_p = center_cols[np]
            any(c -> c != 0, cols_p) || continue
            for nq = 1:ns
                cols_q = center_cols[nq]
                any(c -> c != 0, cols_q) || continue
                rho_nm_up = Matrix{Float64}(undef, length(cols_p), length(cols_q))
                rho_nm_dn = Matrix{Float64}(undef, length(cols_p), length(cols_q))
                mul!(rho_nm_up, Srho_up[np], S[nq]', 1.0, 0.0)
                mul!(rho_nm_dn, Srho_dn[np], S[nq]', 1.0, 0.0)
                Vpq = _slice_vee(V, np, nq)
                for a = 1:length(cols_p)
                    col_a = cols_p[a]
                    col_a == 0 && continue
                    for d = 1:length(cols_q)
                        col_d = cols_q[d]
                        col_d == 0 && continue
                        acc_up = 0.0
                        acc_dn = 0.0
                        for b = 1:length(cols_p), c = 1:length(cols_q)
                            v = Vpq[a, b, c, d]
                            acc_up += rho_nm_up[b, c] * v
                            acc_dn += rho_nm_dn[b, c] * v
                        end
                        Fup[col_a, col_d] -= acc_up
                        Fdn[col_a, col_d] -= acc_dn
                    end
                end
            end
        end

        WL = win.WL
        WLswap = win.WLswap
        for s = 1:ns
            cols = center_cols[s]
            for b = 1:length(cols)
                col_b = cols[b]
                col_b == 0 && continue
                for i = 1:ml
                    acc_up = 0.0
                    acc_dn = 0.0
                    for a = 1:length(cols)
                        col_a = cols[a]
                        col_a == 0 && continue
                        for j = 1:ml
                            v = WL[s][i, j, a, b]
                            acc_up += rhoup[j, col_a] * v
                            acc_dn += rhodn[j, col_a] * v
                        end
                    end
                    Fup[i, col_b] -= acc_up
                    Fdn[i, col_b] -= acc_dn
                end
            end
            for a = 1:length(cols)
                col_a = cols[a]
                col_a == 0 && continue
                for j = 1:ml
                    acc_up = 0.0
                    acc_dn = 0.0
                    for b = 1:length(cols)
                        col_b = cols[b]
                        col_b == 0 && continue
                        for i = 1:ml
                            v = WLswap[s][i, j, a, b]
                            acc_up += rhoup[col_b, i] * v
                            acc_dn += rhodn[col_b, i] * v
                        end
                    end
                    Fup[col_a, j] -= acc_up
                    Fdn[col_a, j] -= acc_dn
                end
            end
        end

        WR = win.WR
        WRswap = win.WRswap
        for s = 1:ns
            cols = center_cols[s]
            for b = 1:length(cols)
                col_b = cols[b]
                col_b == 0 && continue
                for i = 1:mr
                    acc_up = 0.0
                    acc_dn = 0.0
                    for a = 1:length(cols)
                        col_a = cols[a]
                        col_a == 0 && continue
                        for j = 1:mr
                            v = WR[s][i, j, a, b]
                            acc_up += rhoup[ml + lc + j, col_a] * v
                            acc_dn += rhodn[ml + lc + j, col_a] * v
                        end
                    end
                    Fup[ml + lc + i, col_b] -= acc_up
                    Fdn[ml + lc + i, col_b] -= acc_dn
                end
            end
            for a = 1:length(cols)
                col_a = cols[a]
                col_a == 0 && continue
                for j = 1:mr
                    acc_up = 0.0
                    acc_dn = 0.0
                    for b = 1:length(cols)
                        col_b = cols[b]
                        col_b == 0 && continue
                        for i = 1:mr
                            v = WRswap[s][i, j, a, b]
                            acc_up += rhoup[col_b, ml + lc + i] * v
                            acc_dn += rhodn[col_b, ml + lc + i] * v
                        end
                    end
                    Fup[col_a, ml + lc + j] -= acc_up
                    Fdn[col_a, ml + lc + j] -= acc_dn
                end
            end
        end
        if timing
            _bench_timing[:direct_CL_CR] += (time_ns() - t0) * 1e-9
        end

        if timing
            t0 = time_ns()
        end
        for i = 1:ml, l = 1:mr
            acc_up = 0.0
            acc_dn = 0.0
            for j = 1:ml, k = 1:mr
                v = VLR[i, j, k, l]
                acc_up += rhoup[j, ml + lc + k] * v
                acc_dn += rhodn[j, ml + lc + k] * v
            end
            Fup[i, ml + lc + l] -= acc_up
            Fdn[i, ml + lc + l] -= acc_dn
        end
        for k = 1:mr, j = 1:ml
            acc_up = 0.0
            acc_dn = 0.0
            for l = 1:mr, i = 1:ml
                v = VRL[k, l, i, j]
                acc_up += rhoup[ml + lc + l, i] * v
                acc_dn += rhodn[ml + lc + l, i] * v
            end
            Fup[ml + lc + k, j] -= acc_up
            Fdn[ml + lc + k, j] -= acc_dn
        end
        if timing
            _bench_timing[:direct_LR] += (time_ns() - t0) * 1e-9
        end
    else
        # Correctness-first path: when any slice is split, exchange falls back to
        # full slice-projection to maintain parity with the projection backend.
        # This is intentionally non-local.
        if timing
            t0 = time_ns()
        end
        superdim = size(Fup, 1)
        for n1 = 1:ns
            dn = dims[n1]
            for m = 1:ns
                dm = dims[m]
                Vnm = _slice_vee(V, n1, m)
                tmp = Matrix{Float64}(undef, dn, superdim)
                for (Fspin, Srho_spin) in ((Fup, Srho_up), (Fdn, Srho_dn))
                    rho_nm = Matrix{Float64}(undef, dn, dm)
                    mul!(rho_nm, Srho_spin[n1], S[m]', 1.0, 0.0)
                    Knm = Matrix{Float64}(undef, dn, dm)
                    for a = 1:dn, d = 1:dm
                        acc = 0.0
                        for b = 1:dn, c = 1:dm
                            acc += rho_nm[b, c] * Vnm[a, b, c, d]
                        end
                        Knm[a, d] = acc
                    end
                    mul!(tmp, Knm, S[m], 1.0, 0.0)
                    mul!(Fspin, S[n1]', tmp, -1.0, 1.0)
                end
            end
        end
        if timing
            _bench_timing[:direct_CL_CR] += (time_ns() - t0) * 1e-9
        end
    end

    if timing
        t0 = time_ns()
    end
    rho_slice = win.rho_slice
    @views for n1 = 1:ns
        Sn = S[n1]
        mul!(rho_slice[n1], Srho[n1], Sn', 1.0, 0.0)
    end
    # J[p,r] is block diagonal in the output slice pair (p,r).
    J = win.K
    tmp_n = win.tmp_n
    @views for n1 = 1:ns
        J1 = J[n1]
        Jvec = reshape(J1, :)
        fill!(Jvec, 0.0)
        dn1 = dims[n1]
        dn1sq = dn1 * dn1
        for n2 = 1:ns
            dn2 = dims[n2]
            Vnm = _slice_vee(V, n1, n2)
            Vpair = reshape(Vnm, dn1sq, dn2 * dn2)
            rvec = reshape(rho_slice[n2], dn2 * dn2)
            mul!(Jvec, Vpair, rvec, 1.0, 1.0)
        end
        tmp = tmp_n[n1]
        mul!(tmp, J1, S[n1], 1.0, 0.0)
        mul!(Fup, S[n1]', tmp, 1.0, 1.0)
        mul!(Fdn, S[n1]', tmp, 1.0, 1.0)
    end
    if timing
        _bench_timing[:exchange] += (time_ns() - t0) * 1e-9
    end

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
