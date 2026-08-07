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

_retire_without_rebuild(::SlicedBasisBackendCached) = true

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

function _exact_symmetric_matrix(A, name)
    size(A, 1) == size(A, 2) || error("$name must be square")
    scale = 1.0; disagreement = 0.0
    @inbounds for j = 1:size(A, 2), i = 1:j
        a, b = A[i, j], A[j, i]
        isfinite(a) && isfinite(b) ||
            error("$name must be finite in roundoff-exact sliced mode")
        scale = max(scale, abs(a), abs(b))
        disagreement = max(disagreement, abs(a - b))
    end
    disagreement <= 128eps(Float64) * scale ||
        error("$name must be symmetric in roundoff-exact sliced mode")
end
function _frozen_backend_check(Hup, Hdn, backend::SlicedBasisBackendCached,
        psiup, psidn, restricted, partition)
    layout, S = backend.layout, nslices(backend.layout)
    partition isa SliceLayout && partition.dims == layout.dims &&
        partition.offs == layout.offs ||
        error("roundoff-exact cached sliced mode requires block_partition == backend.layout")
    9maximum(layout.dims) + 2 <= 40 ||
        error("roundoff-exact cached sliced audit does not support this layout")
    _exact_symmetric_matrix(Hup, "Hup"); _exact_symmetric_matrix(Hdn, "Hdn")
    N = layout.offs[end]
    N*eps(Float64)<0.5 || error("roundoff-exact sliced dimension exceeds arithmetic band")
    for (C, name) in ((psiup, "psiup0"), (psidn, "psidn0"))
        size(C, 1) == N || error("$name has wrong row dimension")
        all(isfinite, C) || error("$name must be finite in roundoff-exact sliced mode")
        g = N * eps(Float64) / (1 - N * eps(Float64))
        norm(C' * C - I, Inf) <= 128g * max(1, size(C, 2)) ||
            error("$name must be orthonormal in roundoff-exact sliced mode")
    end
    vmax, verr = 0.0, 0.0
    for n = 1:S, m = 1:S
        Vnm, Vmn = _slice_vee(backend.V, n, m), _slice_vee(backend.V, m, n)
        all(isfinite, Vnm) || error("sliced interaction must be finite in roundoff-exact mode")
        dn, dm = layout.dims[n], layout.dims[m]
        @inbounds for a = 1:dn, b = 1:dn, c = 1:dm, d = 1:dm
            v = Vnm[a, b, c, d]; vmax = max(vmax, abs(v))
            verr = max(verr, abs(v - Vnm[b, a, c, d]),
                abs(v - Vnm[a, b, d, c]), abs(v - Vmn[c, d, a, b]))
        end
    end
    verr <= 128eps(Float64) * max(1.0, vmax) ||
        error("sliced interaction lacks required physical symmetry")
    Float64
end

_empty_frozen_field(backend::SlicedBasisBackendCached, m) = _SlicedFrozenField(
    zeros(sum(abs2, backend.layout.dims)), zeros(m, m), zeros(m, m))
_frozen_field_empty(f::_SlicedFrozenField) = all(iszero, f.Jpack) &&
    all(iszero, f.kup) && all(iszero, f.kdn)

mutable struct _SlicedFrozenScratch
    pair_offs::Vector{Int}
    sources::Vector{Int}
    direct::Vector{Float64}
    source_density::Vector{Float64}
    density1::Matrix{Float64}
    density2::Matrix{Float64}
    kernel::Matrix{Float64}
    work::Matrix{Float64}
    square::Matrix{Float64}
    Gup::Matrix{Float64}
    Gdn::Matrix{Float64}
end
function _SlicedFrozenScratch(backend::SlicedBasisBackendCached)
    layout = backend.layout
    pair_offs = _pair_offs(layout, 1:nslices(layout))
    d = maximum(layout.dims)
    _SlicedFrozenScratch(pair_offs, zeros(Int, nslices(layout)),
        zeros(pair_offs[end]), zeros(pair_offs[end]), zeros(d, d), zeros(d, d),
        zeros(d, d), zeros(0, d), zeros(0, 0), zeros(0, 0), zeros(0, 0))
end
function _sliced_frozen_scratch!(p, backend::SlicedBasisBackendCached, n)
    p.scratch isa _SlicedFrozenScratch ||
        (p.scratch = _SlicedFrozenScratch(backend))
    scratch = p.scratch::_SlicedFrozenScratch
    d = maximum(backend.layout.dims)
    if size(scratch.square, 1) < n
        scratch.work = zeros(n, d)
        scratch.square = zeros(n, n)
        scratch.Gup = zeros(n, n)
        scratch.Gdn = zeros(n, n)
    end
    scratch
end

function _sliced_frozen_direct!(scratch, backend, Cup, Cdn)
    dims, offs, S = backend.layout.dims, backend.layout.offs,
        nslices(backend.layout)
    fill!(scratch.direct, 0); nsources = 0
    for m = 1:S
        rm = offs[m] + 1:offs[m + 1]
        cup, cdn = @view(Cup[rm, :]), @view(Cdn[rm, :])
        (any(x -> !iszero(x), cup) || any(x -> !iszero(x), cdn)) || continue
        nsources += 1; scratch.sources[nsources] = m
        T = reshape(@view(scratch.source_density[
            _pair_range(scratch.pair_offs, m)]), dims[m], dims[m])
        extra = @view scratch.density2[1:dims[m], 1:dims[m]]
        mul!(T, cup, transpose(cup))
        mul!(extra, cdn, transpose(cdn)); T .+= extra
    end
    for n = 1:S, source = 1:nsources
        m = scratch.sources[source]
        Jn = reshape(@view(scratch.direct[
            _pair_range(scratch.pair_offs, n)]), dims[n], dims[n])
        T = reshape(@view(scratch.source_density[
            _pair_range(scratch.pair_offs, m)]), dims[m], dims[m])
        Vnm = _slice_vee(backend.V, n, m)
        @inbounds for a = 1:dims[n], b = 1:dims[n], c = 1:dims[m], d = 1:dims[m]
            Jn[a, b] += Vnm[a, b, c, d] * T[c, d]
        end
    end
    scratch.direct
end
function _sliced_frozen_exchange_projected!(Kbar, scratch, backend, C, ra, phi)
    fill!(Kbar, 0); size(C, 2) == 0 && return Kbar
    layout, dims, offs = backend.layout, backend.layout.dims, backend.layout.offs
    k = size(phi, 2)
    for n in _active_slices(layout, ra), m in _active_slices(layout, ra)
        dn, dm = dims[n], dims[m]
        rn, rm = offs[n] + 1:offs[n + 1], offs[m] + 1:offs[m + 1]
        D = @view scratch.density1[1:dn, 1:dm]
        K = @view scratch.kernel[1:dn, 1:dm]
        mul!(D, @view(C[rn, :]), transpose(@view(C[rm, :])))
        fill!(K, 0); Vnm = _slice_vee(backend.V, n, m)
        @inbounds for a = 1:dn, d = 1:dm, b = 1:dn, c = 1:dm
            K[a, d] += Vnm[a, b, c, d] * D[b, c]
        end
        phin, phim = _slice_phi(layout, ra, phi, n),
            _slice_phi(layout, ra, phi, m)
        work = @view scratch.work[1:k, 1:dm]
        projected = @view scratch.square[1:k, 1:k]
        mul!(work, transpose(phin), K); mul!(projected, work, phim)
        Kbar .+= projected
    end
    Kbar
end
function _sliced_frozen_cross_exchange!(scratch, backend, Cf, Cp, ra)
    (size(Cf, 2) == 0 || size(Cp, 2) == 0) && return 0.0
    _sliced_frozen_cross_exchange(backend, Cf, Cp, ra)
end
function _sliced_frozen_direct_trace!(scratch, backend, Cup, Cdn, J)
    offs, dims = backend.layout.offs, backend.layout.dims; value = 0.0
    all(==(size(scratch.density1, 1)), dims) ||
        return _sliced_frozen_direct_trace(backend, Cup, Cdn, J)
    for n = 1:nslices(backend.layout)
        d = dims[n]; rn = offs[n] + 1:offs[n + 1]
        T = @view scratch.density1[1:d, 1:d]
        extra = @view scratch.density2[1:d, 1:d]
        mul!(T, @view(Cup[rn, :]), transpose(@view(Cup[rn, :])))
        mul!(extra, @view(Cdn[rn, :]), transpose(@view(Cdn[rn, :])))
        T .+= extra
        value += dot(T, reshape(@view(J[
            _pair_range(scratch.pair_offs, n)]), d, d))
    end
    value
end

function _sliced_frozen_direct(backend, Cup, Cdn)
    dims, offs, S = backend.layout.dims, backend.layout.offs, nslices(backend.layout)
    poffs = _pair_offs(backend.layout, 1:S); J = zeros(poffs[end])
    sources = [m for m=1:S if any(x -> !iszero(x), @view(Cup[offs[m]+1:offs[m+1],:])) ||
        any(x -> !iszero(x), @view(Cdn[offs[m]+1:offs[m+1],:]))]
    T = [(@view Cup[offs[m]+1:offs[m+1],:]) * (@view Cup[offs[m]+1:offs[m+1],:])' +
        (@view Cdn[offs[m]+1:offs[m+1],:]) * (@view Cdn[offs[m]+1:offs[m+1],:])'
        for m in sources]
    for n = 1:S, (i, m) in enumerate(sources)
        Jn = reshape(@view(J[_pair_range(poffs, n)]), dims[n], dims[n])
        Vnm = _slice_vee(backend.V, n, m)
        @inbounds for a = 1:dims[n], b = 1:dims[n], c = 1:dims[m], d = 1:dims[m]
            Jn[a, b] += Vnm[a, b, c, d] * T[i][c, d]
        end
    end
    J
end
function _sliced_frozen_exchange_projected(backend, C, ra, phi)
    layout, dims, offs = backend.layout, backend.layout.dims, backend.layout.offs
    active = _active_slices(layout, ra); Kbar = zeros(size(phi,2),size(phi,2))
    for n in active, m in active
        rn, rm = offs[n] + 1:offs[n + 1], offs[m] + 1:offs[m + 1]
        D = (@view C[rn, :]) * (@view C[rm, :])'; K = zeros(dims[n], dims[m])
        Vnm = _slice_vee(backend.V, n, m)
        @inbounds for a = 1:dims[n], d = 1:dims[m], b = 1:dims[n], c = 1:dims[m]
            K[a, d] += Vnm[a, b, c, d] * D[b, c]
        end
        phin, phim = _slice_phi(layout, ra, phi, n), _slice_phi(layout, ra, phi, m)
        Kbar .+= phin' * K * phim
    end
    Kbar
end
function _sliced_frozen_cross_exchange(backend, Cf, Cp, ra)
    layout, dims, offs = backend.layout, backend.layout.dims, backend.layout.offs; value = 0.0
    for n in _active_slices(layout,ra), m in _active_slices(layout,ra)
        rn,rm=offs[n]+1:offs[n+1],offs[m]+1:offs[m+1]
        Df=(@view Cf[rn,:])*(@view Cf[rm,:])'; Dp=(@view Cp[rn,:])*(@view Cp[rm,:])'
        K=zeros(dims[n],dims[m]); Vnm=_slice_vee(backend.V,n,m)
        @inbounds for a=1:dims[n],d=1:dims[m],b=1:dims[n],c=1:dims[m]
            K[a,d]+=Vnm[a,b,c,d]*Dp[b,c]
        end
        value+=dot(Df,K)
    end
    value
end
function _sliced_frozen_direct_trace(backend, Cup, Cdn, J)
    offs, dims = backend.layout.offs, backend.layout.dims
    poffs = _pair_offs(backend.layout, 1:nslices(backend.layout)); value = 0.0
    for n = 1:nslices(backend.layout)
        rn = offs[n] + 1:offs[n + 1]
        Tn = (@view Cup[rn, :]) * (@view Cup[rn, :])' + (@view Cdn[rn, :]) * (@view Cdn[rn, :])'
        value += dot(Tn, reshape(@view(J[_pair_range(poffs, n)]), dims[n], dims[n]))
    end
    value
end
function _sliced_frozen_fresh(backend, p, Cup, Cdn, phi, ra)
    J = _sliced_frozen_direct(backend, Cup, Cdn)
    Kup, Kdn = _sliced_frozen_exchange_projected(backend, Cup, ra, phi),
        _sliced_frozen_exchange_projected(backend, Cdn, ra, phi)
    one = dot(Cup, p.Hup * Cup) + dot(Cdn, p.Hdn * Cdn); e = one + 0.5 *
        (_sliced_frozen_direct_trace(backend, Cup, Cdn, J) - _sliced_frozen_cross_exchange(backend, Cup, Cup, ra) -
        _sliced_frozen_cross_exchange(backend, Cdn, Cdn, ra))
    _SlicedFrozenField(J, Kup, Kdn), e
end
function _frozen_absorption_map(backend::SlicedBasisBackendCached, side, oldblock,
        cra, newblock)
    _aligned_absorption(side, oldblock.ra, cra, newblock.raV, backend.layout) ||
        error("roundoff-exact cached absorption is not whole-slice aligned")
    A = _old_basis_map(side, oldblock.vee, cra, newblock.phi)
    A === nothing && error("roundoff-exact cached absorption would use fallback")
    A
end
function _promote_frozen_field(backend::SlicedBasisBackendCached, p, old, oldblock,
        newblock, growth, A, diagnostics)
    ra, f = newblock.ra, old.field
    if size(growth.Cup,2)+size(growth.Cdn,2)==0
        return _SlicedFrozenField(f.Jpack,A'*f.kup*A,A'*f.kdn*A),old.energy
    end
    scratch = _sliced_frozen_scratch!(p, backend, size(newblock.phi, 2))
    JP = _sliced_frozen_direct!(scratch, backend, growth.Cup, growth.Cdn)
    KPup, KPdn = zeros(size(newblock.phi, 2), size(newblock.phi, 2)),
        zeros(size(newblock.phi, 2), size(newblock.phi, 2))
    _sliced_frozen_exchange_projected!(KPup, scratch, backend, growth.Cup, ra,
        newblock.phi)
    _sliced_frozen_exchange_projected!(KPdn, scratch, backend, growth.Cdn, ra,
        newblock.phi)
    Cfu, Cfd = _ledger_columns(old, _frozen_size(p), :up),
        _ledger_columns(old, _frozen_size(p), :dn)
    fdirect = _sliced_frozen_direct_trace!(scratch, backend, growth.Cup,
        growth.Cdn, f.Jpack)
    Xup, Xdn = oldblock.phi' * growth.Cup[oldblock.ra, :],
        oldblock.phi' * growth.Cdn[oldblock.ra, :]
    fexchange = dot(Xup, f.kup * Xup) + dot(Xdn, f.kdn * Xdn)
    rdirect = _sliced_frozen_direct_trace!(scratch, backend, Cfu, Cfd, JP)
    rexchange = _sliced_frozen_cross_exchange!(
        scratch, backend, Cfu, growth.Cup, ra) +
        _sliced_frozen_cross_exchange!(scratch, backend, Cfd, growth.Cdn, ra)
    self = 0.5 * (_sliced_frozen_direct_trace!(scratch, backend, growth.Cup,
        growth.Cdn, JP) - _sliced_frozen_cross_exchange!(
            scratch, backend, growth.Cup, growth.Cup, ra) -
        _sliced_frozen_cross_exchange!(scratch, backend, growth.Cdn,
            growth.Cdn, ra))
    one = dot(growth.Cup, p.Hup * growth.Cup) +
        dot(growth.Cdn, p.Hdn * growth.Cdn)
    kup = A' * f.kup * A + KPup
    kdn = A' * f.kdn * A + KPdn
    if diagnostics !== nothing
        diagnostics.frozen_promotions_up += size(growth.Cup, 2)
        diagnostics.frozen_promotions_dn += size(growth.Cdn, 2)
        diagnostics.frozen_direct_orientation_error = max(
            diagnostics.frozen_direct_orientation_error, abs(fdirect - rdirect))
        diagnostics.frozen_exchange_orientation_error = max(
            diagnostics.frozen_exchange_orientation_error, abs(fexchange - rexchange))
    end
    _SlicedFrozenField(f.Jpack + JP, kup, kdn),
        old.energy + one + self + 0.5 * (fdirect + rdirect - fexchange - rexchange)
end
function _sliced_project_frozen_direct(backend, B, J)
    layout = backend.layout; poffs = _pair_offs(layout, 1:nslices(layout)); G = zeros(size(B, 2), size(B, 2))
    for n = 1:nslices(layout)
        rn = layout.offs[n] + 1:layout.offs[n + 1]
        Jn = reshape(@view(J[_pair_range(poffs, n)]), layout.dims[n], layout.dims[n]); Bn = @view B[rn, :]; G .+= Bn' * Jn * Bn
    end
    G
end
function _sliced_project_frozen_direct!(G, scratch, backend, B, J)
    layout = backend.layout; fill!(G, 0); w = size(B, 2)
    for n = 1:nslices(layout)
        d = layout.dims[n]; rn = layout.offs[n] + 1:layout.offs[n + 1]
        Jn = reshape(@view(J[_pair_range(scratch.pair_offs, n)]), d, d)
        Bn = @view B[rn, :]
        work = @view scratch.work[1:w, 1:d]
        projected = @view scratch.square[1:w, 1:w]
        mul!(work, transpose(Bn), Jn); mul!(projected, work, Bn)
        G .+= projected
    end
    G
end
function _frozen_effective(backend::SlicedBasisBackendCached, p, L, R, B)
    N, lf, rf = size(B, 1), L.field, R.field
    scratch = _sliced_frozen_scratch!(p, backend, size(B, 2))
    CLu, CLd = _ledger_columns(L, N, :up), _ledger_columns(L, N, :dn)
    CRu, CRd = _ledger_columns(R, N, :up), _ledger_columns(R, N, :dn)
    @. scratch.direct = lf.Jpack + rf.Jpack
    w = size(B, 2); Gup = @view scratch.Gup[1:w, 1:w]
    Gdn = @view scratch.Gdn[1:w, 1:w]
    _sliced_project_frozen_direct!(Gup, scratch, backend, B, scratch.direct)
    copyto!(Gdn, Gup); ml, mr = size(lf.kup, 1), size(rf.kup, 1)
    Gup[1:ml, 1:ml] -= lf.kup; Gdn[1:ml, 1:ml] -= lf.kdn
    Gup[end - mr + 1:end, end - mr + 1:end] -= rf.kup; Gdn[end - mr + 1:end, end - mr + 1:end] -= rf.kdn
    lr = _sliced_frozen_direct_trace!(scratch, backend, CLu, CLd, rf.Jpack)
    rl = _sliced_frozen_direct_trace!(scratch, backend, CRu, CRd, lf.Jpack)
    Gup, Gdn, L.energy + R.energy + 0.5 * (lr + rl),
        hcat(CLu, CRu), hcat(CLd, CRd)
end
_frozen_window_check(backend::SlicedBasisBackendCached, win) =
    win isa SlicedBasisCachedWindow ||
        error("roundoff-exact cached sliced mode forbids projection windows")
_frozen_audit_context(p, spin, backend::SlicedBasisBackendCached) = spin === :up ?
    _SlicedFockAuditContext(p.Hup, backend, p.psiup, p.psidn, :up) :
    _SlicedFockAuditContext(p.Hdn, backend, p.psiup, p.psidn, :dn)
function _sliced_frozen_payload(p)
    values, seen = 0, IdDict{Any,Nothing}()
    for i = eachindex(p.cuts)
        isassigned(p.cuts, i) || continue
        cut = p.cuts[i]; f = cut.field
        for array in (f.Jpack, f.kup, f.kdn)
            haskey(seen, array) || (seen[array] = nothing; values += length(array))
        end
        values += 3
        node = cut.ledger
        while node !== nothing
            if !haskey(seen, node)
                seen[node] = nothing
                for array in (node.up, node.dn, node.originup, node.origindn)
                    haskey(seen, array) ||
                        (seen[array] = nothing; values += length(array))
                end
            end
            node = node.parent
        end
    end
    values
end
_frozen_payload(p::_FrozenPolicy{T,B}) where {T,B<:SlicedBasisBackendCached} =
    _sliced_frozen_payload(p)

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
