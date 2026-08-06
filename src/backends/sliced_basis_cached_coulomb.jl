struct _SlicedBasisBackendCachedCoulomb{T}
    layout::SliceLayout
    V::T
end
struct _SlicedCoulombBlockState
    ra::UnitRange{Int}
    phi::Matrix{Float64}
    G::Matrix{Float64}
    exterior_slices::UnitRange{Int}
    pair_offs::Vector{Int}
    W::Matrix{Float64}
end
struct _SlicedCoulombScratch
    density::Matrix{Float64}
    physical_work::Matrix{Float64}
    mixed::Matrix{Float64}
    exchange::Matrix{Float64}
    mixed_dn::Matrix{Float64}
    exchange_dn::Matrix{Float64}
    pair_density::Vector{Float64}
    pair_field::Vector{Float64}
    exterior_pair::Vector{Float64}
end
struct _SlicedCoulombWindow{T}
    V::T
    layout::SliceLayout
    ml::Int
    lc::Int
    mr::Int
    center_slices::UnitRange{Int}
    L::_SlicedCoulombBlockState
    R::_SlicedCoulombBlockState
    scratch::_SlicedCoulombScratch
end
@inline _npair(n) = n * (n + 1) ÷ 2
@inline _sym_pair(i, j) = (a = min(i, j); b = max(i, j); b * (b - 1) ÷ 2 + a)
@inline _sym_scale(i, j) = i == j ? 1.0 : sqrt(2.0)
_sym_pair_offs(layout, slices) = [0; cumsum(_npair.(layout.dims[slices]))]
_sym_state_range(state, s) = _pair_range(
    state.pair_offs, s - first(state.exterior_slices) + 1)
function _certify_coulomb(layout, V)
    dims = layout.dims
    for n in eachindex(dims), m in eachindex(dims)
        X, Y = _slice_vee(V, n, m), _slice_vee(V, m, n)
        dn, dm = dims[n], dims[m]
        scale = max(1.0, maximum(abs, X), maximum(abs, Y))
        tol = 256 * eps(Float64) * max(dn, dm)^2 * scale
        err1 = err2 = errx = 0.0
        for d = 1:dm, c = 1:dm, b = 1:dn, a = 1:dn
            value = X[a, b, c, d]
            err1 = max(err1, abs(value - X[b, a, c, d]))
            err2 = max(err2, abs(value - X[a, b, d, c]))
            errx = max(errx, abs(value - Y[c, d, a, b]))
        end
        err1 <= tol || error("sliced interaction fails first-pair Coulomb symmetry")
        err2 <= tol || error("sliced interaction fails second-pair Coulomb symmetry")
        errx <= tol || error("sliced interaction fails pair-exchange Coulomb symmetry")
    end
end
function _SlicedBasisBackendCachedCoulomb(layout::SliceLayout, V6::Array{Float64,6})
    generic = SlicedBasisBackendCached(layout, V6)
    _certify_coulomb(layout, V6)
    _SlicedBasisBackendCachedCoulomb{typeof(generic.V)}(layout, generic.V)
end
function _SlicedBasisBackendCachedCoulomb(layout::SliceLayout, vee::SlicedVeeRagged)
    layout.dims == vee.layout.dims || error("layout does not match ragged interaction")
    _certify_coulomb(layout, vee)
    _SlicedBasisBackendCachedCoulomb{SlicedVeeRagged}(layout, vee)
end
function _pack_sym!(v, X, add = false)
    n = size(X, 1)
    for j = 1:n, i = 1:j
        index = _sym_pair(i, j)
        value = _sym_scale(i, j) * X[i, j]
        v[index] = add ? v[index] + value : value
    end
end
function _unpack_sym!(X, v)
    n = size(X, 1)
    for j = 1:n, i = 1:j
        value = v[_sym_pair(i, j)] / _sym_scale(i, j)
        X[i, j] = X[j, i] = value
    end
end
function _pack_coulomb!(Xhat, X)
    p, q = size(X, 1), size(X, 3)
    for l = 1:q, j = 1:l, k = 1:p, i = 1:k
        Xhat[_sym_pair(i, k), _sym_pair(j, l)] =
            _sym_scale(i, k) * _sym_scale(j, l) * X[i, k, j, l]
    end
end
@inline function _sym_value(X, i, k, j, l)
    X[_sym_pair(i, k), _sym_pair(j, l)] /
        (_sym_scale(i, k) * _sym_scale(j, l))
end
function _add_sym_transform!(dest, X, L, R, stage, input, output, work;
        add_transpose = false)
    p, rl = size(L); q, rr = size(R)
    Kq, Kl, Kr = _npair(q), _npair(rl), _npair(rr)
    S = view(stage, 1:Kl, 1:Kq)
    Xin = view(input, 1:p, 1:p)
    Yout = view(output, 1:rl, 1:rl)
    for col = 1:Kq
        _unpack_sym!(Xin, view(X, :, col))
        _two_sided!(Yout, Xin, L, work)
        _pack_sym!(view(S, :, col), Yout)
    end
    Xin = view(input, 1:q, 1:q)
    Yout = view(output, 1:rr, 1:rr)
    for row = 1:Kl
        _unpack_sym!(Xin, view(S, row, :))
        _two_sided!(Yout, Xin, R, work)
        for l = 1:rr, j = 1:l
            col = _sym_pair(j, l)
            value = _sym_scale(j, l) * Yout[j, l]
            dest[row, col] += value
            add_transpose && (dest[col, row] += value)
        end
    end
end
function _transform_sym_channels!(dest, source, A, source_cols, input, output, work)
    p, r = size(A)
    Xin, Yout = view(input, 1:p, 1:p), view(output, 1:r, 1:r)
    for (col, oldcol) in enumerate(source_cols)
        _unpack_sym!(Xin, view(source, :, oldcol))
        _two_sided!(Yout, Xin, A, work)
        _pack_sym!(view(dest, :, col), Yout)
    end
end
function _add_sym_left!(dest, X, A, input, output, work)
    p, r = size(A)
    Xin, Yout = view(input, 1:p, 1:p), view(output, 1:r, 1:r)
    for col in axes(X, 2)
        _unpack_sym!(Xin, view(X, :, col))
        _two_sided!(Yout, Xin, A, work)
        _pack_sym!(view(dest, :, col), Yout, true)
    end
end
_sym_work(maxin, maxout) = (zeros(_npair(maxout), _npair(maxin)), zeros(maxin, maxin), zeros(maxout, maxout), zeros(maxin, maxout))
function _build_coulomb_state(ra, raV, phi_in, backend)
    phi = Matrix{Float64}(phi_in); layout, dims = backend.layout, backend.layout.dims
    m = size(phi, 2); P = _build_slice_basis(layout, ra, phi)
    active, exterior = _active_slices(layout, ra), _active_slices(layout, raV)
    offs = _sym_pair_offs(layout, exterior); W = zeros(_npair(m), offs[end])
    maxd = maximum(dims); stage, input, output, work = _sym_work(maxd, m)
    raw = zeros(_npair(maxd), _npair(maxd))
    for n in active, (js, s) in enumerate(exterior)
        dn, ds = dims[n], dims[s]
        X = view(raw, 1:_npair(dn), 1:_npair(ds))
        _pack_coulomb!(X, _slice_vee(backend.V, n, s))
        _add_sym_left!(view(W, :, _pair_range(offs, js)), X, P[n], input, output, work)
    end
    G = zeros(_npair(m), _npair(m))
    for n in active, s in active
        dn, ds = dims[n], dims[s]
        X = view(raw, 1:_npair(dn), 1:_npair(ds))
        _pack_coulomb!(X, _slice_vee(backend.V, n, s))
        _add_sym_transform!(G, X, P[n], P[s], stage, input, output, work)
    end
    _SlicedCoulombBlockState(ra, phi, G, exterior, offs, W)
end
function vee_init_block(side, ra, raV, phi, backend::_SlicedBasisBackendCachedCoulomb)
    _check_side(side); _check_range(ra, "ra"); _check_range(raV, "raV")
    _check_cached_exterior(side, ra, raV, backend.layout.offs[end])
    size(phi, 1) == length(ra) || error("phi has wrong row count for ra")
    _bench_timing_enabled() && (_bench_timing[:cache_init_calls] += 1)
    _build_coulomb_state(ra, raV, phi, backend)
end
function _old_basis_map(side, old::_SlicedCoulombBlockState, cra, phi_new)
    rows = side === :left ? (1:length(old.ra)) : (length(cra) + 1:size(phi_new, 1))
    old_rows = view(phi_new, rows, :); A = old.phi' * old_rows
    err = norm(old_rows - old.phi * A)
    scale = max(1.0, norm(old_rows)) * max(size(old_rows)...)
    err <= 32 * eps(Float64) * scale ? A : nothing
end
function _incremental_coulomb_state(side, old, cra, phi_in, A, raV, backend)
    phi = Matrix{Float64}(phi_in); layout, dims = backend.layout, backend.layout.dims
    newra = side === :left ? (first(old.ra):last(cra)) : (first(cra):last(old.ra))
    absorbed, exterior = _active_slices(layout, cra), _active_slices(layout, raV)
    offs = _sym_pair_offs(layout, exterior); mnew, mold = size(phi, 2), size(old.phi, 2)
    W = zeros(_npair(mnew), offs[end]); maxd = maximum(dims)
    maxin = max(mold, maxd); stage, input, output, work = _sym_work(maxin, mnew)
    oldcols = (first(_sym_state_range(old, first(exterior))):
        last(_sym_state_range(old, last(exterior))))
    _transform_sym_channels!(W, old.W, A, oldcols, input, output, work)
    raw = zeros(_npair(maxd), _npair(maxd))
    for n in absorbed, (js, s) in enumerate(exterior)
        Cn = _slice_phi(layout, newra, phi, n); dn, ds = dims[n], dims[s]
        X = view(raw, 1:_npair(dn), 1:_npair(ds))
        _pack_coulomb!(X, _slice_vee(backend.V, n, s))
        _add_sym_left!(view(W, :, _pair_range(offs, js)), X, Cn, input, output, work)
    end
    G = zeros(_npair(mnew), _npair(mnew))
    _add_sym_transform!(G, old.G, A, A, stage, input, output, work)
    for n in absorbed
        Cn = _slice_phi(layout, newra, phi, n)
        _add_sym_transform!(G, view(old.W, :, _sym_state_range(old, n)), A, Cn,
            stage, input, output, work; add_transpose = true)
        for s in absorbed
            Cs = _slice_phi(layout, newra, phi, s); dn, ds = dims[n], dims[s]
            X = view(raw, 1:_npair(dn), 1:_npair(ds))
            _pack_coulomb!(X, _slice_vee(backend.V, n, s))
            _add_sym_transform!(G, X, Cn, Cs, stage, input, output, work)
        end
    end
    _SlicedCoulombBlockState(newra, phi, G, exterior, offs, W)
end
function vee_absorb_block(side, old::_SlicedCoulombBlockState, cra, Phi_old, Phi_C,
        phi_new, raV, backend::_SlicedBasisBackendCachedCoulomb)
    _check_side(side); _check_range(cra, "cra"); _check_range(raV, "raV")
    newra = side === :left ? (first(old.ra):last(cra)) : (first(cra):last(old.ra))
    _check_cached_exterior(side, newra, raV, backend.layout.offs[end])
    A = _aligned_absorption(side, old.ra, cra, raV, backend.layout) ?
        _old_basis_map(side, old, cra, phi_new) : nothing
    if A === nothing
        _bench_timing_enabled() && (_bench_timing[:cache_fallback_calls] += 1)
        return _build_coulomb_state(newra, raV, phi_new, backend)
    end
    _bench_timing_enabled() && (_bench_timing[:cache_incremental_calls] += 1)
    _incremental_coulomb_state(side, old, cra, phi_new, A, raV, backend)
end
function vee_window(L::_SlicedCoulombBlockState, R::_SlicedCoulombBlockState, Cra,
        backend::_SlicedBasisBackendCachedCoulomb)
    _check_range(Cra, "Cra"); N = backend.layout.offs[end]
    aligned = first(L.ra) == 1 && last(L.ra) + 1 == first(Cra) &&
        last(Cra) + 1 == first(R.ra) && last(R.ra) == N &&
        _whole_slice_range(backend.layout, L.ra) &&
        _whole_slice_range(backend.layout, Cra) && _whole_slice_range(backend.layout, R.ra)
    if !aligned
        _bench_timing_enabled() && (_bench_timing[:window_projection_calls] += 1)
        projection = SlicedBasisBackend(backend.layout, backend.V)
        return vee_window(SlicedBlockState(L.ra, L.phi), SlicedBlockState(R.ra, R.phi),
            Cra, projection)
    end
    _bench_timing_enabled() && (_bench_timing[:window_local_calls] += 1)
    ml, mr, maxd = size(L.phi, 2), size(R.phi, 2), maximum(backend.layout.dims)
    maxk = max(ml, mr)
    exterior_pairs = ml < 12 ? 0 : sum(_npair(backend.layout.dims[s])
        for s in _active_slices(backend.layout, R.ra))
    scratch = _SlicedCoulombScratch(zeros(maxd, maxd),
        zeros(maxd, maxk), zeros(maxk, maxd), zeros(maxk, maxd),
        zeros(maxk, maxd), zeros(maxk, maxd), zeros(maxk^2), zeros(maxk^2),
        zeros(exterior_pairs))
    centers = (first(_slice_local(backend.layout, first(Cra))):
        first(_slice_local(backend.layout, last(Cra))))
    _SlicedCoulombWindow(backend.V, backend.layout, ml, length(Cra), mr, centers,
        L, R, scratch)
end
@inline _sector_value(X, x, z, y, w, nx, ny, ::Val{:symmetric}) =
    _sym_value(X, x, z, y, w)
@inline _sector_value(X, x, z, y, w, nx, ny, ::Val{:symmetric_swapped}) =
    _sym_value(X, y, w, x, z)
function _add_sym_sector_both!(Fup, Fdn, rhoup, rhodn, X, xoff, yoff, nx, ny,
        restricted)
    _add_cached_sector!(Fup, Fdn, rhoup, rhodn, X, xoff, yoff, nx, ny,
        Val(:symmetric), restricted)
    _add_cached_sector!(Fup, Fdn, rhoup, rhodn, X, yoff, xoff, ny, nx,
        Val(:symmetric_swapped), restricted)
end
function _pack_total_density!(v, up, dn, off, n, ::Val{restricted}) where restricted
    for j = 1:n, i = 1:j
        value = restricted ? 2up[off + i, off + j] :
            up[off + i, off + j] + dn[off + i, off + j]
        v[_sym_pair(i, j)] = _sym_scale(i, j) * value
    end
end
function _add_coulomb_cross_direct_ordered!(Fup, Fdn, rhoup, rhodn, X, Y, xoff,
        yoff, layout, scratch, ::Val{restricted}) where restricted
    kx, ky = size(X.phi, 2), size(Y.phi, 2); Kx = _npair(kx)
    _pack_total_density!(scratch.pair_density, rhoup, rhodn, xoff, kx,
        Val(restricted)); fill!(view(scratch.pair_field, 1:Kx), 0)
    yrange = yoff + 1:yoff + ky
    for s in _active_slices(layout, Y.ra)
        P = _slice_phi(layout, Y.ra, Y.phi, s); d = size(P, 1)
        density = view(scratch.density, 1:d, 1:d)
        work = view(scratch.physical_work, 1:d, 1:ky)
        mul!(work, P, view(rhoup, yrange, yrange)); mul!(density, work, transpose(P))
        if restricted
            density .*= 2
        else
            mul!(work, P, view(rhodn, yrange, yrange))
            mul!(density, work, transpose(P), 1.0, 1.0)
        end
        W = view(X.W, :, _sym_state_range(X, s))
        for b = 1:d, a = 1:b, I = 1:Kx
            A = _sym_pair(a, b)
            scratch.pair_field[I] += W[I, A] * _sym_scale(a, b) * density[a, b]
        end
        for b = 1:d, a = 1:b
            A = _sym_pair(a, b); value = 0.0
            for I = 1:Kx
                value = muladd(W[I, A], scratch.pair_density[I], value)
            end
            value /= _sym_scale(a, b)
            density[a, b] = density[b, a] = value
        end
        mul!(work, density, P)
        mul!(view(Fup, yrange, yrange), transpose(P), work, 1.0, 1.0)
        restricted || mul!(view(Fdn, yrange, yrange), transpose(P), work, 1.0, 1.0)
    end
    for k = 1:kx, i = 1:k
        I = _sym_pair(i, k); value = scratch.pair_field[I] / _sym_scale(i, k)
        Fup[xoff + i, xoff + k] += value
        i != k && (Fup[xoff + k, xoff + i] += value)
        if !restricted
            Fdn[xoff + i, xoff + k] += value
            i != k && (Fdn[xoff + k, xoff + i] += value)
        end
    end
end
function _add_coulomb_cross_direct_mul!(Fup, Fdn, rhoup, rhodn, X, Y, xoff, yoff,
        layout, scratch, ::Val{restricted}) where restricted
    kx, ky = size(X.phi, 2), size(Y.phi, 2); Kx = _npair(kx)
    _pack_total_density!(scratch.pair_density, rhoup, rhodn, xoff, kx,
        Val(restricted))
    yrange = yoff + 1:yoff + ky
    slices = _active_slices(layout, Y.ra); qoff = 0
    for s in slices
        P = _slice_phi(layout, Y.ra, Y.phi, s); d = size(P, 1)
        density = view(scratch.density, 1:d, 1:d)
        work = view(scratch.physical_work, 1:d, 1:ky)
        mul!(work, P, view(rhoup, yrange, yrange)); mul!(density, work, transpose(P))
        if restricted
            density .*= 2
        else
            mul!(work, P, view(rhodn, yrange, yrange))
            mul!(density, work, transpose(P), 1.0, 1.0)
        end
        qrange = qoff + 1:qoff + _npair(d)
        _pack_sym!(view(scratch.exterior_pair, qrange), density); qoff = last(qrange)
    end
    cols = first(_sym_state_range(X, first(slices))):last(
        _sym_state_range(X, last(slices)))
    W = view(X.W, :, cols); q = view(scratch.exterior_pair, 1:qoff)
    mul!(view(scratch.pair_field, 1:Kx), W, q)
    mul!(q, transpose(W), view(scratch.pair_density, 1:Kx))
    qoff = 0
    for s in slices
        P = _slice_phi(layout, Y.ra, Y.phi, s); d = size(P, 1)
        density = view(scratch.density, 1:d, 1:d)
        work = view(scratch.physical_work, 1:d, 1:ky)
        qrange = qoff + 1:qoff + _npair(d)
        _unpack_sym!(density, view(q, qrange)); qoff = last(qrange)
        mul!(work, density, P)
        mul!(view(Fup, yrange, yrange), transpose(P), work, 1.0, 1.0)
        restricted || mul!(view(Fdn, yrange, yrange), transpose(P), work, 1.0, 1.0)
    end
    for k = 1:kx, i = 1:k
        I = _sym_pair(i, k); value = scratch.pair_field[I] / _sym_scale(i, k)
        Fup[xoff + i, xoff + k] += value
        i != k && (Fup[xoff + k, xoff + i] += value)
        if !restricted
            Fdn[xoff + i, xoff + k] += value
            i != k && (Fdn[xoff + k, xoff + i] += value)
        end
    end
end
function _add_coulomb_cross_direct!(Fup, Fdn, rhoup, rhodn, X, Y, xoff, yoff,
        layout, scratch, restricted)
    size(X.phi, 2) < 12 && return _add_coulomb_cross_direct_ordered!(Fup, Fdn,
        rhoup, rhodn, X, Y, xoff, yoff, layout, scratch, restricted)
    _add_coulomb_cross_direct_mul!(Fup, Fdn, rhoup, rhodn, X, Y, xoff, yoff,
        layout, scratch, restricted)
end
@inline function _exchange_scatter!(fup, fdn, mup, mdn, i, k, a, b, value,
        ::Val{restricted}) where restricted
    fup[i, b] = muladd(value, mup[k, a], fup[i, b])
    fup[i, a] = muladd(value, mup[k, b], fup[i, a])
    fup[k, b] = muladd(value, mup[i, a], fup[k, b])
    fup[k, a] = muladd(value, mup[i, b], fup[k, a])
    if !restricted
        fdn[i, b] = muladd(value, mdn[k, a], fdn[i, b])
        fdn[i, a] = muladd(value, mdn[k, b], fdn[i, a])
        fdn[k, b] = muladd(value, mdn[i, a], fdn[k, b])
        fdn[k, a] = muladd(value, mdn[i, b], fdn[k, a])
    end
end
@inline function _exchange_scatter_physical_diagonal!(fup, fdn, mup, mdn, i, k,
        a, value, ::Val{restricted}) where restricted
    fup[i, a] = muladd(value, mup[k, a], fup[i, a])
    fup[k, a] = muladd(value, mup[i, a], fup[k, a])
    if !restricted
        fdn[i, a] = muladd(value, mdn[k, a], fdn[i, a])
        fdn[k, a] = muladd(value, mdn[i, a], fdn[k, a])
    end
end
@inline function _exchange_scatter_pair_diagonal!(fup, fdn, mup, mdn, k, a, b,
        value, ::Val{restricted}) where restricted
    fup[k, b] = muladd(value, mup[k, a], fup[k, b])
    b != a && (fup[k, a] = muladd(value, mup[k, b], fup[k, a]))
    if !restricted
        fdn[k, b] = muladd(value, mdn[k, a], fdn[k, b])
        b != a && (fdn[k, a] = muladd(value, mdn[k, b], fdn[k, a]))
    end
end
function _add_coulomb_cross_exchange_ordered!(F, rho, X, Y, xoff, yoff, layout,
        scratch)
    kx, ky = size(X.phi, 2), size(Y.phi, 2)
    xrange, yrange = xoff + 1:xoff + kx, yoff + 1:yoff + ky
    for s in _active_slices(layout, Y.ra)
        P = _slice_phi(layout, Y.ra, Y.phi, s); d = size(P, 1)
        mixed = view(scratch.mixed, 1:kx, 1:d)
        field = view(scratch.exchange, 1:kx, 1:d)
        mul!(mixed, view(rho, xrange, yrange), transpose(P))
        W = view(X.W, :, _sym_state_range(X, s))
        for b = 1:d, i = 1:kx
            value = 0.0
            for a = 1:d, k = 1:kx
                value = muladd(_sym_value(W, i, k, a, b), mixed[k, a], value)
            end
            field[i, b] = value
        end
        mul!(view(F, xrange, yrange), field, P, -1.0, 1.0)
        mul!(view(F, yrange, xrange), transpose(P), transpose(field), -1.0, 1.0)
    end
end
function _add_coulomb_cross_exchange!(Fup, Fdn, rhoup, rhodn, X, Y, xoff,
        yoff, layout, scratch, restricted_value::Val{restricted}) where restricted
    kx, ky = size(X.phi, 2), size(Y.phi, 2)
    if kx < 10
        _add_coulomb_cross_exchange_ordered!(Fup, rhoup, X, Y, xoff, yoff,
            layout, scratch)
        restricted || _add_coulomb_cross_exchange_ordered!(Fdn, rhodn, X, Y,
            xoff, yoff, layout, scratch)
        return
    end
    xrange, yrange = xoff + 1:xoff + kx, yoff + 1:yoff + ky
    for s in _active_slices(layout, Y.ra)
        P = _slice_phi(layout, Y.ra, Y.phi, s); d = size(P, 1)
        mixedup = view(scratch.mixed, 1:kx, 1:d)
        fieldup = view(scratch.exchange, 1:kx, 1:d)
        mixeddn = view(scratch.mixed_dn, 1:kx, 1:d)
        fielddn = view(scratch.exchange_dn, 1:kx, 1:d)
        mul!(mixedup, view(rhoup, xrange, yrange), transpose(P)); fill!(fieldup, 0)
        if !restricted
            mul!(mixeddn, view(rhodn, xrange, yrange), transpose(P)); fill!(fielddn, 0)
        end
        W = view(X.W, :, _sym_state_range(X, s))
        A = 0; rt2 = sqrt(2.0)
        for b = 1:d, a = 1:b
            A += 1; I = 0
            for k = 1:kx
                for i = 1:k - 1
                    I += 1
                    value = W[I, A] / (a == b ? rt2 : 2.0)
                    if a == b
                        _exchange_scatter_physical_diagonal!(fieldup, fielddn,
                            mixedup, mixeddn, i, k, a, value, restricted_value)
                    else
                        _exchange_scatter!(fieldup, fielddn, mixedup, mixeddn,
                            i, k, a, b, value, restricted_value)
                    end
                end
                I += 1; value = W[I, A] / (a == b ? 1.0 : rt2)
                _exchange_scatter_pair_diagonal!(fieldup, fielddn, mixedup,
                    mixeddn, k, a, b, value, restricted_value)
            end
        end
        mul!(view(Fup, xrange, yrange), fieldup, P, -1.0, 1.0)
        mul!(view(Fup, yrange, xrange), transpose(P), transpose(fieldup), -1.0, 1.0)
        if !restricted
            mul!(view(Fdn, xrange, yrange), fielddn, P, -1.0, 1.0)
            mul!(view(Fdn, yrange, xrange), transpose(P), transpose(fielddn),
                -1.0, 1.0)
        end
    end
end
function _add_coulomb_cross!(Fup, Fdn, rhoup, rhodn, X, Y, xoff, yoff, layout,
        scratch, restricted::Val{r}) where r
    _add_coulomb_cross_direct!(Fup, Fdn, rhoup, rhodn, X, Y, xoff, yoff, layout,
        scratch, restricted)
    _add_coulomb_cross_exchange!(Fup, Fdn, rhoup, rhodn, X, Y, xoff, yoff,
        layout, scratch, restricted)
end
function _add_coulomb_fock!(Fup, Fdn, rhoup, rhodn, win, restricted)
    ml, lc, mr = win.ml, win.lc, win.mr; roff = ml + lc; dims = win.layout.dims
    _add_cached_sector!(Fup, Fdn, rhoup, rhodn, win.L.G, 0, 0, ml, ml,
        Val(:symmetric), restricted)
    coff = ml
    for s in win.center_slices
        ds = dims[s]; _add_sym_sector_both!(Fup, Fdn, rhoup, rhodn,
            view(win.L.W, :, _sym_state_range(win.L, s)), 0, coff, ml, ds, restricted)
        coff += ds
    end
    _add_coulomb_cross!(Fup, Fdn, rhoup, rhodn, win.L, win.R, 0, roff,
        win.layout, win.scratch, restricted)
    xoff = ml
    for sx in win.center_slices
        yoff = ml
        for sy in win.center_slices
            _add_cached_sector!(Fup, Fdn, rhoup, rhodn, _slice_vee(win.V, sx, sy),
                xoff, yoff, dims[sx], dims[sy], Val(:pair_grouped), restricted)
            yoff += dims[sy]
        end
        xoff += dims[sx]
    end
    coff = ml
    for s in win.center_slices
        ds = dims[s]; _add_sym_sector_both!(Fup, Fdn, rhoup, rhodn,
            view(win.R.W, :, _sym_state_range(win.R, s)), roff, coff, mr, ds,
            restricted); coff += ds
    end
    _add_cached_sector!(Fup, Fdn, rhoup, rhodn, win.R.G, roff, roff, mr, mr,
        Val(:symmetric), restricted)
end
vee_add_fock_r!(F, rho, win::_SlicedCoulombWindow) = begin
    win.ml + win.lc + win.mr == size(F, 1) || error("window dimensions do not match F")
    _add_coulomb_fock!(F, F, rho, rho, win, Val(true)) end
vee_add_fock!(Fup, Fdn, rhoup, rhodn, win::_SlicedCoulombWindow) = begin
    win.ml + win.lc + win.mr == size(Fup, 1) || error("window dimensions do not match F")
    _add_coulomb_fock!(Fup, Fdn, rhoup, rhodn, win, Val(false)) end
solve_hfdmrg(H, backend::_SlicedBasisBackendCachedCoulomb, psiup0; kwargs...) = solve_hfdmrg_core(H, backend, psiup0, psiup0; restricted = true, kwargs...)
function solve_hfdmrg(H, backend::_SlicedBasisBackendCachedCoulomb, psiup0, psidn0; kwargs...)
    solve_hfdmrg_core(H, backend, psiup0, psidn0; restricted = false, kwargs...) end
function solve_hfdmrg(Hup, Hdn, backend::_SlicedBasisBackendCachedCoulomb,
        psiup0, psidn0; kwargs...)
    Hup === Hdn && return solve_hfdmrg(Hup, backend, psiup0, psidn0; kwargs...)
    solve_hfdmrg_core_split(Hup, Hdn, backend, psiup0, psidn0; kwargs...)
end
