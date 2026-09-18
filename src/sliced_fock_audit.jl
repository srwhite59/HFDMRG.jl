const _AUDIT_LIMITS = (capacity = 12, subspace = 40, restarts = 8, actions = 128, probes = 12)
const _AUDIT_RANK_FACTOR, _AUDIT_GEOMETRY_FACTOR = 128.0, 512.0
mutable struct _DensityFockAuditContext{T,TH,TV,TC}
    H::TH; V::TV; C::TC
    j::Vector{T}; diag::Vector{T}; scale::T
    actions::Int; rhs::Int
end
function _DensityFockAuditContext(H, V, C, Cother)
    T, N = eltype(V), size(H, 1)
    q = vec(sum(abs2, C; dims = 2) + sum(abs2, Cother; dims = 2))
    j, d, s = V * q, zeros(T, N), zero(T)
    for p = 1:N
        row = zero(T)
        for r = 1:N
            f = H[p, r] - V[p, r] * dot(@view(C[p, :]), @view(C[r, :]))
            p == r && (f += j[p]; d[p] = f)
            row += abs(f)
        end
        s = max(s, row)
    end
    all(isfinite, j) && all(isfinite, d) && isfinite(s) ||
        error("nonfinite density physical Fock stream")
    _DensityFockAuditContext(H, V, C, j, d, s, 0, 0)
end
function _audit_action!(Y, c::_DensityFockAuditContext, X)
    c.actions += 1; c.rhs += size(X, 2); all(isfinite, X) || error("nonfinite physical Fock action input")
    mul!(Y, c.H, X); Y .+= c.j .* X
    tmp = similar(Y)
    for i = axes(c.C, 2)
        ci = @view c.C[:, i]
        mul!(tmp, c.V, ci .* X); Y .-= ci .* tmp
    end
    all(isfinite, Y) || error("nonfinite physical Fock action output")
    Y
end
mutable struct _SlicedFockAuditContext{TH,TV,TC,TB}
    H::TH; V::TV; layout::SliceLayout; C::TC
    Jpack::Vector{Float64}; poffs::Vector{Int}; diag::Vector{Float64}
    scale::Float64; blocks::TB; setup_reason::Symbol
    setup_seconds::Float64; setup_bytes::Int
    neighbor_values::Int; setup_pairs::Int
    actions::Int; rhs::Int
end
function _sliced_exchange(Vnm, Cn, Cm)
    dn, dm = size(Cn, 1), size(Cm, 1); D, K = Cn * Cm', zeros(Float64, dn, dm)
    @inbounds for a = 1:dn, d = 1:dm, b = 1:dn, cc = 1:dm
        K[a, d] += Vnm[a, b, cc, d] * D[b, cc]
    end
    K
end
function _sliced_hartree_data(layout, V, Cup, Cdn)
    dims, offs, S = layout.dims, layout.offs, nslices(layout)
    poffs = [0; cumsum(dims .^ 2)]; Jpack = zeros(Float64, poffs[end])
    Ttotal = [(@view Cup[offs[m] + 1:offs[m + 1], :]) *
        (@view Cup[offs[m] + 1:offs[m + 1], :])' +
        (@view Cdn[offs[m] + 1:offs[m + 1], :]) *
        (@view Cdn[offs[m] + 1:offs[m + 1], :])' for m = 1:S]
    for n = 1:S, m = 1:S
        Jn = reshape(@view(Jpack[poffs[n] + 1:poffs[n + 1]]), dims[n], dims[n])
        Vnm = _slice_vee(V, n, m)
        @inbounds for a = 1:dims[n], b = 1:dims[n], cc = 1:dims[m], d = 1:dims[m]
            Jn[a, b] += Vnm[a, b, cc, d] * Ttotal[m][cc, d]
        end
    end
    Jpack, poffs
end
function _sliced_audit_data(H, backend, Cup, Cdn, spin)
    layout, V = backend.layout, backend.V
    dims, offs, S = layout.dims, layout.offs, nslices(layout)
    N, C = offs[end], spin === :up ? Cup : Cdn
    Jpack, poffs = _sliced_hartree_data(layout, V, Cup, Cdn)
    diagonal = [zeros(Float64, d, d) for d in dims]
    upper = [zeros(Float64, dims[n], dims[n + 1]) for n = 1:S - 1]
    lower = [zeros(Float64, dims[n + 1], dims[n]) for n = 1:S - 1]
    rows, diag = zeros(Float64, N), zeros(Float64, N)
    for n = 1:S, m = 1:S
        rn, rm = offs[n] + 1:offs[n + 1], offs[m] + 1:offs[m + 1]
        K = _sliced_exchange(_slice_vee(V, n, m), @view(C[rn, :]), @view(C[rm, :]))
        Jn = reshape(@view(Jpack[poffs[n] + 1:poffs[n + 1]]), dims[n], dims[n])
        for a = 1:dims[n], d = 1:dims[m]
            f = H[offs[n] + a, offs[m] + d] - K[a, d]
            n == m && (f += Jn[a, d])
            rows[offs[n] + a] += abs(f)
            n == m && a == d && (diag[offs[n] + a] = f)
            n == m && (diagonal[n][a, d] = f)
            m == n + 1 && (upper[n][a, d] = f)
            n == m + 1 && (lower[m][a, d] = f)
        end
    end
    blocks = (; diagonal, upper, lower)
    finite = all(isfinite, Jpack) && all(isfinite, rows) &&
        all(B -> all(isfinite, B), vcat(diagonal, upper, lower))
    scale = maximum(rows; init = 0.0); asymmetry = 0.0
    if finite
        for B in diagonal; asymmetry = max(asymmetry, opnorm(B - B', Inf)); end
        for n = 1:S - 1
            asymmetry = max(asymmetry, opnorm(upper[n] - lower[n]', Inf))
        end
    end
    g = N * eps(Float64) / (1 - N * eps(Float64))
    reason = 9maximum(dims) + 2 > 40 ? :unsupported_layout :
        !finite ? :nonfinite_blocks :
        asymmetry > 128g * max(1, scale) ? :nonhermitian_blocks : :ok
    if reason === :ok
        for n = 1:S; diagonal[n] = (diagonal[n] + diagonal[n]') / 2; end
        for n = 1:S - 1
            upper[n] = (upper[n] + lower[n]') / 2; lower[n] = upper[n]'
        end
    end
    values = sum(length, diagonal) + sum(length, upper) + sum(length, lower)
    (; C, Jpack, poffs, diag, scale, blocks, reason, values, pairs = S^2)
end
function _SlicedFockAuditContext(H, backend, Cup, Cdn, spin)
    timed = @timed _sliced_audit_data(H, backend, Cup, Cdn, spin); d = timed.value
    _SlicedFockAuditContext(H, backend.V, backend.layout, d.C, d.Jpack, d.poffs,
        d.diag, d.scale, d.blocks, d.reason, timed.time, timed.bytes,
        d.values, d.pairs, 0, 0)
end
function _sliced_occupied_action!(Y, H, V, layout, C, Jpack, poffs, X)
    all(isfinite, X) || error("nonfinite physical Fock action input")
    mul!(Y, H, X)
    dims, offs, S = layout.dims, layout.offs, nslices(layout)
    for n = 1:S
        rn = offs[n] + 1:offs[n + 1]
        Jn = reshape(@view(Jpack[poffs[n] + 1:poffs[n + 1]]), dims[n], dims[n])
        mul!(@view(Y[rn, :]), Jn, @view(X[rn, :]), 1.0, 1.0)
        for m = 1:S
            rm = offs[m] + 1:offs[m + 1]
            K = _sliced_exchange(_slice_vee(V, n, m), @view(C[rn, :]), @view(C[rm, :]))
            mul!(@view(Y[rn, :]), K, @view(X[rm, :]), -1.0, 1.0)
        end
    end
    all(isfinite, Y) || error("nonfinite physical Fock action output")
    Y
end
function _audit_action!(Y, c::_SlicedFockAuditContext, X)
    c.actions += 1; c.rhs += size(X, 2)
    _sliced_occupied_action!(Y, c.H, c.V, c.layout, c.C, c.Jpack, c.poffs, X)
end
_q(C, X) = X - C * (C' * X)
_vnorm(X) = size(X, 2) == 0 ? 0.0 : opnorm(X)
_audit_field(x, name, default) = hasproperty(x, name) ? getproperty(x, name) : default
function _fixed_orth(C, X)
    size(X, 2) == 0 && return similar(C, size(C, 1), 0)
    Y = _q(C, Matrix(X)); Y = _q(C, Y)
    all(isfinite, Y) || return similar(C, size(C, 1), 0)
    F = try svd(Y) catch; return similar(C, size(C, 1), 0) end; isempty(F.S) && return similar(C, size(C, 1), 0)
    tol = _AUDIT_RANK_FACTOR * eps(eltype(C)) * max(size(Y)...) * max(1, F.S[1])
    F.U[:, 1:count(>(tol), F.S)]
end
function _mix64(x::UInt64)
    z = x + 0x9e3779b97f4a7c15; z = xor(z, z >> 30) * 0xbf58476d1ce4e5b9
    z = xor(z, z >> 27) * 0x94d049bb133111eb
    xor(z, z >> 31)
end
function _signprobe(T, N, run, j)
    key = run === :hot ? 0x243f6a8885a308d3 : 0x13198a2e03707344
    [((_mix64(xor(key, _mix64(UInt64(j)), _mix64(UInt64(i)))) >> 63) == 0 ?
        -one(T) : one(T)) / sqrt(T(N)) for i = 1:N]
end
function _act(ctx, X, limits)
    ctx.actions < limits.actions || return nothing, :action_cap
    try; _audit_action!(similar(X), ctx, X), :ok
    catch; nothing, :action_failure end
end
function _geometry(C, V, W)
    g = size(C, 1) * eps(eltype(C)) / (1 - size(C, 1) * eps(eltype(C)))
    virtuality = _vnorm(C' * W); search = _vnorm(V' * W)
    gram = size(W, 2) == 0 ? 0.0 : opnorm(W' * W - I)
    (; ok = all(isfinite, (virtuality, search, gram)) &&
        max(virtuality, search, gram) <= _AUDIT_GEOMETRY_FACTOR * g,
        virtuality, search, gram, gate = _AUDIT_GEOMETRY_FACTOR * g)
end
function _relative_corrections(C, V, X)
    Y = Matrix(X)
    for _ = 1:2; Y = _q(C, Y); Y -= V * (V' * Y); end
    all(isfinite, Y) || return nothing, :nonfinite_correction
    F = try svd(Y) catch; return nothing, :correction_failure end; isempty(F.S) && return nothing, :deficient_correction
    tol = _AUDIT_RANK_FACTOR * eps(eltype(C)) * max(size(Y)...) * F.S[1]
    r = count(>(tol), F.S); r > 0 || return nothing, :deficient_correction
    W = F.U[:, 1:r]; for _ = 1:2; W = _q(C, W); W -= V * (V' * W); end
    W = try Matrix(qr(W).Q[:, 1:r]) catch; return nothing, :correction_failure end; geometry = _geometry(C, V, W)
    geometry.ok ? (W, :ok) : (nothing, :correction_geometry)
end
function _structured_probes(ctx, C, run; modes = 9, hashes = 2)
    dims, offs, S = ctx.layout.dims, ctx.layout.offs, nslices(ctx.layout)
    dmax, N = maximum(dims), size(C, 1); ncols = modes * dmax + hashes
    ncols <= 40 || return similar(C, N, 0)
    X = zeros(eltype(C), N, ncols); col = 0
    for a = 1:dmax, k = 1:modes
        col += 1
        for s = 1:S
            a > dims[s] && continue
            x = run === :hot ? s / (S + 1) : (s - 0.5) / S
            X[offs[s] + a, col] = run === :hot ? sinpi(k * x) : cospi((k - 0.5) * x)
        end
    end
    for j = 1:hashes; X[:, col + j] = _signprobe(eltype(C), N, run, j); end
    X
end
function _sliced_start(ctx, C, run, limits; probes = _structured_probes(ctx, C, run))
    nv = size(C, 1) - size(C, 2); X = _fixed_orth(C, probes)
    size(X, 2) >= min(4, nv) || return (; ok = false, reason = :deficient_coarse_start)
    size(X, 2) <= limits.subspace || return (; ok = false, reason = :subspace_cap)
    Y, reason = _act(ctx, X, limits); Y === nothing && return (; ok = false, reason)
    Y = _q(C, Y); E = try eigen(Symmetric(X' * Y)) catch; return (; ok = false, reason = :coarse_failure) end; order = sortperm(E.values)
    k = min(4, nv, size(X, 2)); S = E.vectors[:, order[1:k]]
    (; ok = true, V = X * S, AV = Y * S, coarse_rank = size(X, 2))
end
_push_seed(C, B, x) = ((X = _fixed_orth(C, hcat(B, x))); (X, size(X, 2) > size(B, 2)))
function _density_start(ctx, C, Rocc, warm, run, tau, limits)
    T, N, nv = eltype(C), size(C, 1), size(C, 1) - size(C, 2)
    run === :hot && !all(isfinite, warm) && return (; ok = false, reason = :nonfinite_warm)
    B = run === :hot ? _fixed_orth(C, warm) : similar(C, N, 0)
    E = try svd(Rocc) catch; return (; ok = false, reason = :residual_failure) end; rr = count(>(tau / 16), E.S)
    rr > limits.capacity && return (; ok = false, reason = :residual_seed_cap)
    for i = 1:rr; B, ok = _push_seed(C, B, @view(E.U[:, i:i])); ok || return (; ok = false, reason = :deficient_residual); end
    required = size(B, 2); target = min(nv, max(4, required + 4))
    target <= limits.capacity || return (; ok = false, reason = :seed_cap)
    nh = 0
    for j = 1:limits.probes
        B, ok = _push_seed(C, B, _signprobe(T, N, run, j)); nh += ok
        (nh == min(2, nv - required) || size(B, 2) == nv) && break
    end
    nh == min(2, nv - required) || return (; ok = false, reason = :deficient_hash)
    order = sortperm(eachindex(ctx.diag); by = i -> (ctx.diag[i], run === :hot ? i : -i))
    run === :cold && (order = order[[1 + mod(k * (N - 1), N) for k = 0:N - 1]])
    nc = 0
    for i in order
        e = zeros(T, N); e[i] = 1; B, ok = _push_seed(C, B, e); nc += ok
        (nc >= min(2, nv - required - nh) && size(B, 2) >= target || size(B, 2) == nv) && break
    end
    size(B, 2) == target || size(B, 2) == nv || return (; ok = false, reason = :deficient_start)
    size(B, 2) <= limits.subspace || return (; ok = false, reason = :subspace_cap)
    AV, reason = _act(ctx, B, limits); AV === nothing && return (; ok = false, reason)
    (; ok = true, V = B, AV = _q(C, AV), coarse_rank = size(B, 2))
end
struct _DiagonalAuditPreconditioner{T}
    diag::Vector{T}; delta::T
end
struct _NeighborAuditFactor{TF,TU,TM,TO}
    factors::TF; upper::TU; multipliers::TM; offs::TO
    min_rcond::Float64; values::Int
end
function _reciprocal_condition(A)
    try
        s = svdvals(A); isempty(s) || !isfinite(s[1]) || s[1] == 0 ? 0.0 : s[end] / s[1]
    catch; 0.0
    end
end
function _neighbor_factor(blocks, offs, shift; factorize = lu)
    isfinite(shift) || return nothing, :nonfinite_shift
    S = length(blocks.diagonal); factors = Vector{Any}(undef, S)
    multipliers = Vector{Matrix{Float64}}(undef, S - 1); minrc = Inf
    try
        schur = blocks.diagonal[1] - shift * I
        for n = 1:S
            schur = (schur + schur') / 2; rc = _reciprocal_condition(schur)
            minrc = min(minrc, rc); rc >= sqrt(eps(Float64)) || return nothing, :factor_condition
            factors[n] = factorize(schur)
            n == S && break
            multipliers[n] = transpose(factors[n] \ transpose(blocks.lower[n]))
            all(isfinite, multipliers[n]) || return nothing, :nonfinite_factor
            schur = blocks.diagonal[n + 1] - shift * I - multipliers[n] * blocks.upper[n]
        end
    catch
        return nothing, :factorization_failure
    end
    values = sum(length, blocks.diagonal) + sum(length, multipliers)
    _NeighborAuditFactor(factors, blocks.upper, multipliers, offs, minrc, values), :ok
end
function _precondition(p::_DiagonalAuditPreconditioner, R, theta)
    X = similar(R)
    for j = axes(R, 2)
        den = p.diag .- theta[j]
        den = [abs(x) < p.delta ? copysign(p.delta, iszero(x) ? one(x) : x) : x for x in den]
        X[:, j] = R[:, j] ./ den
    end
    all(isfinite, X) ? X : nothing
end
function _precondition(p::_NeighborAuditFactor, R, theta)
    Y, S = copy(R), length(p.factors)
    for n = 2:S
        rn, rp = p.offs[n] + 1:p.offs[n + 1], p.offs[n - 1] + 1:p.offs[n]
        Y[rn, :] -= p.multipliers[n - 1] * Y[rp, :]
    end
    X = similar(R); rn = p.offs[S] + 1:p.offs[S + 1]
    X[rn, :] = p.factors[S] \ Y[rn, :]
    for n = S - 1:-1:1
        rn, rr = p.offs[n] + 1:p.offs[n + 1], p.offs[n + 1] + 1:p.offs[n + 2]
        X[rn, :] = p.factors[n] \ (Y[rn, :] - p.upper[n] * X[rr, :])
    end
    all(isfinite, X) ? X : nothing
end
function _normalized_ritz(C, Z, AZ)
    for i = axes(Z, 2)
        z, az = _q(C, @view(Z[:, i])), _q(C, @view(AZ[:, i]))
        nz = norm(z); isfinite(nz) && nz > 0 || return nothing
        Z[:, i], AZ[:, i] = z / nz, az / nz
    end
    _geometry(C, similar(C, size(C, 1), 0), Z).ok ? (Z, AZ) : nothing
end
function _refine(ctx, C, V, AV, occmax, tau, preconditioner, limits)
    nv, restarts, maxspace, largest = size(C, 1) - size(C, 2), 0, size(V, 2), size(V, 2)
    while true
        all(isfinite, V) && all(isfinite, AV) || return (; outcome = :indeterminate, reason = :nonfinite_iteration, restarts, maxspace, largest)
        size(V, 2) <= limits.subspace || return (; outcome = :indeterminate, reason = :subspace_cap, restarts, maxspace, largest)
        E = try eigen(Symmetric(V' * AV)) catch; return (; outcome = :indeterminate, reason = :iteration_failure, restarts, maxspace, largest) end; order = sortperm(E.values); S = E.vectors[:, order]
        normalized = _normalized_ritz(C, V * S, AV * S)
        normalized === nothing && return (; outcome = :indeterminate, reason = :ritz_geometry, restarts, maxspace, largest)
        Z, AZ = normalized; theta = [dot(@view(Z[:, i]), @view(AZ[:, i])) for i = axes(Z, 2)]
        residuals = [norm(@view(AZ[:, i]) - theta[i] * @view(Z[:, i])) for i = eachindex(theta)]
        all(isfinite, theta) && all(isfinite, residuals) || return (; outcome = :indeterminate, reason = :nonfinite_iteration, restarts, maxspace, largest)
        order = sortperm(theta); theta, residuals = theta[order], residuals[order]
        Z, AZ = Z[:, order], AZ[:, order]
        witness = findfirst(<(occmax - tau), theta)
        witness === nothing || return (; outcome = :inverted, reason = :rayleigh_witness,
            witness = theta[witness], residual = residuals[witness], restarts, maxspace, largest)
        lower, upper = theta - residuals, theta + residuals
        m = count(<=(theta[1] + 4tau), theta); boundary = m + 1
        wholecluster = size(V, 2) == nv && m == nv
        boundaryok = boundary <= length(theta) && residuals[boundary] <= tau / 8 &&
            lower[boundary] > maximum(@view upper[1:m])
        complete = maximum(@view residuals[1:m]) <= tau / 8 && (wholecluster || boundaryok)
        if complete
            m <= limits.capacity || return (; outcome = :indeterminate, reason = :cluster_cap, restarts, maxspace, largest)
            lower[1] >= occmax - tau || return (; outcome = :indeterminate, reason = :threshold_overlap, restarts, maxspace, largest)
            gap = wholecluster ? Inf : lower[boundary] - maximum(@view upper[1:m])
            return (; outcome = :candidate, reason = :converged, interval = (lower[1], upper[1]),
                values = theta[1:m], cluster = Z[:, 1:m], residual = maximum(@view residuals[1:m]),
                gap, restarts, maxspace, largest)
        end
        k = min(4, length(theta)); R = AZ[:, 1:k] - Z[:, 1:k] * Diagonal(theta[1:k])
        P = try _precondition(preconditioner, R, theta[1:k]) catch; nothing end
        P === nothing && return (; outcome = :indeterminate, reason = :preconditioner_failure, restarts, maxspace, largest)
        W, reason = _relative_corrections(C, V, P)
        W === nothing && return (; outcome = :indeterminate, reason, restarts, maxspace, largest)
        if size(V, 2) + size(W, 2) > limits.subspace
            restarts += 1
            restarts <= limits.restarts || return (; outcome = :indeterminate, reason = :restart_cap, restarts, maxspace, largest)
            V = _fixed_orth(C, Z[:, 1:min(12, length(theta))])
            size(V, 2) > 0 && _geometry(C, similar(C, size(C, 1), 0), V).ok ||
                return (; outcome = :indeterminate, reason = :restart_geometry, restarts, maxspace, largest)
            AV, reason = _act(ctx, V, limits); AV === nothing || (AV = _q(C, AV))
        else
            AW, reason = _act(ctx, W, limits); AV = AW === nothing ? nothing : hcat(AV, _q(C, AW)); V = hcat(V, W)
        end
        AV === nothing && return (; outcome = :indeterminate, reason, restarts, maxspace, largest)
        maxspace = max(maxspace, size(V, 2)); largest = max(largest, size(V, 2), size(W, 2))
    end
end
_projector_distance(A, B) = max(opnorm(B - A * (A' * B)), opnorm(A - B * (B' * A)))
function _combine_streams(hot, cold, tau, g)
    hot.outcome === cold.outcome || return (; outcome = :indeterminate, reason = :outcome_disagreement)
    hot.outcome === :inverted && return (; outcome = :inverted, reason = :confirmed_witness,
        witness = min(hot.witness, cold.witness))
    (hot.outcome === :candidate && cold.outcome === :candidate) ||
        return (; outcome = :indeterminate, reason = :stream_failure)
    length(hot.values) == length(cold.values) || return (; outcome = :indeterminate, reason = :cluster_dimension)
    maximum(abs, hot.values - cold.values) <= tau / 4 || return (; outcome = :indeterminate, reason = :cluster_values)
    gap = min(hot.gap, cold.gap); (gap > 0 || isinf(gap)) || return (; outcome = :indeterminate, reason = :projector_gap)
    gate = isinf(gap) ? 128g : 4(hot.residual + cold.residual) / gap + 128g
    gate < 1 / 8 || return (; outcome = :indeterminate, reason = :projector_gap)
    pdiff = try _projector_distance(hot.cluster, cold.cluster) catch; Inf end
    isfinite(pdiff) && pdiff <= gate || return (; outcome = :indeterminate, reason = :projector_disagreement, projector_difference = pdiff)
    (; outcome = :ordered, reason = :converged, interval = hot.interval,
        values = hot.values, cluster = hot.cluster, projector_difference = pdiff,
        projector_gate = gate)
end
function _stream(ctx::_DensityFockAuditContext, C, R, warm, run, occmax, tau, p, limits)
    a0, r0 = ctx.actions, ctx.rhs; start = _density_start(ctx, C, R, warm, run, tau, limits)
    start.ok || return (; outcome = :indeterminate, reason = start.reason, actions_used = ctx.actions - a0, rhs_used = ctx.rhs - r0)
    result = _refine(ctx, C, start.V, start.AV, occmax, tau, p, limits)
    merge(result, (; actions_used = ctx.actions - a0, rhs_used = ctx.rhs - r0,
        largest = max(_audit_field(result, :largest, 0), start.coarse_rank), coarse_rank = start.coarse_rank))
end
function _stream(ctx::_SlicedFockAuditContext, C, R, warm, run, occmax, tau, p, limits)
    a0, r0 = ctx.actions, ctx.rhs; start = _sliced_start(ctx, C, run, limits)
    start.ok || return (; outcome = :indeterminate, reason = start.reason, actions_used = ctx.actions - a0, rhs_used = ctx.rhs - r0)
    result = _refine(ctx, C, start.V, start.AV, occmax, tau, p, limits)
    merge(result, (; actions_used = ctx.actions - a0, rhs_used = ctx.rhs - r0,
        largest = max(_audit_field(result, :largest, 0), start.coarse_rank), coarse_rank = start.coarse_rank))
end
_prepare(ctx::_DensityFockAuditContext, occmax, tau) = (; ok = true,
    preconditioner = _DiagonalAuditPreconditioner(ctx.diag,
        max(tau, sqrt(eps(eltype(ctx.diag))) * max(1, ctx.scale))), info = (;))
function _prepare(ctx::_SlicedFockAuditContext, occmax, tau)
    ctx.setup_reason === :ok || return (; ok = false, reason = ctx.setup_reason)
    sneighbor = maximum(opnorm, vcat(ctx.blocks.upper, ctx.blocks.lower); init = 0.0)
    delta2 = max(1.0, sneighbor); shift = occmax - sqrt(delta2)
    timed = @timed _neighbor_factor(ctx.blocks, ctx.layout.offs, shift)
    factor, reason = timed.value; factor === nothing && return (; ok = false, reason)
    info = (; sneighbor, delta2, shift, min_rcond = factor.min_rcond,
        factor_seconds = timed.time, factor_bytes = timed.bytes,
        factor_values = factor.values, setup_seconds = ctx.setup_seconds,
        setup_bytes = ctx.setup_bytes, neighbor_values = ctx.neighbor_values,
        setup_pairs = ctx.setup_pairs, nearest_blocks = 3nslices(ctx.layout) - 2)
    (; ok = true, preconditioner = factor, info)
end
function _occupation_order(ctx, C, warm = similar(C, size(C, 1), 0); limits = _AUDIT_LIMITS)
    T, N, n = eltype(C), size(C, 1), size(C, 2); emptywarm = similar(C, N, 0)
    N * eps(T) < 0.5 || error("physical dimension exceeds exact audit roundoff model")
    all(isfinite, C) || error("nonfinite occupied determinant")
    g = N * eps(T) / (1 - N * eps(T)); tau = 128g * max(1, ctx.scale)
    norm(C' * C - I) <= 128g * max(1, n) || error("nonorthonormal audit determinant")
    ctx isa _SlicedFockAuditContext && ctx.setup_reason !== :ok && return (;
        outcome = :indeterminate, reason = ctx.setup_reason, warm = emptywarm, tau,
        actions = ctx.actions, rhs = ctx.rhs)
    (n == 0 || n == N) && return (; outcome = :ordered, reason = :vacuous,
        warm = emptywarm, tau, occupied_residual = zero(T), actions = ctx.actions, rhs = ctx.rhs)
    FC, reason = _act(ctx, C, limits)
    FC === nothing && return (; outcome = :indeterminate, reason, warm = emptywarm,
        tau, actions = ctx.actions, rhs = ctx.rhs)
    R = _q(C, FC); rocc = _vnorm(R); occmax = eigmax(Symmetric(C' * FC))
    rocc <= tau || return (; outcome = :indeterminate, reason = :occupied_residual,
        warm = emptywarm, tau, occupied = occmax, occupied_residual = rocc,
        actions = ctx.actions, rhs = ctx.rhs)
    prep = _prepare(ctx, occmax, tau)
    prep.ok || return (; outcome = :indeterminate, reason = prep.reason, warm = emptywarm,
        tau, occupied = occmax, occupied_residual = rocc, actions = ctx.actions, rhs = ctx.rhs)
    th = @timed(_stream(ctx, C, R, warm, :hot, occmax, tau, prep.preconditioner, limits)); tc = @timed(_stream(ctx, C, R, emptywarm, :cold, occmax, tau, prep.preconditioner, limits))
    hot, cold = th.value, tc.value
    combined = _combine_streams(hot, cold, tau, g)
    resultwarm = combined.outcome === :ordered ? combined.cluster : emptywarm
    (; combined..., warm = resultwarm, tau, occupied = occmax,
        occupied_residual = rocc, actions = ctx.actions, rhs = ctx.rhs, hot, cold,
        restarts = _audit_field(hot, :restarts, 0) + _audit_field(cold, :restarts, 0),
        maxspace = max(_audit_field(hot, :maxspace, 0), _audit_field(cold, :maxspace, 0)),
        largest = max(_audit_field(hot, :largest, 0), _audit_field(cold, :largest, 0)),
        stream_seconds = (hot = th.time, cold = tc.time), stream_bytes = (hot = th.bytes, cold = tc.bytes),
        prep.info...)
end
