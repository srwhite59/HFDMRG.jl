const _SlicedHistoryBackend = Union{SlicedBasisBackend,SlicedBasisBackendCached}

_sliced_history_finite(V::Array{Float64,6}) = all(isfinite, V)
function _sliced_history_finite(vee::SlicedVeeRagged)
    all(block -> all(isfinite, block), Iterators.flatten(vee.V))
end

function _sliced_history_policy_valid(p)
    p.bootstrap_sweeps > 0 && length(p.captures) >= 2 &&
        issorted(p.captures) && allunique(p.captures) &&
        last(p.captures) == p.bootstrap_sweeps && first(p.captures) > 0 &&
        p.fifo_length >= length(p.captures) && p.cleanup_sweeps > 0 &&
        p.max_cycles >= 0 && p.target >= 0 && p.diis_iterations > 0 &&
        p.diis_memory > 1 && 0 <= p.minimum_overlap <= 1
end

function _history_inputs(H, backend::_SlicedHistoryBackend, C, policy)
    H isa Matrix{Float64} && C isa Matrix{Float64} ||
        error("sliced history RHF requires Matrix{Float64} H and C0")
    N, n = size(C)
    size(H) == (N, N) && backend.layout.offs[end] == N && 0 < n < N ||
        error("inconsistent sliced history RHF sizes")
    all(isfinite, H) && all(isfinite, C) && _sliced_history_finite(backend.V) ||
        error("nonfinite sliced history RHF input")
    norm(H - H', Inf) <= 64eps(Float64) * N * max(1, norm(H, Inf)) ||
        error("sliced history RHF H must be symmetric")
    norm(C' * C - I) <= 32eps(Float64) * sqrt(N * n) ||
        error("sliced history RHF orbitals must be orthonormal")
    _sliced_history_policy_valid(policy) || error("invalid history RHF policy")
    nothing
end

function _history_uhf_inputs(Hup, Hdn, backend::_SlicedHistoryBackend,
        Cup, Cdn, policy)
    all(A -> A isa Matrix{Float64}, (Hup, Hdn, Cup, Cdn)) ||
        error("sliced history UHF requires Matrix{Float64} inputs")
    N, nup, ndn = size(Hup, 1), size(Cup, 2), size(Cdn, 2)
    size(Hup) == size(Hdn) == (N, N) && backend.layout.offs[end] == N &&
        size(Cup, 1) == size(Cdn, 1) == N && 0 < nup + ndn &&
        nup < N && ndn < N || error("inconsistent sliced history UHF sizes")
    all(A -> all(isfinite, A), (Hup, Hdn, Cup, Cdn)) &&
        _sliced_history_finite(backend.V) || error("nonfinite sliced history UHF input")
    for H in (Hup, Hdn)
        norm(H - H', Inf) <= 64eps(Float64) * N * max(1, norm(H, Inf)) ||
            error("sliced history UHF H must be symmetric")
    end
    for C in (Cup, Cdn)
        _history_gram(C) <= 32eps(Float64) * sqrt(length(C)) ||
            error("sliced history UHF orbitals must be orthonormal")
    end
    _sliced_history_policy_valid(policy) || error("invalid history UHF policy")
    nothing
end

function _sliced_history_plan(layout, r, matrix_budget)
    try
        r2 = Base.checked_mul(r, r)
        r4 = Base.checked_mul(r2, r2)
        N2 = Base.checked_mul(layout.offs[end], layout.offs[end])
        pairs = Base.checked_mul.(layout.dims, layout.dims)
        S2 = foldl(Base.checked_add, pairs; init = 0)
        payload = Base.checked_mul(S2, S2)
        limit = min(payload, Base.checked_mul(matrix_budget, N2))
        base = Base.checked_add(Base.checked_mul(S2, r2), r4)
        perrow = Base.checked_add(S2, r2)
        available = Base.checked_sub(limit, base)
        batch = min(max(256, maximum(pairs)), available ÷ perrow)
        r4 <= N2 && batch >= maximum(pairs) || return nothing
        (; r2, r4, pairs, S2, batch)
    catch error
        error isa OverflowError || rethrow()
        nothing
    end
end

function _sliced_history_contract(backend::_SlicedHistoryBackend, Q, matrix_budget)
    layout = backend.layout
    N, r = size(Q)
    N == layout.offs[end] || error("history basis does not match sliced layout")
    plan = _sliced_history_plan(layout, r, matrix_budget)
    plan === nothing && return nothing
    poffs = [0; cumsum(plan.pairs)]
    R = Matrix{Float64}(undef, plan.S2, plan.r2)
    for s = 1:nslices(layout)
        _pair_map_mat!(@view(R[poffs[s] + 1:poffs[s + 1], :]),
            @view(Q[layout.offs[s] + 1:layout.offs[s + 1], :]))
    end
    Vrows = zeros(Float64, plan.batch, plan.S2)
    Yrows = zeros(Float64, plan.batch, plan.r2)
    G = zeros(Float64, plan.r2, plan.r2)
    firsts = 1
    while firsts <= nslices(layout)
        lasts, rows = firsts, plan.pairs[firsts]
        while lasts < nslices(layout) && rows + plan.pairs[lasts + 1] <= plan.batch
            lasts += 1
            rows += plan.pairs[lasts]
        end
        Vb, Yb = @view(Vrows[1:rows, :]), @view(Yrows[1:rows, :])
        for s = firsts:lasts, t = 1:nslices(layout)
            rr = poffs[s] - poffs[firsts] + 1:poffs[s + 1] - poffs[firsts]
            cc = poffs[t] + 1:poffs[t + 1]
            copyto!(@view(Vb[rr, cc]),
                reshape(_slice_vee(backend.V, s, t), length(rr), length(cc)))
        end
        mul!(Yb, Vb, R)
        mul!(G, transpose(@view(R[poffs[firsts] + 1:poffs[lasts + 1], :])),
            Yb, 1.0, 1.0)
        firsts = lasts + 1
    end
    reshape(G, r, r, r, r)
end

function _history_model(H, backend::_SlicedHistoryBackend, Q)
    G = _sliced_history_contract(backend, Q, 2)
    G === nothing ? nothing : (; H = Q' * H * Q, G)
end

function _history_uhf_model(Hup, Hdn, backend::_SlicedHistoryBackend, Q)
    G = _sliced_history_contract(backend, Q, 4)
    G === nothing && return nothing
    Hqup = Q' * Hup * Q
    (; Hup = Hqup, Hdn = Hdn === Hup ? Hqup : Q' * Hdn * Q, G)
end

function _history_model_metrics(model::NamedTuple{(:H, :G)}, C)
    D, r = C * C', size(C, 1)
    F = model.H + 2reshape(reshape(model.G, r^2, r^2) * vec(D), r, r)
    for b = 1:r, a = 1:r, d = 1:r, c = 1:r
        F[a, b] -= model.G[a, c, d, b] * D[c, d]
    end
    FC = F * C; R, K = FC - C * (C' * FC), F * D - D * F
    (; F, energy = dot(D, F + model.H), residual = norm(R),
       commutator = K, kresidual = norm(K))
end

function _history_uhf_model_metrics(model::NamedTuple{(:Hup, :Hdn, :G)}, Cup, Cdn)
    Dup, Ddn, r = Cup * Cup', Cdn * Cdn', size(Cup, 1)
    Fup, Fdn = copy(model.Hup), copy(model.Hdn)
    for b = 1:r, a = 1:r, d = 1:r, c = 1:r
        direct = model.G[a, b, c, d] * (Dup[c, d] + Ddn[c, d])
        Fup[a, b] += direct - model.G[a, c, d, b] * Dup[c, d]
        Fdn[a, b] += direct - model.G[a, c, d, b] * Ddn[c, d]
    end
    FupC, FdnC = Fup * Cup, Fdn * Cdn
    Rup, Rdn = FupC - Cup * (Cup' * FupC), FdnC - Cdn * (Cdn' * FdnC)
    Kup, Kdn = Fup * Dup - Dup * Fup, Fdn * Ddn - Ddn * Fdn
    (; Fup, Fdn, energy = 0.5dot(Dup, Fup + model.Hup) + 0.5dot(Ddn, Fdn + model.Hdn),
       residual_up = norm(Rup), residual_dn = norm(Rdn), Kup, Kdn,
       kresidual = hypot(norm(Kup), norm(Kdn)) / sqrt(2))
end

function _history_physical(H, backend::_SlicedHistoryBackend, C)
    D, F = C * C', copy(H)
    sliced_add_fock_r!(F, D, backend.V, backend.layout)
    FC = F * C
    (; energy = dot(D, F + H), residual = norm(FC - C * (C' * FC)),
       gram = norm(C' * C - I))
end

function _history_uhf_physical(Hup, Hdn, backend::_SlicedHistoryBackend, Cup, Cdn)
    Dup, Ddn = Cup * Cup', Cdn * Cdn'
    Fup, Fdn = copy(Hup), copy(Hdn)
    sliced_add_fock_uhf!(Fup, Fdn, Dup, Ddn, backend.V, backend.layout)
    FupC, FdnC = Fup * Cup, Fdn * Cdn
    Rup, Rdn = FupC - Cup * (Cup' * FupC), FdnC - Cdn * (Cdn' * FdnC)
    rup, rdn = norm(Rup), norm(Rdn)
    sz = (size(Cup, 2) - size(Cdn, 2)) / 2
    (; energy = 0.5dot(Dup, Fup + Hup) + 0.5dot(Ddn, Fdn + Hdn),
       residual_up = rup, residual_dn = rdn, residual = max(rup, rdn),
       residual_rms = hypot(rup, rdn) / sqrt(2), gram_up = _history_gram(Cup),
       gram_dn = _history_gram(Cdn),
       s2 = sz^2 + (size(Cup, 2) + size(Cdn, 2)) / 2 - norm(Cup' * Cdn)^2)
end
