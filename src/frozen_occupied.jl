struct _FrozenRecord{T}; origin::Matrix{T}; work::Matrix{T}; side::Symbol; ra::UnitRange{Int}
    interval::Tuple{T,T}; discarded::T; rung::T; end
struct _FrozenLedger{T}
    generation::Int
    parent::Union{Nothing,_FrozenLedger{T}}
    up::Matrix{T}
    dn::Matrix{T}
    originup::Matrix{T}
    origindn::Matrix{T}
    uprecords::Vector{_FrozenRecord{T}}
    dnrecords::Vector{_FrozenRecord{T}}
end
struct _FrozenLeakageSchedule
    thresholds::NTuple{4,Float64}
    sweeps_per_rung::Int
    function _FrozenLeakageSchedule(sweeps_per_rung::Integer = 2)
        sweeps_per_rung > 0 || error("frozen schedule requires positive sweeps per rung")
        new((1e-4, 1e-6, 1e-8, 1e-10), sweeps_per_rung)
    end
end
struct _DensityFrozenField{T}
    j::Vector{T}
    kup::Matrix{T}
    kdn::Matrix{T}
end
struct _SlicedFrozenField
    Jpack::Vector{Float64}
    kup::Matrix{Float64}
    kdn::Matrix{Float64}
end
struct _FrozenCut{T,F}
    ledger::Union{Nothing,_FrozenLedger{T}}
    nup::Int
    ndn::Int
    field::F
    energy::T
    epoch::Int
end

mutable struct _FrozenPolicy{T,B}
    Hup
    Hdn
    backend::B
    restricted::Bool
    psiup::Matrix{T}
    psidn::Matrix{T}
    activeup::Matrix{T}
    activedn::Matrix{T}
    currentup::Matrix{T}
    currentdn::Matrix{T}
    protectedup::Matrix{T}
    protecteddn::Matrix{T}
    cuts::Vector{Any}
    generation::Int
    max_live::Int
    initializing::Bool
    rebuilt::Bool
    warmup::Matrix{T}
    warmdn::Matrix{T}
    schedule::_FrozenLeakageSchedule; rung::Int; rung_sweeps::Int; off_sweeps::Int
    draining::Bool
    scratch::Any
end
_frozen_size(p) = size(p.Hup, 1)
_retire_without_rebuild(backend) = false
_retire_without_rebuild(::DensityDensityBackend) = true
_frozen_backend_check(Hup, Hdn, backend, psiup, psidn, restricted, partition) =
    error("frozen_occupied=:roundoff_exact is unsupported for $(typeof(backend))")
function _frozen_backend_check(Hup, Hdn, backend::DensityDensityBackend,
        psiup, psidn, restricted, partition)
    V, T = backend.V, eltype(backend.V)
    T <: AbstractFloat ||
        error("roundoff-exact density interactions require AbstractFloat elements")
    all(isfinite, V) || error("roundoff-exact density interaction must be finite")
    scale = max(one(T), maximum(abs, V))
    maximum(abs, V .- V') <= 128eps(T) * scale ||
        error("roundoff-exact density interaction must be symmetric")
    T
end
function _frozen_policy(Hup, Hdn, backend, psiup, psidn, restricted, partition, schedule::_FrozenLeakageSchedule)
    T = _frozen_backend_check(Hup, Hdn, backend, psiup, psidn, restricted, partition)
    N = size(Hup, 1)
    z = zeros(T, N, 0)
    _FrozenPolicy(Hup, Hdn, backend, restricted, Matrix{T}(psiup), Matrix{T}(psidn),
        Matrix{T}(psiup), Matrix{T}(psidn), copy(z), copy(z), copy(z), copy(z),
        Any[], 0, 0, true, false, copy(z), copy(z), schedule, 1, 0, 0,
        false, nothing)
end
function _ledger_data(cut::_FrozenCut{T}, N, spin) where T
    @timeg "frozen ledger materialization" begin
    n = spin === :up ? cut.nup : cut.ndn
    C, O, p, node = zeros(T, N, n), zeros(T, N, n), n, cut.ledger
    while node !== nothing
        x = spin === :up ? node.up : node.dn
        o = spin === :up ? node.originup : node.origindn
        k = size(x, 2)
        if k > 0
            C[:, p - k + 1:p] = x
            O[:, p - k + 1:p] = o
        end
        p -= k
        node = node.parent
    end
    C, O
    end
end
_ledger_columns(cut, N, spin) = _ledger_data(cut, N, spin)[1]
_empty_frozen_field(backend::DensityDensityBackend, m) =
    _DensityFrozenField(zeros(eltype(backend.V), size(backend.V, 1)),
        zeros(eltype(backend.V), m, m), zeros(eltype(backend.V), m, m))
function _empty_frozen_cut(p::_FrozenPolicy{T}, m) where T
    _FrozenCut(nothing, 0, 0, _empty_frozen_field(p.backend, m), zero(T), p.rung)
end
function _frozen_begin_blocks!(p, nblocks, mleft, mright)
    p.cuts = Vector{Any}(undef, nblocks)
    p.cuts[1] = _empty_frozen_cut(p, mleft)
    p.cuts[end] = _empty_frozen_cut(p, mright)
    p.activeup, p.activedn = p.psiup, p.psidn
    p.currentup = similar(p.psiup, size(p.psiup, 1), 0)
    p.currentdn = similar(p.psidn, size(p.psidn, 1), 0)
    p.initializing = true
end
function _orthcols(X)
    size(X, 2) == 0 && return X
    F = svd(X)
    tol = 128eps(eltype(X)) * max(size(X)...) * max(1, F.S[1])
    r = count(>(tol), F.S)
    F.U[:, 1:r]
end
function _active_from_det(C, Cf)
    @timeg "frozen active recovery" begin
    size(Cf, 2) == 0 && return C
    F = svd(Cf' * C; full = true)
    tol = 128eps(eltype(C)) * max(size(C)...) * max(1, isempty(F.S) ? 1 : F.S[1])
    r = count(>(tol), F.S)
    C * F.V[:, r + 1:end]
    end
end
function _frozen_valid(C, CL, CR)
    T, N = eltype(C), size(C, 1)
    g = N * eps(T) / (1 - N * eps(T))
    Cf = hcat(CL, CR)
    scale = 128g * max(1, size(Cf, 2), norm(C, Inf))
    norm(Cf' * Cf - I) <= scale &&
        norm(Cf - C * (C' * Cf)) <= scale &&
        norm(CL' * CR) <= scale
end
function _all_frozen(p, spin)
    deepest = Dict{Int,Tuple{Int,Any}}()
    for i = eachindex(p.cuts)
        isassigned(p.cuts, i) || continue
        cut, node, depth = p.cuts[i], p.cuts[i].ledger, 0
        node === nothing && continue
        generation = node.generation; while node !== nothing
            node.generation == generation || error("mixed frozen ledger generation")
            depth += 1; node = node.parent
        end
        get(deepest, generation, (-1, nothing))[1] < depth &&
            (deepest[generation] = (depth, cut))
    end
    for i = eachindex(p.cuts)
        isassigned(p.cuts, i) || continue
        cut = p.cuts[i]; cut.ledger === nothing && continue
        best = deepest[cut.ledger.generation][2]
        cut.nup <= best.nup && cut.ndn <= best.ndn || error("nonnested frozen ledger counts")
        node = best.ledger
        while node !== nothing && node !== cut.ledger; node = node.parent; end
        node === cut.ledger || error("nonnested frozen ledger parent chain")
    end
    N, T = _frozen_size(p), eltype(p.psiup)
    cols = [_ledger_columns(deepest[g][2], N, spin) for g in sort!(collect(keys(deepest)))]
    isempty(cols) && return zeros(T, N, 0)
    _orthcols(hcat(cols...))
end
function _live_generations(p, in_flight = nothing)
    ids = Set(p.cuts[i].ledger.generation for i = eachindex(p.cuts)
        if isassigned(p.cuts, i) && p.cuts[i].ledger !== nothing)
    in_flight === nothing || push!(ids, in_flight.generation)
    length(ids)
end
function _reachable_records(p, spin)
    seen = IdDict{Any,Nothing}()
    records = _FrozenRecord{eltype(p.psiup)}[]
    for i = eachindex(p.cuts)
        isassigned(p.cuts, i) || continue
        node = p.cuts[i].ledger
        while node !== nothing
            if !haskey(seen, node)
                seen[node] = nothing
                append!(records, spin === :up ? node.uprecords : node.dnrecords)
            end
            node = node.parent
        end
    end
    records
end
function _cut_records(cut, spin)
    nodes, node = Any[], cut.ledger
    while node !== nothing
        push!(nodes, node)
        node = node.parent
    end
    records = _FrozenRecord{typeof(cut.energy)}[]
    for item in Iterators.reverse(nodes)
        append!(records, spin === :up ? item.uprecords : item.dnrecords)
    end
    records
end
function _record_columns(records, name, T, N)
    isempty(records) && return zeros(T, N, 0)
    hcat((getproperty(record, name) for record in records)...)
end
function _protected_record(record, protected)
    size(protected, 2) == 0 && return false
    W = record.work
    T, N = eltype(W), size(W, 1)
    err = norm(W - protected * (protected' * W))
    gate = 256N * eps(T) * max(1, norm(W))
    overlap = opnorm(protected' * W)
    err <= gate && return true
    overlap <= 256N * eps(T) ||
        error("protected subspace partially overlaps a frozen record")
    false
end
function _protection_addition(records, mask, T, N)
    selected = records[mask]
    isempty(selected) && return zeros(T, N, 0)
    _orthcols(_record_columns(selected, :work, T, N))
end
function _commit_protection!(p, up, dn)
    size(up, 2) > 0 && (p.protectedup = _orthcols(hcat(p.protectedup, up)))
    if p.restricted
        p.protecteddn = copy(p.protectedup)
    elseif size(dn, 2) > 0
        p.protecteddn = _orthcols(hcat(p.protecteddn, dn))
    end
end
_protect_stale!(p, up, dn, diagnostics) =
    error("slot-scoped frozen state fails closed on stale containment")
function _thin_projector_error(A, B)
    size(A, 1) == size(B, 1) || return Inf
    overlap = A' * B
    hypot(norm(B - A * overlap), norm(A - B * overlap'))
end
function _record_frozen_rebuild_error!(diagnostics, p, Cup, Cdn)
    diagnostics === nothing && return
    diagnostics.frozen_rebuild_projector_error = max(
        diagnostics.frozen_rebuild_projector_error,
        _thin_projector_error(p.psiup, Cup), _thin_projector_error(p.psidn, Cdn))
end
function _canonical_subspace(X)
    N, k = size(X); k == 0 && return X
    Q, R = zeros(eltype(X), N, k), Matrix(X); for j = 1:k
        pivot = argmax(vec(sum(abs2, R; dims = 2))); v = R * copy(@view R[pivot, :]); v -= Q[:, 1:j-1] * (Q[:, 1:j-1]' * v)
        nv = norm(v); nv > 128eps(eltype(X)) * max(size(X)...) * max(1, norm(X)) || error("frozen cluster lost rank during deterministic gauge")
        v ./= nv; v[argmax(abs.(v))] < 0 && (v .*= -1); Q[:, j] = v; R -= v * (v' * R)
    end
    Q
end
function _occupied_spectrum(C, ra)
    M = C[ra, :]' * C[ra, :]; all(isfinite, M) || error("nonfinite occupied spectrum")
    T, n = eltype(C), size(C, 2); scale, asym = one(T), zero(T); for j = 1:n, i = 1:j
        scale = max(scale, abs(M[i, j]), abs(M[j, i])); asym = max(asym, abs(M[i, j] - M[j, i]))
        M[i, j] = M[j, i] = (M[i, j] + M[j, i]) / 2
    end
    g = length(ra) * eps(T) / (1 - length(ra) * eps(T)); band = 128g * max(1, n, scale)
    asym <= band || error("occupied spectrum is nonsymmetric")
    E = eigen(Symmetric(M)); any(x -> x < -band || x > 1 + band, E.values) &&
        error("occupied spectrum lies resolvably outside [0,1]")
    clamp.(E.values, zero(T), one(T)), E.vectors, band
end
function _classify_frozen(p, C, ra, protected, side, spin)
    N, T = size(C, 1), eltype(C)
    empty = zeros(T, N, 0)
    p.draining && return C, empty, empty, _FrozenRecord{T}[]
    size(C, 2) == 0 && return C, empty, empty, _FrozenRecord{T}[]
    values, vectors, band = _occupied_spectrum(C, ra)
    tau = T(p.schedule.thresholds[p.rung])
    groups, firsti = UnitRange{Int}[], 1; for i = 2:length(values); values[i] - values[i - 1] > 2band && (push!(groups, firsti:i - 1); firsti = i); end
    push!(groups, firsti:length(values)); keep = trues(length(values))
    work, origins = Matrix{T}[], Matrix{T}[]
    records = _FrozenRecord{T}[]
    Pprot = size(protected, 2) == 0 ? protected :
        _orthcols(C * (C' * protected))
    for group in groups
        vals = @view values[group]; minimum(vals) > T(0.5) + band || continue
        leakage = one(T) .- vals; maximum(leakage) <= tau || continue
        origin = _canonical_subspace(C * @view(vectors[:, group]))
        size(Pprot, 2) > 0 && opnorm(Pprot' * origin) >= 1 - 128band && continue
        candidate = zeros(T, N, length(group)); candidate[ra, :] = origin[ra, :]; candidate = _canonical_subspace(candidate)
        any(W -> opnorm(W' * candidate) > 256N * eps(T), work) && continue
        record = _FrozenRecord(origin, candidate, side, ra,
            extrema(leakage), sum(leakage), tau)
        push!(work, candidate); push!(origins, origin); push!(records, record)
        keep[group] .= false
    end
    newwork = isempty(work) ? zeros(T, N, 0) : hcat(work...); neworigin = isempty(origins) ? zeros(T, N, 0) : hcat(origins...)
    C * @view(vectors[:, keep]), newwork, neworigin, records
end
function _demote_last!(records, work, origins)
    isempty(records) && error("structural demotion requires a frozen record")
    record = pop!(records)
    k = size(record.work, 2)
    work[:, 1:end-k], origins[:, 1:end-k], record.origin,
        _orthcols(hcat(record.origin, record.work))
end
function _release_records(C, current, records, protected)
    T, N = eltype(C), size(C, 1)
    drop = [_protected_record(record, protected) for record in records]
    keep = records[.!drop]
    dropped = records[drop]
    isempty(dropped) && return current, _active_from_det(C, current), keep,
        dropped, 0.0
    W = _record_columns(dropped, :work, T, N)
    keptcurrent = _active_from_det(current, W)
    active = _active_from_det(C, keptcurrent)
    candidate = hcat(keptcurrent, active)
    err = max(_thin_projector_error(candidate, C),
        norm(W - active * (active' * W)))
    gate = 256N * eps(T) * max(1, size(C, 2), norm(C))
    err <= gate || error("slot-scoped release lost determinant containment")
    keptcurrent, active, keep, dropped, err
end
function _frozen_growth!(p, side, oldblock, cra, oldidx, environment_cutoff)
    oldcut = p.cuts[oldidx]
    N = _frozen_size(p)
    if p.initializing
        currentup = _ledger_columns(oldcut, N, :up)
        currentdn = _ledger_columns(oldcut, N, :dn)
    else
        currentup, currentdn = p.currentup, p.currentdn
    end
    oldup, olddn = _cut_records(oldcut, :up), _cut_records(oldcut, :dn)
    currentup, activeup, retainedup, droppedup, errorup =
        _release_records(p.psiup, currentup, oldup, p.protectedup)
    if p.restricted
        currentdn, activedn, retaineddn, droppeddn, errordn =
            currentup, activeup, copy(retainedup), copy(droppedup), errorup
    else
        currentdn, activedn, retaineddn, droppeddn, errordn =
            _release_records(p.psidn, currentdn, olddn, p.protecteddn)
    end
    ra = side === :left ? (oldblock.ra[1]:last(cra)) :
        (first(cra):oldblock.ra[end])
    activeup, Cup, Oup, recordsup = _classify_frozen(
        p, activeup, ra, p.protectedup, side, :up)
    if p.restricted
        activedn, Cdn, Odn, recordsdn = activeup, Cup, Oup, copy(recordsup)
    else
        activedn, Cdn, Odn, recordsdn = _classify_frozen(
            p, activedn, ra, p.protecteddn, side, :dn)
    end
    protectedup, protecteddn = p.protectedup, p.protecteddn
    X = hcat(activeup[ra, :], activedn[ra, :])
    B = side === :left ?
        cat(oldblock.phi, Matrix{eltype(X)}(I, length(cra), length(cra)); dims = (1, 2)) :
        cat(Matrix{eltype(X)}(I, length(cra), length(cra)), oldblock.phi; dims = (1, 2))
    Y = B' * X
    F = svd(Y)
    k = findlast(>(environment_cutoff), F.S)
    if k === nothing
        Cf = hcat(p.currentup[ra, :], Cup[ra, :],
            p.currentdn[ra, :], Cdn[ra, :])
        FC = svd(Cf' * B; full = true)
        tol = 128eps(eltype(B)) * max(size(B)...) * max(1, isempty(FC.S) ? 1 : FC.S[1])
        r = count(>(tol), FC.S)
        if r == size(B, 2)
            if p.restricted
                Cup, Oup, c, protect = _demote_last!(recordsup, Cup, Oup)
                Cdn, Odn, _, _ = _demote_last!(recordsdn, Cdn, Odn)
                activeup = hcat(activeup, c); activedn = activeup
                protectedup = _orthcols(hcat(protectedup, protect))
                protecteddn = copy(protectedup)
            elseif size(Cup, 2) > 0
                Cup, Oup, c, protect = _demote_last!(recordsup, Cup, Oup)
                activeup = hcat(activeup, c)
                protectedup = _orthcols(hcat(protectedup, protect))
            else
                Cdn, Odn, c, protect = _demote_last!(recordsdn, Cdn, Odn)
                activedn = hcat(activedn, c)
                protecteddn = _orthcols(hcat(protecteddn, protect))
            end
            X = hcat(activeup[ra, :], activedn[ra, :])
            FX = svd(B' * X); kx = findlast(>(environment_cutoff), FX.S); kx === nothing && error("structural demotion failed to restore active rank")
            O = FX.U[:, 1:kx]
        else
            O = FC.V[:, r + 1:r + 1]
        end
    else
        O = F.U[:, 1:k]
    end
    currentup = hcat(currentup, Cup)
    currentdn = hcat(currentdn, Cdn)
    psiup = hcat(currentup, activeup)
    psidn = p.restricted ? psiup : hcat(currentdn, activedn)
    (; O, Cup, Cdn, Oup, Odn, recordsup, recordsdn, retainedup, retaineddn,
        droppedup, droppeddn, release_error=max(errorup, errordn),
        rebased=!isempty(droppedup) || !isempty(droppeddn), currentup,
        currentdn, activeup, activedn, psiup, psidn, protectedup, protecteddn)
end
function _density_frozen_fields(p, Cup, Cdn)
    V = p.backend.V
    Dup, Ddn = Cup * Cup', Cdn * Cdn'
    j = V * (diag(Dup) + diag(Ddn))
    j, Matrix(Diagonal(j) - V .* Dup), Matrix(Diagonal(j) - V .* Ddn)
end

function _frozen_absorption_map(backend, side, oldblock, cra, newblock)
    nold = length(oldblock.ra)
    side === :left ? oldblock.phi' * newblock.phi[1:nold, :] :
        oldblock.phi' * newblock.phi[end - nold + 1:end, :]
end
function _promote_frozen_field(backend::DensityDensityBackend, p, old, oldblock,
        newblock, growth, A, diagnostics)
    V, ra = backend.V, newblock.ra
    Pup, Pdn = growth.Cup * growth.Cup', growth.Cdn * growth.Cdn'
    jP, GPup, GPdn = _density_frozen_fields(p, growth.Cup, growth.Cdn)
    Cfu, Cfd = _ledger_columns(old, _frozen_size(p), :up),
        _ledger_columns(old, _frozen_size(p), :dn)
    _, Gfup, Gfdn = _density_frozen_fields(p, Cfu, Cfd)
    heffup, heffdn = p.Hup + Gfup, p.Hdn + Gfdn
    e = old.energy + tr(Pup * heffup) + tr(Pdn * heffdn) +
        0.5tr(Pup * GPup) + 0.5tr(Pdn * GPdn)
    f = old.field
    kup = A' * f.kup * A + newblock.phi' * (V[ra, ra] .* Pup[ra, ra]) * newblock.phi
    kdn = A' * f.kdn * A + newblock.phi' * (V[ra, ra] .* Pdn[ra, ra]) * newblock.phi
    _DensityFrozenField(f.j + jP, kup, kdn), e
end
function _density_frozen_fresh(backend::DensityDensityBackend, p, Cup, Cdn,
        phi, ra)
    V = backend.V
    Dup, Ddn = Cup * Cup', Cdn * Cdn'
    j, Gup, Gdn = _density_frozen_fields(p, Cup, Cdn)
    kup = phi' * (V[ra, ra] .* Dup[ra, ra]) * phi
    kdn = phi' * (V[ra, ra] .* Ddn[ra, ra]) * phi
    e = dot(Cup, p.Hup * Cup) + dot(Cdn, p.Hdn * Cdn) +
        0.5dot(Cup, Gup * Cup) + 0.5dot(Cdn, Gdn * Cdn)
    _DensityFrozenField(j, kup, kdn), e
end
_frozen_fresh(backend::DensityDensityBackend, p, Cup, Cdn, phi, ra) =
    _density_frozen_fresh(backend, p, Cup, Cdn, phi, ra)
_frozen_fresh(backend, p, Cup, Cdn, phi, ra) =
    _sliced_frozen_fresh(backend, p, Cup, Cdn, phi, ra)
function _frozen_store_growth!(p::_FrozenPolicy{T}, idx, oldidx, side, oldblock,
        newblock, growth, diagnostics) where T
    old = p.cuts[oldidx]
    allup = vcat(growth.retainedup, growth.recordsup)
    alldn = vcat(growth.retaineddn, growth.recordsdn)
    Cup = _record_columns(allup, :work, T, _frozen_size(p))
    Cdn = _record_columns(alldn, :work, T, _frozen_size(p))
    Oup = _record_columns(allup, :origin, T, _frozen_size(p))
    Odn = _record_columns(alldn, :origin, T, _frozen_size(p))
    firstnode = old.ledger === nothing
    hasnew = !isempty(growth.recordsup) || !isempty(growth.recordsdn)
    newroot = growth.rebased || firstnode && hasnew
    generation = newroot ? p.generation + 1 :
        old.ledger === nothing ? p.generation : old.ledger.generation
    if growth.rebased
        ledger = isempty(allup) && isempty(alldn) ? nothing :
            _FrozenLedger(generation, nothing, Cup, Cdn, Oup, Odn, allup, alldn)
        field, e = ledger === nothing ?
            (_empty_frozen_field(p.backend, size(newblock.phi, 2)), zero(T)) :
            @timeg "frozen field promotion" _frozen_fresh(
                p.backend, p, Cup, Cdn, newblock.phi, newblock.ra)
        if diagnostics !== nothing
            diagnostics.frozen_promotions_up += size(growth.Cup, 2)
            diagnostics.frozen_promotions_dn += size(growth.Cdn, 2)
        end
    else
        ledger = hasnew ? _FrozenLedger(generation,
            firstnode ? nothing : old.ledger, growth.Cup, growth.Cdn,
            growth.Oup, growth.Odn, growth.recordsup, growth.recordsdn) :
            old.ledger
        cra = side === :left ? (last(oldblock.ra) + 1:last(newblock.ra)) :
            (first(newblock.ra):first(oldblock.ra) - 1)
        A = _frozen_absorption_map(p.backend, side, oldblock, cra, newblock)
        field, e = @timeg "frozen field promotion" begin
            _promote_frozen_field(
                p.backend, p, old, oldblock, newblock, growth, A, diagnostics)
        end
    end
    nup, ndn = size(Cup, 2), size(Cdn, 2)
    p.max_live = max(p.max_live, _live_generations(p, ledger))
    epoch = p.draining ? length(p.schedule.thresholds) + 1 : p.rung
    cut = _FrozenCut(ledger, nup, ndn, field, e, epoch)
    newroot && (p.generation = generation)
    p.cuts[idx] = cut
    p.currentup, p.currentdn = growth.currentup, growth.currentdn
    p.activeup, p.activedn = growth.activeup, growth.activedn
    p.psiup, p.psidn = growth.psiup, growth.psidn
    p.protectedup, p.protecteddn = growth.protectedup, growth.protecteddn
    if diagnostics !== nothing
        diagnostics.frozen_window_entries_up += length(growth.droppedup)
        diagnostics.frozen_window_entries_dn += length(growth.droppeddn)
        diagnostics.frozen_window_entry_error = max(
            diagnostics.frozen_window_entry_error, growth.release_error)
    end
    if _detailed(diagnostics)
        push!(diagnostics.frozen_cuts, (; ordinal=length(diagnostics.frozen_cuts)+1,
            side, cut=idx, physical_range=newblock.ra, active_rank=size(newblock.phi,2),
            frozen_up=cut.nup, frozen_dn=cut.ndn,
            rung=epoch,
            threshold=p.draining ? :off : p.schedule.thresholds[p.rung],
            leakage_up=getproperty.(growth.recordsup, :interval),
            leakage_dn=getproperty.(growth.recordsdn, :interval)))
    end
end

function _window_basis(N, left, cra, right)
    @timeg "frozen window basis" begin
    B = zeros(promote_type(eltype(left.phi), eltype(right.phi)), N,
        left.m + length(cra) + right.m)
    B[left.ra, 1:left.m] = left.phi
    B[cra, left.m + 1:left.m + length(cra)] =
        Matrix{eltype(B)}(I, length(cra), length(cra))
    B[right.ra, end - right.m + 1:end] = right.phi
    B
    end
end

function _allowed_basis(Cf, B)
    @timeg "frozen allowed basis" begin
    size(Cf, 2) == 0 && return Matrix{eltype(B)}(I, size(B, 2), size(B, 2))
    A = Cf' * B
    F = svd(A; full = true)
    g = max(size(A)...) * eps(eltype(B)) / (1 - max(size(A)...) * eps(eltype(B)))
    tol = 128g * max(1, norm(A))
    r = count(>(tol), F.S)
    Z = F.V[:, r + 1:end]
    norm(A * Z) <= tol * max(1, sqrt(size(Z, 2))) ||
        error("failed to construct spin-specific frozen complement")
    Z
    end
end

function _frozen_effective(p, left, right, B)
    L, R = p.cuts[left], p.cuts[right]
    _frozen_effective(p.backend, p, L, R, B)
end
function _frozen_effective(backend::DensityDensityBackend, p, L, R, B)
    N = size(B, 1)
    CLu, CLd = _ledger_columns(L, N, :up), _ledger_columns(L, N, :dn)
    CRu, CRd = _ledger_columns(R, N, :up), _ledger_columns(R, N, :dn)
    lf, rf = L.field, R.field
    j = lf.j + rf.j
    Gup = B' * (j .* B)
    Gdn = copy(Gup)
    ml, mr = size(lf.kup, 1), size(rf.kup, 1)
    Gup[1:ml, 1:ml] -= lf.kup
    Gdn[1:ml, 1:ml] -= lf.kdn
    Gup[end - mr + 1:end, end - mr + 1:end] -= rf.kup
    Gdn[end - mr + 1:end, end - mr + 1:end] -= rf.kdn
    qL = diag(CLu * CLu' + CLd * CLd')
    qR = diag(CRu * CRu' + CRd * CRd')
    ef = L.energy + R.energy + 0.5(dot(qL, rf.j) + dot(qR, lf.j))
    Gup, Gdn, ef, hcat(CLu, CRu), hcat(CLd, CRd)
end

function _frozen_local!(p, li, ri, left, cra, right, Hup, Hdn, win, lambda,
        scf_cutoff, relative, diagnostics)
    _frozen_window_check(p.backend, win)
    B = _window_basis(_frozen_size(p), left, cra, right)
    Gup, Gdn, ef, Cfu, Cfd =
        @timeg "frozen effective field" _frozen_effective(p, li, ri, B)
    Lu, _ = _ledger_data(p.cuts[li], size(B, 1), :up)
    Ld, _ = _ledger_data(p.cuts[li], size(B, 1), :dn)
    Ru, _ = _ledger_data(p.cuts[ri], size(B, 1), :up)
    Rd, _ = _ledger_data(p.cuts[ri], size(B, 1), :dn)
    staleup = !_frozen_valid(p.psiup, Lu, Ru)
    staledn = !p.restricted && !_frozen_valid(p.psidn, Ld, Rd)
    (staleup || staledn) &&
        return (; stale = true, up = staleup ? Cfu : Cfu[:, 1:0],
            dn = staledn ? Cfd : Cfd[:, 1:0])
    Heffup, Heffdn = Hup + Gup, Hdn + Gdn
    Zup, Zdn = _allowed_basis(Cfu, B), _allowed_basis(Cfd, B)
    Aup = _active_from_det(p.psiup, Cfu)
    Adn = p.restricted ? Aup : _active_from_det(p.psidn, Cfd)
    nup, ndn = size(Aup, 2), size(Adn, 2)
    if diagnostics !== nothing
        diagnostics.frozen_zero_up_windows += nup == 0
        diagnostics.frozen_zero_dn_windows += ndn == 0
        diagnostics.frozen_zero_both_windows += nup == 0 && ndn == 0
    end
    nup <= size(Zup, 2) && ndn <= size(Zdn, 2) ||
        error("active occupation exceeds its spin-specific allowed space")
    Uup, Udn = Zup' * (B' * Aup), Zdn' * (B' * Adn)
    rup, rdn = Uup * Uup', Udn * Udn'
    psiup, psidn = Zup * Uup, Zdn * Udn
    energies, rises = Float64[], Int[]
    energy, lastenergy, converged, nit = ef, 1e10, false, 0
    if nup == 0 && ndn == 0
        converged = true
    else
        rhoup, rhodn = Zup * rup * Zup', Zdn * rdn * Zdn'
        Fup = copy(Heffup)
        if p.restricted
            @timeg "Fock build" vee_add_fock_r!(Fup, rhoup, win)
        else
            Fdn = copy(Heffdn)
            @timeg "Fock build" vee_add_fock!(Fup, Fdn, rhoup, rhodn, win)
        end
        for s = 1:4
            if p.restricted
                nup > 0 && (Uup = eigsym(Zup' * Fup * Zup)[2][:, 1:nup])
                psiup = Zup * Uup
                rup = (1 - lambda) * rup + lambda * (Uup * Uup')
                rhoup = Zup * rup * Zup'
                Fup = copy(Heffup)
                @timeg "Fock build" vee_add_fock_r!(Fup, rhoup, win)
                energy = ef + tr(rhoup * (Fup + Heffup))
                psidn = psiup
            else
                nup > 0 && (Uup = eigsym(Zup' * Fup * Zup)[2][:, 1:nup])
                ndn > 0 && (Udn = eigsym(Zdn' * Fdn * Zdn)[2][:, 1:ndn])
                psiup, psidn = Zup * Uup, Zdn * Udn
                rup = (1 - lambda) * rup + lambda * (Uup * Uup')
                rdn = (1 - lambda) * rdn + lambda * (Udn * Udn')
                rhoup, rhodn = Zup * rup * Zup', Zdn * rdn * Zdn'
                Fup, Fdn = copy(Heffup), copy(Heffdn)
                @timeg "Fock build" vee_add_fock!(Fup, Fdn, rhoup, rhodn, win)
                energy = ef + 0.5tr(rhoup * (Fup + Heffup)) +
                    0.5tr(rhodn * (Fdn + Heffdn))
            end
            nit = s
            push!(energies, energy)
            energy > lastenergy && (push!(rises, s); lambda *= 0.5)
            converged = relative ? _rel_converged(lastenergy, energy, scf_cutoff) :
                abs(lastenergy - energy) < scf_cutoff
            (converged || s == 4) && break
            lastenergy = energy
        end
    end
    p.currentup, p.currentdn = Cfu, Cfd
    p.activeup, p.activedn = B * psiup, p.restricted ? B * psiup : B * psidn
    p.psiup = hcat(Cfu, p.activeup)
    p.psidn = p.restricted ? p.psiup : hcat(Cfd, p.activedn)
    _frozen_valid(p.psiup, Cfu, Cfu[:, 1:0]) ||
        error("reassembled alpha determinant failed frozen containment")
    p.restricted || _frozen_valid(p.psidn, Cfd, Cfd[:, 1:0]) ||
        error("reassembled beta determinant failed frozen containment")
    (; stale = false, psiup, psidn, energy, nit, converged, lambda, rises, energies)
end

_frozen_residual_norm(F, C, Cf) = opnorm(F * Cf - C * (C' * (F * Cf)))

_frozen_window_check(backend, win) = nothing
_frozen_audit_context(p, spin, backend::DensityDensityBackend) = spin === :up ?
    _DensityFockAuditContext(p.Hup, backend.V, p.psiup, p.psidn) :
    _DensityFockAuditContext(p.Hdn, backend.V, p.psidn, p.psiup)

_frozen_payload(p) = 0
_frozen_field_empty(f::_DensityFrozenField) = all(iszero, f.j) &&
    all(iszero, f.kup) && all(iszero, f.kdn)
function _frozen_drained(p)
    all(eachindex(p.cuts)) do i
        !isassigned(p.cuts, i) || begin
            c = p.cuts[i]
            c.ledger === nothing && c.nup == 0 && c.ndn == 0 &&
                iszero(c.energy) && _frozen_field_empty(c.field)
        end
    end
end
function _residual_gate(ctx, C, X)
    Y = similar(X); _audit_action!(Y, ctx, X); N, T = size(C, 1), eltype(C); g = N * eps(T) / (1 - N * eps(T)); value, gate = _vnorm(_q(C, Y)), 128g * max(1, ctx.scale)
    value, gate, value <= gate
end
function _audit_result(ctx, C, warm, enabled)
    enabled || return (; outcome = :deferred, reason = :occupied_residual, warm, actions = ctx.actions, rhs = ctx.rhs)
    result = _occupation_order(ctx, C, warm); (result.reason === :action_failure || hasproperty(result, :hot) && (result.hot.reason === :action_failure || result.cold.reason === :action_failure)) && error("matrix-free occupation action failed")
    result
end
function _record_passes(record, tau); values, _, band = _occupied_spectrum(record.origin, record.ra)
    minimum(values) > 0.5 + band && maximum(1 .- values) <= tau; end
_thaw_columns(records, mask) = sum((size(records[i].work, 2) for i in eachindex(records) if mask[i]); init = 0)
function _frozen_audit!(p, diagnostics)
    @timeg "frozen sweep audit" begin
    p.rung_sweeps += 1
    if p.draining
        drained = _frozen_drained(p)
        _detailed(diagnostics) && push!(diagnostics.frozen_events,
            (; rung=:off, rung_sweep=p.rung_sweeps, threshold=:off,
                drained, retained_payload=_frozen_payload(p)))
        return drained ? :drained : :continue
    end
    boundary = p.rung_sweeps == p.schedule.sweeps_per_rung
    Cfu, Cfd = _all_frozen(p, :up), _all_frozen(p, :dn)
    diagnostics === nothing || (diagnostics.frozen_max_generations = p.max_live; diagnostics.frozen_live_generations = _live_generations(p))
    up, dn = _reachable_records(p, :up), _reachable_records(p, :dn)
    staleup = !_frozen_valid(p.psiup, Cfu, Cfu[:, 1:0]); staledn = !p.restricted && !_frozen_valid(p.psidn, Cfd, Cfd[:, 1:0])
    cup, cdn = @timeg "frozen audit context" begin
        (_frozen_audit_context(p, :up, p.backend),
            _frozen_audit_context(p, :dn, p.backend))
    end
    rfu, gfu, passfu = _residual_gate(cup, p.psiup, Cfu); rfd, gfd, passfd = p.restricted ? (rfu, gfu, passfu) :
        _residual_gate(cdn, p.psidn, Cfd)
    passfu &= !staleup; passfd &= !staledn
    thawup = fill(!passfu && !isempty(up), length(up))
    thawdn = fill(!passfd && !isempty(dn), length(dn))
    diagnostics === nothing || passfu || isempty(up) || (diagnostics.frozen_residual_thaws += 1);
    diagnostics === nothing || passfd || isempty(dn) || (diagnostics.frozen_residual_thaws += 1)
    upresult = dnresult = nothing
    if boundary
        rcu, gcu, passcu = _residual_gate(cup, p.psiup, p.psiup); rcd, gcd, passcd = p.restricted ? (rcu, gcu, passcu) :
            _residual_gate(cdn, p.psidn, p.psidn)
        upresult, dnresult = @timeg "frozen occupation audit" begin
            upvalue = _audit_result(cup, p.psiup, p.warmup, passcu)
            dnvalue = p.restricted ? upvalue :
                _audit_result(cdn, p.psidn, p.warmdn, passcd)
            (upvalue, dnvalue)
        end
        p.warmup = upresult.warm; p.warmdn = dnresult.warm
        for (result, records, mask, spin) in ((upresult, up, thawup, :up),
                (dnresult, dn, thawdn, :dn))
            thaw = result.outcome === :inverted || (result.outcome === :indeterminate && result.reason !== :occupied_residual)
            thaw && fill!(mask, true)
            if thaw && diagnostics !== nothing && !isempty(records)
                result.outcome === :inverted ? (diagnostics.frozen_inversion_thaws += 1) :
                    (diagnostics.frozen_indeterminate_thaws += 1)
            end
        end
        nexttau = p.rung < length(p.schedule.thresholds) ? p.schedule.thresholds[p.rung + 1] : -1.0
        scheduledup = [nexttau < 0 || !_record_passes(c, nexttau) for c in up]
        scheduleddn = [nexttau < 0 || !_record_passes(c, nexttau) for c in dn]
        diagnostics === nothing || (diagnostics.frozen_scheduled_thaws +=
            _thaw_columns(up, scheduledup) + _thaw_columns(dn, scheduleddn))
        thawup .|= scheduledup; thawdn .|= scheduleddn
    end
    payload = diagnostics === nothing ? 0 : _frozen_payload(p); diagnostics === nothing || (diagnostics.frozen_retained_values = max(diagnostics.frozen_retained_values, payload))
    if _detailed(diagnostics)
        push!(diagnostics.frozen_events, (; rung=p.rung, rung_sweep=p.rung_sweeps, threshold=p.schedule.thresholds[p.rung],
            frozen_residual=(rfu,rfd), frozen_gate=(gfu,gfd), complete_residual=boundary ? (rcu,rcd) : nothing, complete_gate=boundary ? (gcu,gcd) : nothing,
            occupation=boundary ? ((upresult.outcome,upresult.reason),
                (dnresult.outcome,dnresult.reason)) : nothing, retained_payload=payload))
    end
    _retire_without_rebuild(p.backend) ||
        error("slot-scoped retirement requires a nonrebuilding backend")
    @timeg "frozen transition handling" begin
        T, N = eltype(p.psiup), _frozen_size(p)
        protectup = _protection_addition(up, thawup, T, N)
        protectdn = p.restricted ? protectup :
            _protection_addition(dn, thawdn, T, N)
        _commit_protection!(p, protectup, protectdn)
        if boundary
            p.rung_sweeps = 0
            if p.rung == length(p.schedule.thresholds)
                p.draining = true
            else
                p.rung += 1
            end
        end
    end
    :continue
    end
end
function _frozen_final_audit!(p, Cup, Cdn, diagnostics)
    @timeg "frozen terminal audit" begin
    p.psiup, p.psidn = Cup, p.restricted ? Cup : Cdn
    cup, cdn = @timeg "frozen audit context" begin
        (_frozen_audit_context(p, :up, p.backend),
            _frozen_audit_context(p, :dn, p.backend))
    end
    ru, gu, pu = _residual_gate(cup,p.psiup,p.psiup); rd, gd, pd = p.restricted ? (ru,gu,pu) : _residual_gate(cdn,p.psidn,p.psidn)
    up, dn = @timeg "frozen occupation audit" begin
        upvalue = _audit_result(cup,p.psiup,p.warmup,pu)
        dnvalue = p.restricted ? upvalue :
            _audit_result(cdn,p.psidn,p.warmdn,pd)
        (upvalue, dnvalue)
    end
    p.warmup, p.warmdn = up.warm, dn.warm
    _detailed(diagnostics) && push!(diagnostics.frozen_events, (; rung=:off, rung_sweep=p.off_sweeps, threshold=:off, complete_residual=(ru,rd), complete_gate=(gu,gd), occupation=((up.outcome,up.reason),(dn.outcome,dn.reason))))
    pu && pd && up.outcome === :ordered && dn.outcome === :ordered
    end
end
