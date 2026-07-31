struct _FrozenLedger{T}
    generation::Int
    parent::Union{Nothing,_FrozenLedger{T}}
    up::Matrix{T}
    dn::Matrix{T}
end
struct _FrozenCut{T}
    ledger::Union{Nothing,_FrozenLedger{T}}
    nup::Int
    ndn::Int
    j::Vector{T}
    kup::Matrix{T}
    kdn::Matrix{T}
    energy::T
end
mutable struct _FrozenPolicy{T}
    Hup
    Hdn
    V
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
end
function _frozen_policy(Hup, Hdn, backend, psiup, psidn, restricted)
    backend isa DensityDensityBackend ||
        error("frozen_occupied=:roundoff_exact is supported only by density-density F1")
    V, T = backend.V, eltype(backend.V)
    T <: AbstractFloat ||
        error("roundoff-exact density interactions require AbstractFloat elements")
    all(isfinite, V) || error("roundoff-exact density interaction must be finite")
    scale = max(one(T), maximum(abs, V))
    maximum(abs, V .- V') <= 128eps(T) * scale ||
        error("roundoff-exact density interaction must be symmetric")
    N = size(Hup, 1)
    z = zeros(T, N, 0)
    _FrozenPolicy(Hup, Hdn, V, restricted, Matrix{T}(psiup), Matrix{T}(psidn),
        Matrix{T}(psiup), Matrix{T}(psidn), copy(z), copy(z), copy(z), copy(z),
        Any[], 0, 0, true, false)
end
function _ledger_columns(cut::_FrozenCut{T}, N, spin) where T
    n = spin === :up ? cut.nup : cut.ndn
    C, p, node = zeros(T, N, n), n, cut.ledger
    while node !== nothing
        x = spin === :up ? node.up : node.dn
        k = size(x, 2)
        k > 0 && (C[:, p - k + 1:p] = x)
        p -= k
        node = node.parent
    end
    C
end
function _empty_frozen_cut(p::_FrozenPolicy{T}, m) where T
    _FrozenCut{T}(nothing, 0, 0, zeros(T, size(p.V, 1)),
        zeros(T, m, m), zeros(T, m, m), zero(T))
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
    size(Cf, 2) == 0 && return C
    F = svd(Cf' * C; full = true)
    tol = 128eps(eltype(C)) * max(size(C)...) * max(1, isempty(F.S) ? 1 : F.S[1])
    r = count(>(tol), F.S)
    C * F.V[:, r + 1:end]
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
    cols = [_ledger_columns(p.cuts[i], size(p.V, 1), spin)
        for i = eachindex(p.cuts) if isassigned(p.cuts, i)]
    _orthcols(hcat(cols...))
end
function _live_generations(p, in_flight = nothing)
    ids = Set(p.cuts[i].ledger.generation for i = eachindex(p.cuts)
        if isassigned(p.cuts, i) && p.cuts[i].ledger !== nothing)
    in_flight === nothing || push!(ids, in_flight.generation)
    length(ids)
end
function _protect_stale!(p, up, dn, diagnostics)
    p.protectedup = _orthcols(hcat(p.protectedup, p.psiup * (p.psiup' * up)))
    p.protecteddn = _orthcols(hcat(p.protecteddn, p.psidn * (p.psidn' * dn)))
    p.rebuilt = true
    diagnostics === nothing || (diagnostics.frozen_rebuilds += 1)
end
function _classify_frozen(C, ra, protected)
    N, n = size(C)
    n == 0 && return C, similar(C, N, 0)
    E = eigen(Symmetric(C[ra, :]' * C[ra, :]))
    T = eltype(C)
    g = length(ra) * eps(T) / (1 - length(ra) * eps(T))
    tl = 128g * max(1, n)
    any(x -> x < -tl || x > 1 + tl, E.values) &&
        error("restricted occupied spectrum lies resolvably outside [0,1]")
    rotated = C * E.vectors
    Pprot = size(protected, 2) == 0 ? protected :
        _orthcols(C * (C' * protected))
    frozen = Int[]
    outside = setdiff(1:N, ra)
    for i = 1:n
        tail = norm(@view rotated[outside, i])
        tt = 128eps(T) * max(1, n) * sqrt(max(1, length(outside))) *
             max(1, norm(C, Inf))
        protected_i = size(Pprot, 2) > 0 &&
            norm(Pprot' * @view(rotated[:, i])) >= 1 - tl
        1 - clamp(E.values[i], zero(T), one(T)) <= tl && tail <= tt &&
            !protected_i && push!(frozen, i)
    end
    keep = setdiff(1:n, frozen)
    Cp = zeros(T, N, length(frozen))
    for (j, i) in enumerate(frozen)
        Cp[ra, j] = rotated[ra, i]
        Cp[:, j] ./= norm(@view Cp[:, j])
    end
    rotated[:, keep], Cp
end
function _frozen_growth!(p, side, oldblock, cra, oldidx, environment_cutoff)
    oldcut = p.cuts[oldidx]
    N = size(p.V, 1)
    if p.initializing
        p.currentup = _ledger_columns(oldcut, N, :up)
        p.currentdn = _ledger_columns(oldcut, N, :dn)
        p.activeup = _active_from_det(p.psiup, p.currentup)
        p.activedn = p.restricted ? p.activeup :
            _active_from_det(p.psidn, p.currentdn)
    end
    ra = side === :left ? (oldblock.ra[1]:last(cra)) :
        (first(cra):oldblock.ra[end])
    p.activeup, Cup = _classify_frozen(p.activeup, ra, p.protectedup)
    if p.restricted
        p.activedn, Cdn = p.activeup, Cup
    else
        p.activedn, Cdn = _classify_frozen(p.activedn, ra, p.protecteddn)
    end
    X = hcat(p.activeup[ra, :], p.activedn[ra, :])
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
                c = Cup[:, end:end]
                Cup, Cdn = Cup[:, 1:end - 1], Cdn[:, 1:end - 1]
                p.activeup = hcat(p.activeup, c)
                p.activedn = p.activeup
            elseif size(Cup, 2) > 0
                p.activeup = hcat(p.activeup, Cup[:, end:end])
                Cup = Cup[:, 1:end - 1]
            else
                p.activedn = hcat(p.activedn, Cdn[:, end:end])
                Cdn = Cdn[:, 1:end - 1]
            end
            X = hcat(p.activeup[ra, :], p.activedn[ra, :])
            O = svd(B' * X).U[:, 1:1]
        else
            O = FC.V[:, r + 1:r + 1]
        end
    else
        O = F.U[:, 1:k]
    end
    p.currentup = hcat(p.currentup, Cup)
    p.currentdn = hcat(p.currentdn, Cdn)
    p.psiup = hcat(p.currentup, p.activeup)
    p.psidn = p.restricted ? p.psiup : hcat(p.currentdn, p.activedn)
    (; O, Cup, Cdn)
end
function _frozen_fields(p, Cup, Cdn)
    Dup, Ddn = Cup * Cup', Cdn * Cdn'
    j = p.V * (diag(Dup) + diag(Ddn))
    j, Matrix(Diagonal(j) - p.V .* Dup), Matrix(Diagonal(j) - p.V .* Ddn)
end

function _frozen_store_growth!(p::_FrozenPolicy{T}, idx, oldidx, side, oldblock,
        newblock, growth, diagnostics) where T
    old = p.cuts[oldidx]
    firstnode = old.ledger === nothing
    firstnode && (p.generation += 1)
    generation = firstnode ? p.generation : old.ledger.generation
    ledger = firstnode || size(growth.Cup, 2) + size(growth.Cdn, 2) > 0 ?
        _FrozenLedger(generation, firstnode ? nothing : old.ledger,
            growth.Cup, growth.Cdn) : old.ledger
    ra, nold = newblock.ra, length(oldblock.ra)
    A = side === :left ? oldblock.phi' * newblock.phi[1:nold, :] :
        oldblock.phi' * newblock.phi[end - nold + 1:end, :]
    Pup, Pdn = growth.Cup * growth.Cup', growth.Cdn * growth.Cdn'
    jP, GPup, GPdn = _frozen_fields(p, growth.Cup, growth.Cdn)
    Cfu, Cfd = _ledger_columns(old, size(p.V, 1), :up),
        _ledger_columns(old, size(p.V, 1), :dn)
    _, Gfup, Gfdn = _frozen_fields(p, Cfu, Cfd)
    heffup, heffdn = p.Hup + Gfup, p.Hdn + Gfdn
    e = old.energy + tr(Pup * heffup) + tr(Pdn * heffdn) +
        0.5tr(Pup * GPup) + 0.5tr(Pdn * GPdn)
    kup = A' * old.kup * A +
        newblock.phi' * (p.V[ra, ra] .* Pup[ra, ra]) * newblock.phi
    kdn = A' * old.kdn * A +
        newblock.phi' * (p.V[ra, ra] .* Pdn[ra, ra]) * newblock.phi
    p.max_live = max(p.max_live, _live_generations(p, ledger))
    cut = _FrozenCut(ledger, old.nup + size(growth.Cup, 2),
        old.ndn + size(growth.Cdn, 2), old.j + jP, kup, kdn, e)
    p.cuts[idx] = cut
    if _detailed(diagnostics)
        Cu, Cd = _ledger_columns(cut, size(p.V, 1), :up),
            _ledger_columns(cut, size(p.V, 1), :dn)
        Du, Dd = Cu * Cu', Cd * Cd'
        jf = p.V * (diag(Du) + diag(Dd))
        ef = tr(Du * p.Hup) + tr(Dd * p.Hdn) +
            0.5dot(diag(Du) + diag(Dd), jf) -
            0.5sum(p.V .* (Du .* Du + Dd .* Dd))
        diagnostics.frozen_recurrence_error = max(diagnostics.frozen_recurrence_error,
            norm(cut.j - jf), norm(cut.kup - newblock.phi' *
                (p.V[ra, ra] .* Du[ra, ra]) * newblock.phi),
            norm(cut.kdn - newblock.phi' *
                (p.V[ra, ra] .* Dd[ra, ra]) * newblock.phi), abs(cut.energy - ef))
    end
end

function _window_basis(N, left, cra, right)
    B = zeros(promote_type(eltype(left.phi), eltype(right.phi)), N,
        left.m + length(cra) + right.m)
    B[left.ra, 1:left.m] = left.phi
    B[cra, left.m + 1:left.m + length(cra)] =
        Matrix{eltype(B)}(I, length(cra), length(cra))
    B[right.ra, end - right.m + 1:end] = right.phi
    B
end

function _allowed_basis(Cf, B)
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

function _frozen_effective(p, left, right, B)
    L, R = p.cuts[left], p.cuts[right]
    N = size(B, 1)
    CLu, CLd = _ledger_columns(L, N, :up), _ledger_columns(L, N, :dn)
    CRu, CRd = _ledger_columns(R, N, :up), _ledger_columns(R, N, :dn)
    j = L.j + R.j
    Gup = B' * (j .* B)
    Gdn = copy(Gup)
    ml, mr = size(L.kup, 1), size(R.kup, 1)
    Gup[1:ml, 1:ml] -= L.kup
    Gdn[1:ml, 1:ml] -= L.kdn
    Gup[end - mr + 1:end, end - mr + 1:end] -= R.kup
    Gdn[end - mr + 1:end, end - mr + 1:end] -= R.kdn
    qL = diag(CLu * CLu' + CLd * CLd')
    qR = diag(CRu * CRu' + CRd * CRd')
    ef = L.energy + R.energy + 0.5(dot(qL, R.j) + dot(qR, L.j))
    Gup, Gdn, ef, hcat(CLu, CRu), hcat(CLd, CRd)
end

function _frozen_local!(p, li, ri, left, cra, right, Hup, Hdn, win, lambda,
        scf_cutoff, relative, diagnostics)
    B = _window_basis(size(p.V, 1), left, cra, right)
    Gup, Gdn, ef, Cfu, Cfd = _frozen_effective(p, li, ri, B)
    Lu, Ld = _ledger_columns(p.cuts[li], size(B, 1), :up),
        _ledger_columns(p.cuts[li], size(B, 1), :dn)
    Ru, Rd = _ledger_columns(p.cuts[ri], size(B, 1), :up),
        _ledger_columns(p.cuts[ri], size(B, 1), :dn)
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
        for s = 1:4
            rhoup, rhodn = Zup * rup * Zup', Zdn * rdn * Zdn'
            Fup, Fdn = copy(Heffup), copy(Heffdn)
            if p.restricted
                vee_add_fock_r!(Fup, rhoup, win)
                nup > 0 && (Uup = eigsym(Zup' * Fup * Zup, diagnostics)[2][:, 1:nup])
                psiup = Zup * Uup
                rup = (1 - lambda) * rup + lambda * (Uup * Uup')
                rhoup = Zup * rup * Zup'
                Fup = copy(Heffup)
                vee_add_fock_r!(Fup, rhoup, win)
                energy = ef + tr(rhoup * (Fup + Heffup))
                psidn = psiup
            else
                vee_add_fock!(Fup, Fdn, rhoup, rhodn, win)
                nup > 0 && (Uup = eigsym(Zup' * Fup * Zup, diagnostics)[2][:, 1:nup])
                ndn > 0 && (Udn = eigsym(Zdn' * Fdn * Zdn, diagnostics)[2][:, 1:ndn])
                psiup, psidn = Zup * Uup, Zdn * Udn
                rup = (1 - lambda) * rup + lambda * (Uup * Uup')
                rdn = (1 - lambda) * rdn + lambda * (Udn * Udn')
                rhoup, rhodn = Zup * rup * Zup', Zdn * rdn * Zdn'
                Fup, Fdn = copy(Heffup), copy(Heffdn)
                vee_add_fock!(Fup, Fdn, rhoup, rhodn, win)
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

function _full_focks(p)
    Dup, Ddn = p.psiup * p.psiup', p.psidn * p.psidn'
    j = p.V * (diag(Dup) + diag(Ddn))
    p.Hup + Diagonal(j) - p.V .* Dup, p.Hdn + Diagonal(j) - p.V .* Ddn
end
function _frozen_full_energy(p)
    Fup, Fdn = _full_focks(p)
    Dup, Ddn = p.psiup * p.psiup', p.psidn * p.psidn'
    0.5tr(Dup * (Fup + p.Hup)) + 0.5tr(Ddn * (Fdn + p.Hdn))
end

function _aufbau_inverted(F, C)
    n, N = size(C, 2), size(C, 1)
    (n == 0 || n == N) && return false
    U = svd(C; full = true).U[:, n + 1:end]
    eo = eigmax(Symmetric(C' * F * C))
    ev = eigmin(Symmetric(U' * F * U))
    g = N * eps(eltype(F)) / (1 - N * eps(eltype(F)))
    eo - ev > 128g * max(1, norm(F, Inf))
end

function _frozen_audit!(p, prospective, diagnostics)
    Fup, Fdn = _full_focks(p)
    N, T = size(p.V, 1), eltype(p.V)
    g = N * eps(T) / (1 - N * eps(T))
    Cfu, Cfd = _all_frozen(p, :up), _all_frozen(p, :dn)
    staleup = !_frozen_valid(p.psiup, Cfu, Cfu[:, 1:0])
    staledn = !p.restricted && !_frozen_valid(p.psidn, Cfd, Cfd[:, 1:0])
    if staleup || staledn
        _protect_stale!(p, staleup ? Cfu : Cfu[:, 1:0],
            staledn ? Cfd : Cfd[:, 1:0], diagnostics)
        return true, false
    end
    rup, rdn = zeros(T, N, 0), zeros(T, N, 0)
    for (F, C, Cf, spin) in ((Fup, p.psiup, Cfu, :up),
            (Fdn, p.psidn, Cfd, :dn))
        tau = 128g * max(1, norm(F, Inf))
        for c in eachcol(Cf)
            r = F * c - C * (C' * (F * c))
            if norm(r) > tau
                spin === :up ? (rup = hcat(rup, c)) : (rdn = hcat(rdn, c))
            end
        end
    end
    invup = invdn = false
    if prospective
        timed = @timed begin
            invup = _aufbau_inverted(Fup, p.psiup)
            invdn = p.restricted ? invup : _aufbau_inverted(Fdn, p.psidn)
        end
        if diagnostics !== nothing
            diagnostics.aufbau_calls += 1
            diagnostics.aufbau_seconds += timed.time
            diagnostics.aufbau_bytes += timed.bytes
        end
        invup && size(Cfu, 2) > 0 && (rup = hcat(rup, Cfu))
        invdn && size(Cfd, 2) > 0 && (rdn = hcat(rdn, Cfd))
    end
    rebuild = size(rup, 2) + size(rdn, 2) > 0
    if rebuild
        p.protectedup = _orthcols(hcat(p.protectedup, rup))
        p.protecteddn = _orthcols(hcat(p.protecteddn, rdn))
        p.rebuilt = true
        diagnostics === nothing || (diagnostics.frozen_rebuilds += 1)
    end
    diagnostics === nothing ||
        (diagnostics.frozen_max_generations = p.max_live;
         diagnostics.frozen_live_generations = _live_generations(p))
    unresolved = prospective && p.rebuilt &&
        ((invup && size(Cfu, 2) == 0) || (invdn && size(Cfd, 2) == 0))
    rebuild, prospective && !rebuild && !unresolved
end
