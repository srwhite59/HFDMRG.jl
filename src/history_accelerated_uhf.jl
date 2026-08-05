struct _HistoryUHFModel
    Hup::Matrix{Float64}
    Hdn::Matrix{Float64}
    G::Array{Float64,4}
end
_history_gram(C) = isempty(C) ? 0.0 : norm(C' * C - I)
_history_overlap(C, trial) = isempty(C) ? 1.0 : minimum(svdvals(C' * trial))
function _history_uhf_inputs(Hup, Hdn, V, Cup, Cdn, policy)
    all(A -> A isa Matrix{Float64}, (Hup, Hdn, V, Cup, Cdn)) ||
        error("history UHF requires Matrix{Float64} inputs")
    N, nup, ndn = size(Hup, 1), size(Cup, 2), size(Cdn, 2)
    size(Hup) == size(Hdn) == size(V) == (N, N) && size(Cup, 1) == size(Cdn, 1) == N &&
        0 < nup + ndn && nup < N && ndn < N || error("inconsistent history UHF sizes")
    all(A -> all(isfinite, A), (Hup, Hdn, V, Cup, Cdn)) ||
        error("nonfinite history UHF input")
    for A in (Hup, Hdn, V)
        norm(A - A', Inf) <= 64eps(Float64) * N * max(1, norm(A, Inf)) ||
            error("history UHF matrices must be symmetric")
    end
    for C in (Cup, Cdn)
        _history_gram(C) <= 32eps(Float64) * sqrt(length(C)) ||
            error("history UHF orbitals must be orthonormal")
    end
    p = policy
    p.bootstrap_sweeps > 0 && length(p.captures) >= 2 && issorted(p.captures) &&
        allunique(p.captures) && last(p.captures) == p.bootstrap_sweeps &&
        first(p.captures) > 0 && p.fifo_length >= length(p.captures) &&
        p.cleanup_sweeps > 0 && p.max_cycles >= 0 && p.target >= 0 &&
        p.diis_iterations > 0 && p.diis_memory > 1 && 0 <= p.minimum_overlap <= 1 ||
        error("invalid history UHF policy")
    nothing
end
function _history_uhf_physical(Hup, Hdn, V, Cup, Cdn)
    Dup, Ddn = Cup * Cup', Cdn * Cdn'
    direct = V * (diag(Dup) + diag(Ddn))
    Fup, Fdn = Hup + Diagonal(direct) - V .* Dup, Hdn + Diagonal(direct) - V .* Ddn
    FupC, FdnC = Fup * Cup, Fdn * Cdn
    Rup, Rdn = FupC - Cup * (Cup' * FupC), FdnC - Cdn * (Cdn' * FdnC)
    rup, rdn = norm(Rup), norm(Rdn)
    sz = (size(Cup, 2) - size(Cdn, 2)) / 2
    (; energy = 0.5dot(Dup, Fup + Hup) + 0.5dot(Ddn, Fdn + Hdn),
       residual_up = rup, residual_dn = rdn, residual = max(rup, rdn),
       residual_rms = hypot(rup, rdn) / sqrt(2), gram_up = _history_gram(Cup), gram_dn = _history_gram(Cdn),
       s2 = sz^2 + (size(Cup, 2) + size(Cdn, 2)) / 2 - norm(Cup' * Cdn)^2)
end
function _history_uhf_basis(queue)
    Cup, Cdn = queue[end]
    N, nup, ndn = size(Cup, 1), size(Cup, 2), size(Cdn, 2)
    anchor_svd = svd(hcat(Cup, Cdn); full = false)
    anchor_cut = 64eps(Float64) * max(N, nup + ndn) * max(1, anchor_svd.S[1])
    nu = count(>(anchor_cut), anchor_svd.S)
    A = anchor_svd.U[:, 1:nu]
    parts = Matrix{Float64}[]
    for (up, dn) in queue[1:end-1]
        nup > 0 && push!(parts, up)
        ndn > 0 && push!(parts, dn)
    end
    q, h = (nup > 0) + (ndn > 0), length(queue) - 1
    Z = hcat(parts...) / sqrt(q * h)
    Z -= A * (A' * Z); Z -= A * (A' * Z)
    factor = svd(Z; full = false)
    values = factor.S .^ 2
    numerical_cut = isempty(factor.S) ? 0.0 :
        64eps(Float64) * max(N, q * h * max(nup, ndn)) * max(1, factor.S[1])
    knum = count(>(numerical_cut), factor.S)
    kabs = min(knum, count(>(2e-11), values))
    total, ktail = sum(values), 0
    if total > 0
        for k = 0:length(values)
            tail = k == length(values) ? 0.0 : sum(@view values[k + 1:end])
            sqrt(tail / total) <= 1e-4 && (ktail = k; break)
        end
    end
    k = min(knum, max(kabs, ktail))
    X = k == 0 ? copy(A) : hcat(A, factor.U[:, 1:k])
    Q = Matrix(qr(X).Q)[:, 1:size(X, 2)]
    gate = 32eps(Float64) * sqrt(N * size(Q, 2))
    cup_error, cdn_error = norm(Cup - Q * (Q' * Cup)), norm(Cdn - Q * (Q' * Cdn))
    tail = total == 0 ? 0.0 : sqrt(sum(@view values[k + 1:end]) / total)
    (; Q, rank = size(Q, 2), extra = k, knum, kabs, ktail, tail,
       valid = cup_error <= gate && cdn_error <= gate && norm(Q' * Q - I) <= gate)
end
function _history_uhf_model(Hup, Hdn, V, Q)
    N, r = size(Q)
    _history_workspace(N, r) || return nothing
    Hqup = Q' * Hup * Q; Hqdn = Hdn === Hup ? Hqup : Q' * Hdn * Q
    P = reshape([Q[p, a] * Q[p, b] for p = 1:N, a = 1:r, b = 1:r], N, r^2); G = reshape(P' * (V * P), r, r, r, r)
    _HistoryUHFModel(Hqup, Hqdn, G)
end
function _history_uhf_model_metrics(model, Cup, Cdn)
    Dup, Ddn, r = Cup * Cup', Cdn * Cdn', size(Cup, 1)
    Fup, Fdn = copy(model.Hup), copy(model.Hdn)
    for b = 1:r, a = 1:r, d = 1:r, c = 1:r
        direct = model.G[a, b, c, d] * (Dup[c, d] + Ddn[c, d])
        Fup[a, b] += direct - model.G[a, c, d, b] * Dup[c, d]
        Fdn[a, b] += direct - model.G[a, c, d, b] * Ddn[c, d]
    end
    Fup = (Fup + Fup') / 2; Fdn = (Fdn + Fdn') / 2
    FupC, FdnC = Fup * Cup, Fdn * Cdn
    Rup, Rdn = FupC - Cup * (Cup' * FupC), FdnC - Cdn * (Cdn' * FdnC)
    Kup, Kdn = Fup * Dup - Dup * Fup, Fdn * Ddn - Ddn * Fdn
    (; Fup, Fdn, energy = 0.5dot(Dup, Fup + model.Hup) + 0.5dot(Ddn, Fdn + model.Hdn),
       residual_up = norm(Rup), residual_dn = norm(Rdn), Kup, Kdn,
       kresidual = hypot(norm(Kup), norm(Kdn)) / sqrt(2))
end
_history_uhf_converged(m) = norm(m.Kup) <= 1e-10 * max(1, norm(m.Fup)) &&
    norm(m.Kdn) <= 1e-10 * max(1, norm(m.Fdn))
_history_occupied(F, n) = n == 0 ? zeros(size(F, 1), 0) : eigen(Symmetric(F)).vectors[:, 1:n]
function _history_uhf_diis(model, cup0, cdn0, policy)
    cup, cdn = copy(cup0), copy(cdn0)
    first = _history_uhf_model_metrics(model, cup, cdn)
    bestup, bestdn, best_energy, best_k = copy(cup), copy(cdn), first.energy, first.kresidual
    focks, errors = Matrix{Float64}[], Matrix{Float64}[]
    consecutive_failures = pulay_rank = 0
    for iteration = 1:policy.diis_iterations
        metrics = _history_uhf_model_metrics(model, cup, cdn)
        finite = all(isfinite, (metrics.energy, metrics.kresidual)) &&
            all(isfinite, metrics.Fup) && all(isfinite, metrics.Fdn)
        finite || return (; cup = bestup, cdn = bestdn, status = :numerical_failure,
            iterations = iteration, failures = consecutive_failures, pulay_rank)
        tau = 128eps(Float64) * max(1, abs(best_energy), abs(metrics.energy))
        if metrics.energy < best_energy - tau ||
                (abs(metrics.energy - best_energy) <= tau && metrics.kresidual < best_k)
            bestup, bestdn = copy(cup), copy(cdn)
            best_energy, best_k = metrics.energy, metrics.kresidual
        end
        _history_uhf_converged(metrics) && return (; cup = bestup, cdn = bestdn,
            status = :converged_best, iterations = iteration, failures = 0, pulay_rank)
        push!(focks, hcat(metrics.Fup, metrics.Fdn))
        push!(errors, hcat(metrics.Kup, metrics.Kdn))
        length(focks) > policy.diis_memory && (popfirst!(focks); popfirst!(errors))
        F, status, pulay_rank, _, affine_sum = _history_pulay(focks, errors)
        consecutive_failures = status === :failed ? consecutive_failures + 1 : 0
        consecutive_failures == 2 && return (; cup = bestup, cdn = bestdn,
            status = :pulay_failure, iterations = iteration, failures = 2, pulay_rank)
        all(isfinite, F) || return (; cup = bestup, cdn = bestdn,
            status = :numerical_failure, iterations = iteration,
            failures = consecutive_failures, pulay_rank)
        r = size(cup0, 1)
        cup = _history_occupied(Matrix(@view F[:, 1:r]), size(cup0, 2))
        cdn = _history_occupied(Matrix(@view F[:, r + 1:2r]), size(cdn0, 2))
    end
    (; cup = bestup, cdn = bestdn, status = :cap, iterations = policy.diis_iterations,
       failures = consecutive_failures, pulay_rank)
end
function _history_uhf_acceptable(before, Cup, Cdn, trialup, trialdn, audit, status, overlap)
    finite = all(isfinite, trialup) && all(isfinite, trialdn) &&
        all(isfinite, (audit.energy, audit.residual_up, audit.residual_dn,
            audit.gram_up, audit.gram_dn, audit.s2))
    oup, odn = finite ? (_history_overlap(Cup, trialup),
        _history_overlap(Cdn, trialdn)) : (-Inf, -Inf)
    tau = 128eps(Float64) * max(1, abs(before.energy), abs(audit.energy))
    valid = finite && status ∉ (:pulay_failure, :numerical_failure) &&
        audit.gram_up <= 32eps(Float64) * sqrt(length(Cup)) &&
        audit.gram_dn <= 32eps(Float64) * sqrt(length(Cdn)) &&
        audit.energy <= before.energy + tau && oup >= overlap && odn >= overlap
    valid, oup, odn
end
function _history_uhf_segment(Hup, Hdn, V, Cup, Cdn, sweeps, captures; solver...)
    saved = Tuple{Matrix{Float64},Matrix{Float64}}[]
    observer = info -> begin
        info.sweep in captures && push!(saved, (copy(info.psiup), copy(info.psidn)))
        false
    end
    timed = @timed if Hup === Hdn
        solve_hfdmrg(Hup, V, Cup, Cdn; solver..., maxiter = sweeps, cutoff = 0.0, observer)
    else
        solve_hfdmrg(Hup, Hdn, V, Cup, Cdn; solver..., maxiter = sweeps, cutoff = 0.0, observer)
    end
    timed.value[1], timed.value[2], saved, timed
end
_solve_history_uhf(H, V, Cup, Cdn; kwargs...) = _solve_history_uhf(H, H, V, Cup, Cdn; kwargs...)
function _solve_history_uhf(Hup, Hdn, V, Cup0, Cdn0;
        _policy = _HistoryRHFPolicy(), nblockcenter = 1, blocksize = 200,
        block_partition = nothing, scf_cutoff = 1e-11,
        environment_cutoff = 1e-10, verbose = false)
    _history_uhf_inputs(Hup, Hdn, V, Cup0, Cdn0, _policy)
    solver = (; nblockcenter, blocksize, block_partition, scf_cutoff,
        environment_cutoff, verbose)
    Cup, Cdn, queue, bootstrap = _history_uhf_segment(Hup, Hdn, V, Cup0, Cdn0,
        _policy.bootstrap_sweeps, _policy.captures; solver...)
    length(queue) == length(_policy.captures) || error("incomplete UHF history capture")
    audit = @timed _history_uhf_physical(Hup, Hdn, V, Cup, Cdn)
    events = NamedTuple[]
    sweep_time, sweep_bytes = bootstrap.time, bootstrap.bytes
    basis_time = setup_time = diis_time = audit_time = 0.0
    basis_bytes = setup_bytes = diis_bytes = 0; audit_bytes = audit.bytes
    audit_time += audit.time
    accepted = rejected = noops = resources = cycles = 0
    while audit.value.residual > _policy.target && cycles < _policy.max_cycles
        cycles += 1
        incoming_energy = audit.value.energy
        basis_timed = @timed _history_uhf_basis(queue); basis = basis_timed.value
        basis_time += basis_timed.time; basis_bytes += basis_timed.bytes
        status, diis_status, trialup, trialdn = :no_op, :not_run, Cup, Cdn
        proposal = (; energy = NaN, residual_up = NaN, residual_dn = NaN,
            gram_up = NaN, gram_dn = NaN, s2 = NaN)
        overlap_up = overlap_dn = NaN
        diis_iterations = diis_failures = diis_rank = 0
        if basis.extra == 0
            noops += 1
        elseif !basis.valid
            status = :rejected; rejected += 1
        else
            setup_timed = @timed _history_uhf_model(Hup, Hdn, V, basis.Q)
            setup_time += setup_timed.time; setup_bytes += setup_timed.bytes
            if setup_timed.value === nothing
                status = :resource_fallback; resources += 1
            else
                diis_timed = @timed _history_uhf_diis(setup_timed.value,
                    basis.Q' * Cup, basis.Q' * Cdn, _policy)
                diis_time += diis_timed.time; diis_bytes += diis_timed.bytes
                diis = diis_timed.value
                diis_status, diis_iterations, diis_failures, diis_rank =
                    diis.status, diis.iterations, diis.failures, diis.pulay_rank
                candidate_up, candidate_dn = basis.Q * diis.cup, basis.Q * diis.cdn
                candidate_audit = @timed _history_uhf_physical(
                    Hup, Hdn, V, candidate_up, candidate_dn)
                audit_time += candidate_audit.time; audit_bytes += candidate_audit.bytes
                proposal = candidate_audit.value
                valid, overlap_up, overlap_dn = _history_uhf_acceptable(audit.value,
                    Cup, Cdn, candidate_up, candidate_dn, proposal, diis_status,
                    _policy.minimum_overlap)
                if valid
                    status = :accepted; accepted += 1
                    trialup, trialdn = candidate_up, candidate_dn
                else
                    status = :rejected; rejected += 1
                end
            end
        end
        Cup, Cdn, _, cleanup = _history_uhf_segment(Hup, Hdn, V, trialup, trialdn,
            _policy.cleanup_sweeps, (); solver...)
        sweep_time += cleanup.time; sweep_bytes += cleanup.bytes
        push!(queue, (Cup, Cdn)); length(queue) > _policy.fifo_length && popfirst!(queue)
        audit = @timed _history_uhf_physical(Hup, Hdn, V, Cup, Cdn)
        audit_time += audit.time; audit_bytes += audit.bytes
        push!(events, (; cycle = cycles, status, rank = basis.rank,
            extra_rank = basis.extra, tail = basis.tail, diis_status, diis_iterations,
            diis_failures, diis_rank, proposal_energy = proposal.energy,
            proposal_residual_up = proposal.residual_up,
            proposal_residual_dn = proposal.residual_dn, proposal_gram_up = proposal.gram_up, proposal_gram_dn = proposal.gram_dn, proposal_s2 = proposal.s2, overlap_up, overlap_dn,
            incoming_energy, energy = audit.value.energy, residual = audit.value.residual))
    end
    (; Cup, Cdn, energy = audit.value.energy, residual = audit.value.residual,
       residual_up = audit.value.residual_up, residual_dn = audit.value.residual_dn,
       residual_rms = audit.value.residual_rms, gram_up = audit.value.gram_up,
       gram_dn = audit.value.gram_dn, s2 = audit.value.s2,
       target_reached = audit.value.residual <= _policy.target, terminal_reason = audit.value.residual <= _policy.target ? :target : :cycle_cap,
       total_sweeps = _policy.bootstrap_sweeps + cycles * _policy.cleanup_sweeps, cycles,
       accepted, rejected, noops, resource_fallbacks = resources, events,
       times = (; sweeps = sweep_time, basis = basis_time, setup = setup_time, diis = diis_time, audit = audit_time),
       bytes = (; sweeps = sweep_bytes, basis = basis_bytes, setup = setup_bytes, diis = diis_bytes, audit = audit_bytes))
end
