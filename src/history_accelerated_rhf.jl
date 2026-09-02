Base.@kwdef struct _HistoryRHFPolicy
    bootstrap_sweeps::Int = 6
    captures::Tuple{Vararg{Int}} = (2, 4, 6)
    fifo_length::Int = 6
    cleanup_sweeps::Int = 2
    max_cycles::Int = 25
    target::Float64 = 1e-6
    diis_iterations::Int = 40
    diis_memory::Int = 7
    minimum_overlap::Float64 = 0.99
end

struct _HistoryRHFModel
    H::Matrix{Float64}
    G::Array{Float64,4}
end

function _history_inputs(H, V, C, policy)
    H isa Matrix{Float64} && V isa Matrix{Float64} && C isa Matrix{Float64} ||
        error("history RHF requires Matrix{Float64} H, V, and C0")
    N, n = size(C)
    size(H) == size(V) == (N, N) && 0 < n < N || error("inconsistent history RHF sizes")
    all(isfinite, H) && all(isfinite, V) && all(isfinite, C) || error("nonfinite history RHF input")
    for A in (H, V)
        norm(A - A', Inf) <= 64eps(Float64) * N * max(1, norm(A, Inf)) ||
            error("history RHF matrices must be symmetric")
    end
    norm(C' * C - I) <= 32eps(Float64) * sqrt(N * n) ||
        error("history RHF orbitals must be orthonormal")
    p = policy
    p.bootstrap_sweeps > 0 && length(p.captures) >= 2 &&
        issorted(p.captures) && allunique(p.captures) && last(p.captures) == p.bootstrap_sweeps &&
        first(p.captures) > 0 && p.fifo_length >= length(p.captures) &&
        p.cleanup_sweeps > 0 && p.max_cycles >= 0 && p.target >= 0 &&
        p.diis_iterations > 0 && p.diis_memory > 1 && 0 <= p.minimum_overlap <= 1 ||
        error("invalid history RHF policy")
    nothing
end

function _history_physical(H, V, C)
    D = C * C'
    d = diag(D)
    J = V * d
    F = H + 2Diagonal(J) - V .* D
    FC = F * C
    R = FC - C * (C' * FC)
    (; energy = 2dot(D, H) + 2dot(d, J) - sum(V .* D .* D),
       residual = norm(R), gram = norm(C' * C - I))
end

function _history_acceptable(before, C, candidate, audit, status, minimum_overlap)
    finite = all(isfinite, candidate) &&
        all(isfinite, (audit.energy, audit.residual, audit.gram))
    overlap = finite ? minimum(svdvals(C' * candidate)) : -Inf
    tau = 128eps(Float64) * max(1, abs(before.energy), abs(audit.energy))
    valid = finite && status ∉ (:pulay_failure, :numerical_failure) &&
        audit.gram <= 32eps(Float64) * sqrt(length(C)) &&
        audit.energy <= before.energy + tau && overlap >= minimum_overlap
    valid, overlap
end

function _history_basis(queue)
    C, h = queue[end], length(queue) - 1
    Z = hcat(queue[1:end-1]...) / sqrt(h)
    Z -= C * (C' * Z); Z -= C * (C' * Z)
    factor = svd(Z; full = false)
    values = factor.S .^ 2
    threshold = isempty(factor.S) ? 0.0 :
        64eps(Float64) * max(size(C, 1), h * size(C, 2)) * max(1, factor.S[1])
    knum = count(>(threshold), factor.S)
    kabs = min(knum, count(>(2e-11), values))
    total, ktail = sum(values), 0
    if total > 0
        for k = 0:length(values)
            tail = k == length(values) ? 0.0 : sum(@view values[k + 1:end])
            sqrt(tail / total) <= 1e-4 && (ktail = k; break)
        end
    end
    k = min(knum, max(kabs, ktail))
    X = k == 0 ? copy(C) : hcat(C, factor.U[:, 1:k])
    Q = Matrix(qr(X).Q)[:, 1:size(X, 2)]
    gate = 32eps(Float64) * sqrt(size(C, 1) * size(Q, 2))
    containment, orthogonality = norm(C - Q * (Q' * C)), norm(Q' * Q - I)
    tail = total == 0 ? 0.0 : sqrt(sum(@view values[k + 1:end]) / total)
    (; Q, rank = size(Q, 2), extra = k, knum, kabs, ktail,
       tail, valid = containment <= gate && orthogonality <= gate)
end

function _history_workspace(N, r)
    try
        r2 = Base.checked_mul(r, r)
        Base.checked_add(Base.checked_mul(Base.checked_mul(2, N), r2),
            Base.checked_mul(r2, r2)) <= Base.checked_mul(N, N)
    catch error
        error isa OverflowError || rethrow()
        false
    end
end

function _history_model(H, V, Q)
    N, r = size(Q)
    _history_workspace(N, r) || return nothing
    Hq = Q' * H * Q
    P = reshape([Q[p, a] * Q[p, b] for p = 1:N, a = 1:r, b = 1:r], N, r^2)
    G = reshape(P' * (V * P), r, r, r, r)
    _HistoryRHFModel(Hq, G)
end

function _history_model_metrics(model, C)
    D, r = C * C', size(C, 1)
    F = model.H + 2reshape(reshape(model.G, r^2, r^2) * vec(D), r, r)
    for b = 1:r, a = 1:r, d = 1:r, c = 1:r
        F[a, b] -= model.G[a, c, d, b] * D[c, d]
    end
    F = (F + F') / 2
    FC = F * C
    R, K = FC - C * (C' * FC), F * D - D * F
    (; F, energy = dot(D, F + model.H), residual = norm(R),
       commutator = K, kresidual = norm(K))
end

function _history_pulay(focks, errors)
    m = length(focks)
    m < 2 && return focks[end], :unavailable, 1, 1.0, 1.0
    B = [dot(errors[i], errors[j]) for i = 1:m, j = 1:m]
    scale = maximum(abs, B)
    isfinite(scale) && scale > 0 || return focks[end], :failed, 0, Inf, Inf
    A = zeros(m + 1, m + 1); A[1:m, 1:m] .= B ./ scale
    A[1:m, end] .= -1; A[end, 1:m] .= -1
    rhs = zeros(m + 1); rhs[end] = -1
    factor = svd(A)
    rank = count(>(128eps(Float64) * (m + 1) * factor.S[1]), factor.S)
    rank > 0 || return focks[end], :failed, 0, Inf, Inf
    solution = factor.V[:, 1:rank] * ((factor.U[:, 1:rank]' * rhs) ./ factor.S[1:rank])
    coefficients, s = solution[1:m], sum(solution[1:m])
    cnorm = norm(coefficients)
    valid = all(isfinite, coefficients) && isfinite(cnorm) && cnorm <= 1e6 &&
        abs(s - 1) <= 1024eps(Float64) * max(1, norm(coefficients, 1))
    valid || return focks[end], :failed, rank, cnorm, s
    coefficients ./= s
    F = zeros(size(focks[end]))
    for i = 1:m
        F .+= coefficients[i] .* focks[i]
    end
    F, :used, rank, cnorm, sum(coefficients)
end

function _history_diis(model, c0, policy)
    c = copy(c0)
    first = _history_model_metrics(model, c)
    best_c, best_energy, best_k = copy(c), first.energy, first.kresidual
    focks, errors = Matrix{Float64}[], Matrix{Float64}[]
    consecutive_failures = pulay_rank = 0
    for iteration = 1:policy.diis_iterations
        metrics = _history_model_metrics(model, c)
        finite = all(isfinite, (metrics.energy, metrics.kresidual)) && all(isfinite, metrics.F)
        finite || return (; c = best_c, status = :numerical_failure,
            iterations = iteration, failures = consecutive_failures, pulay_rank)
        tau = 128eps(Float64) * max(1, abs(best_energy), abs(metrics.energy))
        if metrics.energy < best_energy - tau ||
                (abs(metrics.energy - best_energy) <= tau && metrics.kresidual < best_k)
            best_c, best_energy, best_k = copy(c), metrics.energy, metrics.kresidual
        end
        metrics.kresidual <= 1e-10 * max(1, norm(metrics.F)) &&
            return (; c = best_c, status = :converged_best, iterations = iteration,
                failures = 0, pulay_rank)
        push!(focks, metrics.F); push!(errors, metrics.commutator)
        length(focks) > policy.diis_memory && (popfirst!(focks); popfirst!(errors))
        F, status, pulay_rank, _, _ = _history_pulay(focks, errors)
        consecutive_failures = status === :failed ? consecutive_failures + 1 : 0
        consecutive_failures == 2 && return (; c = best_c,
            status = :pulay_failure, iterations = iteration, failures = 2, pulay_rank)
        all(isfinite, F) || return (; c = best_c, status = :numerical_failure,
            iterations = iteration, failures = consecutive_failures, pulay_rank)
        c = eigen(Symmetric(F)).vectors[:, 1:size(c0, 2)]
    end
    (; c = best_c, status = :cap, iterations = policy.diis_iterations,
       failures = consecutive_failures, pulay_rank)
end

function _history_segment(H, V, C, sweeps, captures; solver...)
    saved = Matrix{Float64}[]
    observer = info -> begin
        info.sweep in captures && push!(saved, copy(info.psiup))
        false
    end
    timed = @timed solve_hfdmrg(H, V, C; solver..., maxiter = sweeps,
        cutoff = 0.0, observer)
    timed.value[1], saved, timed
end

function _solve_history_rhf(H, V, C0; _policy = _HistoryRHFPolicy(),
        nblockcenter = 1, blocksize = 200, block_partition = nothing,
        scf_cutoff = 1e-11, environment_cutoff = 1e-10, verbose = false)
    _removed_rhf_route("history-accelerated")
    _history_inputs(H, V, C0, _policy)
    start_ns = time_ns()
    solver = (; nblockcenter, blocksize, block_partition, scf_cutoff,
        environment_cutoff, verbose)
    C, queue, bootstrap = _history_segment(
        H, V, C0, _policy.bootstrap_sweeps, _policy.captures; solver...)
    length(queue) == length(_policy.captures) || error("incomplete history capture")
    audit = @timed _history_physical(H, V, C)
    events = NamedTuple[]
    sweep_time, sweep_bytes = bootstrap.time, bootstrap.bytes
    basis_time = setup_time = diis_time = audit_time = 0.0
    basis_bytes = setup_bytes = diis_bytes = 0
    audit_bytes = audit.bytes
    audit_time += audit.time
    accepted = rejected = noops = resources = 0
    cycles = 0
    while audit.value.residual > _policy.target && cycles < _policy.max_cycles
        cycles += 1
        incoming_energy = audit.value.energy
        basis_timed = @timed _history_basis(queue)
        basis = basis_timed.value
        basis_time += basis_timed.time; basis_bytes += basis_timed.bytes
        status, diis_status, trial = :no_op, :not_run, C
        proposal = (; energy = NaN, residual = NaN, gram = NaN)
        overlap = NaN; diis_iterations = diis_failures = diis_rank = 0
        if basis.extra == 0
            noops += 1
        elseif !basis.valid
            status = :rejected; rejected += 1
        else
            setup_timed = @timed _history_model(H, V, basis.Q)
            setup_time += setup_timed.time; setup_bytes += setup_timed.bytes
            if setup_timed.value === nothing
                status = :resource_fallback; resources += 1
            else
                diis_timed = @timed _history_diis(
                    setup_timed.value, basis.Q' * C, _policy)
                diis_time += diis_timed.time; diis_bytes += diis_timed.bytes
                diis = diis_timed.value
                diis_status, diis_iterations, diis_failures, diis_rank =
                    diis.status, diis.iterations, diis.failures, diis.pulay_rank
                candidate = basis.Q * diis.c
                candidate_audit = @timed _history_physical(H, V, candidate)
                audit_time += candidate_audit.time; audit_bytes += candidate_audit.bytes
                before = audit.value
                proposal = candidate_audit.value
                valid, overlap = _history_acceptable(
                    before, C, candidate, proposal, diis_status, _policy.minimum_overlap)
                if valid
                    status = :accepted; accepted += 1; trial = candidate
                else
                    status = :rejected; rejected += 1
                end
            end
        end
        C, _, cleanup = _history_segment(
            H, V, trial, _policy.cleanup_sweeps, (); solver...)
        sweep_time += cleanup.time; sweep_bytes += cleanup.bytes
        push!(queue, C); length(queue) > _policy.fifo_length && popfirst!(queue)
        audit = @timed _history_physical(H, V, C)
        audit_time += audit.time; audit_bytes += audit.bytes
        push!(events, (; cycle = cycles, status, rank = basis.rank,
            extra_rank = basis.extra, tail = basis.tail, diis_status,
            diis_iterations, diis_failures, diis_rank, proposal_energy = proposal.energy,
            proposal_residual = proposal.residual, proposal_gram = proposal.gram,
            overlap, incoming_energy, energy = audit.value.energy, residual = audit.value.residual,
            elapsed = (time_ns() - start_ns) * 1e-9))
    end
    (; C, energy = audit.value.energy, residual = audit.value.residual,
       target_reached = audit.value.residual <= _policy.target,
       terminal_reason = audit.value.residual <= _policy.target ? :target : :cycle_cap,
       total_sweeps = _policy.bootstrap_sweeps + cycles * _policy.cleanup_sweeps,
       cycles, accepted, rejected, noops, resource_fallbacks = resources, events,
       times = (; sweeps = sweep_time, basis = basis_time, setup = setup_time,
           diis = diis_time, audit = audit_time),
       bytes = (; sweeps = sweep_bytes, basis = basis_bytes, setup = setup_bytes,
           diis = diis_bytes, audit = audit_bytes))
end
