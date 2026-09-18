using Test
using LinearAlgebra
using Random
using HFDMRG

function dense_history_rhf(H, V, C)
    D = C * C'
    d = diag(D)
    J = V * d
    F = H + 2Diagonal(J) - V .* D
    FC = F * C
    (; F, FC, energy = 2dot(D, H) + 2dot(d, J) - sum(V .* D .* D),
       residual = norm(FC - C * (C' * FC)), gram = norm(C' * C - I))
end

function dense_history_uhf(Hup, Hdn, V, Cup, Cdn)
    Dup, Ddn = Cup * Cup', Cdn * Cdn'
    direct = V * (diag(Dup) + diag(Ddn))
    Fup = Hup + Diagonal(direct) - V .* Dup
    Fdn = Hdn + Diagonal(direct) - V .* Ddn
    FupC, FdnC = Fup * Cup, Fdn * Cdn
    Rup = FupC - Cup * (Cup' * FupC)
    Rdn = FdnC - Cdn * (Cdn' * FdnC)
    rup, rdn = norm(Rup), norm(Rdn)
    sz = (size(Cup, 2) - size(Cdn, 2)) / 2
    (; FupC, FdnC,
       energy = 0.5dot(Dup, Fup + Hup) + 0.5dot(Ddn, Fdn + Hdn),
       residual_up = rup, residual_dn = rdn, residual = max(rup, rdn),
       residual_rms = hypot(rup, rdn) / sqrt(2),
       K_alpha = sqrt(2) * rup, K_beta = sqrt(2) * rdn,
       K = sqrt(2) * hypot(rup, rdn), gram_up = norm(Cup' * Cup - I),
       gram_dn = norm(Cdn' * Cdn - I),
       s2 = sz^2 + (size(Cup, 2) + size(Cdn, 2)) / 2 - norm(Cup' * Cdn)^2)
end

@testset "History occupied-action physical audits" begin
    rng = MersenneTwister(91701)
    orthogonal(n, k) = Matrix(qr(randn(rng, n, k)).Q)[:, 1:k]
    N = 100
    A = randn(rng, N, N)
    H = Matrix(Symmetric(A))
    Hdn = H + Diagonal(range(-0.03, 0.02; length = N))
    B = randn(rng, N, N)
    V = 0.02Matrix(Symmetric(B))

    C = orthogonal(N, 3)
    dense = dense_history_rhf(H, V, C)
    audit = HFDMRG._history_physical(H, V, C)
    density = vec(sum(abs2, C; dims = 2))
    action, _ = HFDMRG._history_occupied_action(H, V, C, 2(V * density), C)
    @test isapprox(action, dense.FC; atol = 2e-13, rtol = 2e-14)
    @test isapprox(audit.energy, dense.energy; atol = 2e-12, rtol = 2e-14)
    @test isapprox(audit.residual, dense.residual; atol = 2e-13, rtol = 2e-14)
    @test audit.gram == dense.gram

    gauge = Matrix(qr(randn(rng, 3, 3)).Q)
    rotated = HFDMRG._history_physical(H, V, C * gauge)
    @test isapprox(rotated.energy, audit.energy; atol = 2e-12, rtol = 2e-14)
    @test isapprox(rotated.residual, audit.residual; atol = 2e-13, rtol = 2e-14)

    levels = collect(range(-2.0, 2.0; length = N))
    stationary = HFDMRG._history_physical(Matrix(Diagonal(levels)), zeros(N, N),
        Matrix(I, N, N)[:, 1:3])
    @test stationary.energy == 2sum(levels[1:3])
    @test stationary.residual == stationary.gram == 0

    Cup, Cdn = orthogonal(N, 3), orthogonal(N, 1)
    dense_u = dense_history_uhf(H, Hdn, V, Cup, Cdn)
    audit_u = HFDMRG._history_uhf_physical(H, Hdn, V, Cup, Cdn)
    for field in (:energy, :residual_up, :residual_dn, :residual, :residual_rms,
            :K_alpha, :K_beta, :K, :gram_up, :gram_dn, :s2)
        @test isapprox(getproperty(audit_u, field), getproperty(dense_u, field);
            atol = 3e-12, rtol = 3e-14)
    end
    total_density = vec(sum(abs2, Cup; dims = 2) + sum(abs2, Cdn; dims = 2))
    direct = V * total_density
    up_action, _ = HFDMRG._history_occupied_action(H, V, Cup, direct, Cup)
    dn_action, _ = HFDMRG._history_occupied_action(Hdn, V, Cdn, direct, Cdn)
    @test isapprox(up_action, dense_u.FupC; atol = 3e-13, rtol = 3e-14)
    @test isapprox(dn_action, dense_u.FdnC; atol = 3e-13, rtol = 3e-14)

    up_gauge = Matrix(qr(randn(rng, 3, 3)).Q)
    rotated_u = HFDMRG._history_uhf_physical(H, Hdn, V, Cup * up_gauge, -Cdn)
    @test isapprox(rotated_u.energy, audit_u.energy; atol = 3e-12, rtol = 3e-14)
    @test isapprox(rotated_u.residual_up, audit_u.residual_up;
        atol = 3e-13, rtol = 3e-14)
    @test isapprox(rotated_u.residual_dn, audit_u.residual_dn;
        atol = 3e-13, rtol = 3e-14)
    @test isapprox(rotated_u.s2, audit_u.s2; atol = 3e-13, rtol = 3e-14)

    empty_dn = zeros(N, 0)
    dense_empty = dense_history_uhf(H, Hdn, V, Cup, empty_dn)
    audit_empty = HFDMRG._history_uhf_physical(H, Hdn, V, Cup, empty_dn)
    @test audit_empty.residual_dn == audit_empty.gram_dn == audit_empty.K_beta == 0
    @test isapprox(audit_empty.energy, dense_empty.energy; atol = 3e-12, rtol = 3e-14)
    @test isapprox(audit_empty.residual_up, dense_empty.residual_up;
        atol = 3e-13, rtol = 3e-14)

    candidate = orthogonal(N, 3)
    old_candidate = dense_history_rhf(H, V, candidate)
    new_candidate = HFDMRG._history_physical(H, V, candidate)
    old_decision = HFDMRG._history_acceptable(dense, C, candidate, old_candidate,
        :converged_best, 0.0)
    new_decision = HFDMRG._history_acceptable(audit, C, candidate, new_candidate,
        :converged_best, 0.0)
    @test old_decision == new_decision

    response = orthogonal(N, 1)
    direction = reshape([1.0, 0.0, 0.0], 1, :)
    queue = [Matrix(qr(C + scale * response * direction).Q)[:, 1:3]
        for scale = (0.03, -0.02)]
    push!(queue, C)
    basis = HFDMRG._history_basis(queue)
    policy = HFDMRG._HistoryRHFPolicy(diis_iterations = 12, minimum_overlap = 0.0)
    proposal = basis.Q * HFDMRG._history_diis(
        HFDMRG._history_model(H, V, basis.Q), basis.Q' * C, policy).c
    old_proposal = dense_history_rhf(H, V, proposal)
    new_proposal = HFDMRG._history_physical(H, V, proposal)
    @test HFDMRG._history_acceptable(dense, C, proposal, old_proposal,
        :converged_best, 0.0) == HFDMRG._history_acceptable(
        audit, C, proposal, new_proposal, :converged_best, 0.0)

    up_response, dn_response = orthogonal(N, 1), orthogonal(N, 1)
    uqueue = [(Matrix(qr(Cup + scale * up_response * direction).Q)[:, 1:3],
        Matrix(qr(Cdn - scale * dn_response).Q)[:, 1:1])
        for scale = (0.03, -0.02)]
    push!(uqueue, (Cup, Cdn))
    ubasis = HFDMRG._history_uhf_basis(uqueue)
    udiis = HFDMRG._history_uhf_diis(
        HFDMRG._history_uhf_model(H, Hdn, V, ubasis.Q),
        ubasis.Q' * Cup, ubasis.Q' * Cdn, policy)
    proposal_up, proposal_dn = ubasis.Q * udiis.cup, ubasis.Q * udiis.cdn
    old_u_proposal = dense_history_uhf(H, Hdn, V, proposal_up, proposal_dn)
    new_u_proposal = HFDMRG._history_uhf_physical(H, Hdn, V,
        proposal_up, proposal_dn)
    @test HFDMRG._history_uhf_acceptable(dense_u, Cup, Cdn,
        proposal_up, proposal_dn, old_u_proposal, udiis.status, 0.0) ==
        HFDMRG._history_uhf_acceptable(audit_u, Cup, Cdn,
        proposal_up, proposal_dn, new_u_proposal, udiis.status, 0.0)

    dims = fill(2, 4)
    layout = HFDMRG.SliceLayout(dims)
    factors = [cat([Matrix(Symmetric(randn(rng, 2, 2))) for _ = 1:2]...;
        dims = 3) for _ = eachindex(dims)]
    sliced_v = zeros(2, 2, 2, 2, length(dims), length(dims))
    for s = eachindex(dims), t = eachindex(dims)
        sliced_v[:, :, :, :, s, t] = reshape(
            reshape(factors[s], 4, 2) * reshape(factors[t], 4, 2)', 2, 2, 2, 2)
    end
    backend = HFDMRG.SlicedBasisBackend(layout, sliced_v)
    Ns = sum(dims)
    As = randn(rng, Ns, Ns)
    Hs = Matrix(Symmetric(As))
    Hsdn = Hs + Diagonal(range(-0.01, 0.02; length = Ns))
    Cs, Cus, Cds = orthogonal(Ns, 2), orthogonal(Ns, 2), orthogonal(Ns, 1)
    Ds, Fs = Cs * Cs', copy(Hs)
    HFDMRG.sliced_add_fock_r!(Fs, Ds, sliced_v, layout)
    FCs = Fs * Cs
    sliced_dense = (; energy = dot(Ds, Fs + Hs),
        residual = norm(FCs - Cs * (Cs' * FCs)), gram = norm(Cs' * Cs - I))
    sliced_context = HFDMRG._SlicedFockAuditContext(Hs, backend, Cs, Cs, :up)
    sliced_context_action = HFDMRG._audit_action!(similar(Cs), sliced_context, Cs)
    prepared_j, prepared_offsets = HFDMRG._sliced_hartree_data(
        layout, sliced_v, Cs, Cs)
    @test sliced_context_action ≈ FCs
    @test sliced_context.actions == 1 && sliced_context.rhs == size(Cs, 2)
    @test prepared_j == sliced_context.Jpack && prepared_offsets == sliced_context.poffs
    @test sliced_context.diag ≈ diag(Fs)
    @test sliced_context.scale ≈ opnorm(Fs, Inf)
    @test sliced_context.setup_reason === :ok
    @test all(eachindex(dims)) do slice
        rows = layout.offs[slice] + 1:layout.offs[slice + 1]
        sliced_context.blocks.diagonal[slice] ≈ Fs[rows, rows]
    end
    sliced_audit = HFDMRG._history_physical(Hs, backend, Cs)
    for field in (:energy, :residual, :gram)
        @test isapprox(getproperty(sliced_audit, field), getproperty(sliced_dense, field);
            atol = 3e-11, rtol = 3e-13)
    end
    Dus, Dds = Cus * Cus', Cds * Cds'
    Fus, Fds = copy(Hs), copy(Hsdn)
    HFDMRG.sliced_add_fock_uhf!(Fus, Fds, Dus, Dds, sliced_v, layout)
    FusC, FdsC = Fus * Cus, Fds * Cds
    Rus = FusC - Cus * (Cus' * FusC)
    Rds = FdsC - Cds * (Cds' * FdsC)
    sliced_up_context = HFDMRG._SlicedFockAuditContext(
        Hs, backend, Cus, Cds, :up)
    sliced_dn_context = HFDMRG._SlicedFockAuditContext(
        Hsdn, backend, Cus, Cds, :dn)
    @test HFDMRG._audit_action!(similar(Cus), sliced_up_context, Cus) ≈ FusC
    @test HFDMRG._audit_action!(similar(Cds), sliced_dn_context, Cds) ≈ FdsC
    @test sliced_up_context.Jpack == sliced_dn_context.Jpack
    @test sliced_up_context.setup_reason === sliced_dn_context.setup_reason === :ok
    sliced_u_dense = (; energy = 0.5dot(Dus, Fus + Hs) + 0.5dot(Dds, Fds + Hsdn),
        residual_up = norm(Rus), residual_dn = norm(Rds))
    sliced_u_audit = HFDMRG._history_uhf_physical(Hs, Hsdn, backend, Cus, Cds)
    for field in (:energy, :residual_up, :residual_dn)
        @test isapprox(getproperty(sliced_u_audit, field),
            getproperty(sliced_u_dense, field); atol = 3e-11, rtol = 3e-13)
    end

    Na = 320
    Ha = Matrix(Diagonal(range(-1.0, 1.0; length = Na)))
    Va = [0.01 / (1 + abs(i - j)) for i = 1:Na, j = 1:Na]
    Cua, Cda = orthogonal(Na, 4), orthogonal(Na, 3)
    dense_history_rhf(Ha, Va, Cua)
    HFDMRG._history_physical(Ha, Va, Cua)
    dense_history_uhf(Ha, Ha, Va, Cua, Cda)
    HFDMRG._history_uhf_physical(Ha, Ha, Va, Cua, Cda)
    old_rhf_bytes = @allocated dense_history_rhf(Ha, Va, Cua)
    new_rhf_bytes = @allocated HFDMRG._history_physical(Ha, Va, Cua)
    old_uhf_bytes = @allocated dense_history_uhf(Ha, Ha, Va, Cua, Cda)
    new_uhf_bytes = @allocated HFDMRG._history_uhf_physical(Ha, Ha, Va, Cua, Cda)
    @test 8new_rhf_bytes < old_rhf_bytes
    @test 8new_uhf_bytes < old_uhf_bytes
end
