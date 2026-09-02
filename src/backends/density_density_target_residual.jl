"""Density-density interaction plus a signed residual in an orthonormal target space."""
struct DensityDensityTargetResidualBackend{B,T<:AbstractFloat}
    base::B
    Q::Matrix{T}
    residual_pair::Matrix{T}
    pair_index::Matrix{Int}
    zero_residual::Bool
end

function DensityDensityTargetResidualBackend(V::AbstractMatrix, Q::AbstractMatrix,
        residual_pair::AbstractMatrix)
    N = size(V, 1)
    size(V, 2) == N || error("V must be square")
    size(Q, 1) == N || error("Q must have one row per physical site")
    m = size(Q, 2)
    m > 0 || error("Q must have at least one target column")
    P = m * (m + 1) ÷ 2
    size(residual_pair) == (P, P) || error("residual_pair must have size ($P, $P)")
    T = eltype(V)
    T <: AbstractFloat || error("V, Q, and residual_pair must share an AbstractFloat type")
    (eltype(Q) === T && eltype(residual_pair) === T) ||
        error("V, Q, and residual_pair must share one element type")
    Qm, Rm = Matrix{T}(Q), Matrix{T}(residual_pair)
    all(isfinite, Qm) || error("Q must be finite")
    all(isfinite, Rm) || error("residual_pair must be finite")
    qtol = 128 * eps(T) * max(N, m)
    norm(transpose(Qm) * Qm - I, Inf) <= qtol || error("Q columns must be orthonormal")
    rtol = 128 * eps(T) * max(one(T), maximum(abs, Rm))
    maximum(abs, Rm - transpose(Rm)) <= rtol || error("residual_pair must be symmetric")
    Rm = (Rm + transpose(Rm)) / 2
    pair_index = zeros(Int, m, m)
    a = 0
    for p = 1:m, q = 1:p
        a += 1
        pair_index[p, q] = pair_index[q, p] = a
    end
    DensityDensityTargetResidualBackend(DensityDensityBackend(V), Qm, Rm,
        pair_index, all(iszero, Rm))
end

struct DDTargetBlockState{B,T}
    base::B
    Q::Matrix{T}
end

struct DDTargetScratch{T}
    Da::Matrix{T}
    Db::Matrix{T}
    D::Matrix{T}
    J::Matrix{T}
    Ka::Matrix{T}
    Kb::Matrix{T}
    G::Matrix{T}
    tmp::Matrix{T}
    dpair::Vector{T}
    jpair::Vector{T}
end

struct DDTargetWindow{W,B,T}
    base::W
    backend::B
    Q::Matrix{T}
    scratch::DDTargetScratch{T}
end

function _target_scratch(T, w, m, P)
    z() = zeros(T, m, m)
    DDTargetScratch(z(), z(), z(), z(), z(), z(), z(), zeros(T, w, m),
        zeros(T, P), zeros(T, P))
end

function vee_init_block(side, ra, raV, phi, backend::DensityDensityTargetResidualBackend)
    base = vee_init_block(side, ra, raV, phi, backend.base)
    DDTargetBlockState(base, Matrix(transpose(phi) * view(backend.Q, ra, :)))
end

function vee_absorb_block(side, old::DDTargetBlockState, cra, Phi_old, Phi_C, phi_new,
        raV_new, backend::DensityDensityTargetResidualBackend)
    base = vee_absorb_block(side, old.base, cra, Phi_old, Phi_C, phi_new, raV_new,
        backend.base)
    Qnew = transpose(Phi_old) * old.Q + transpose(Phi_C) * view(backend.Q, cra, :)
    DDTargetBlockState(base, Matrix(Qnew))
end

function vee_window(L::DDTargetBlockState, R::DDTargetBlockState, Cra,
        backend::DensityDensityTargetResidualBackend)
    base = vee_window(L.base, R.base, Cra, backend.base)
    Q = vcat(L.Q, view(backend.Q, Cra, :), R.Q)
    m, w = size(Q, 2), size(Q, 1)
    DDTargetWindow(base, backend, Q,
        _target_scratch(eltype(Q), w, m, size(backend.residual_pair, 1)))
end

function _target_project!(D, rho, Q, tmp)
    mul!(tmp, rho, Q)
    mul!(D, transpose(Q), tmp)
end

function _target_liftadd!(F, G, Q, tmp)
    mul!(tmp, Q, G)
    mul!(F, tmp, transpose(Q), one(eltype(F)), one(eltype(F)))
end

function _target_contract!(J, Ka, Kb, D, Da, Db, win::DDTargetWindow)
    R, index, s = win.backend.residual_pair, win.backend.pair_index, win.scratch
    m = size(D, 1)
    a = 0
    @inbounds for p = 1:m, q = 1:p
        a += 1
        s.dpair[a] = (p == q ? one(eltype(D)) : 2 * one(eltype(D))) * D[p, q]
    end
    mul!(s.jpair, R, s.dpair)
    a = 0
    @inbounds for p = 1:m, q = 1:p
        a += 1
        J[p, q] = J[q, p] = s.jpair[a]
    end
    if Db === nothing
        @inbounds for r = 1:m, p = 1:r
            va = zero(eltype(D))
            for ss = 1:m, q = 1:m
                va = muladd(R[index[p, q], index[r, ss]], Da[q, ss], va)
            end
            Ka[p, r] = Ka[r, p] = va
        end
    else
        @inbounds for r = 1:m, p = 1:r
            va = zero(eltype(D))
            vb = zero(eltype(D))
            for ss = 1:m, q = 1:m
                v = R[index[p, q], index[r, ss]]
                va = muladd(v, Da[q, ss], va)
                vb = muladd(v, Db[q, ss], vb)
            end
            Ka[p, r] = Ka[r, p] = va
            Kb[p, r] = Kb[r, p] = vb
        end
    end
    nothing
end

function vee_add_fock_r!(F, rho, win::DDTargetWindow)
    vee_add_fock_r!(F, rho, win.base)
    s, Q = win.scratch, win.Q
    _target_project!(s.Da, rho, Q, s.tmp)
    _target_contract!(s.J, s.Ka, s.Kb, s.Da, s.Da, nothing, win)
    @. s.G = 2 * s.J - s.Ka
    _target_liftadd!(F, s.G, Q, s.tmp)
    nothing
end

function vee_add_fock!(Fup, Fdn, rhoup, rhodn, win::DDTargetWindow)
    vee_add_fock!(Fup, Fdn, rhoup, rhodn, win.base)
    s, Q = win.scratch, win.Q
    _target_project!(s.Da, rhoup, Q, s.tmp)
    _target_project!(s.Db, rhodn, Q, s.tmp)
    @. s.D = s.Da + s.Db
    _target_contract!(s.J, s.Ka, s.Kb, s.D, s.Da, s.Db, win)
    @. s.G = s.J - s.Ka
    _target_liftadd!(Fup, s.G, Q, s.tmp)
    @. s.G = s.J - s.Kb
    _target_liftadd!(Fdn, s.G, Q, s.tmp)
    nothing
end

function solve_hfdmrg(H, backend::DensityDensityTargetResidualBackend, psiup0; kwargs...)
    _removed_rhf_route("target-residual-backend")
end

function solve_hfdmrg(H, backend::DensityDensityTargetResidualBackend, psiup0, psidn0;
        kwargs...)
    backend.zero_residual && return solve_hfdmrg(H, backend.base.V, psiup0, psidn0; kwargs...)
    solve_hfdmrg_core(H, backend, psiup0, psidn0; restricted = false, kwargs...)
end

function solve_hfdmrg(Hup, Hdn, backend::DensityDensityTargetResidualBackend, psiup0,
        psidn0; kwargs...)
    Hup === Hdn && return solve_hfdmrg(Hup, backend, psiup0, psidn0; kwargs...)
    backend.zero_residual && return solve_hfdmrg(Hup, Hdn, backend.base.V,
        psiup0, psidn0; kwargs...)
    solve_hfdmrg_core_split(Hup, Hdn, backend, psiup0, psidn0; kwargs...)
end
