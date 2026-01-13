"""
Cached sliced-basis backend (correctness-first).

Uses per-slice basis coefficients to build superblock Fock contributions without
embedding into the full orbital basis.
"""
struct SlicedBasisBackendCached
    layout::SliceLayout
    V6::Array{Float64,6}
    nj::Int
    ns::Int
end

struct SlicedCachedBlockState
    ra::UnitRange{Int}
    phi::Matrix{Float64}
    Vijkl::Array{Float64,4}
    W::Vector{Array{Float64,4}}
    # Wswap[s][i,j,a,b] uses the swapped orientation V6[a,b,c,d,s,n] to avoid
    # assuming additional symmetries when building block-center direct terms.
    Wswap::Vector{Array{Float64,4}}
    P::Vector{Matrix{Float64}}
end

struct SlicedBasisCachedWindow
    S::Vector{Matrix{Float64}}
    V6::Array{Float64,6}
    nj::Int
    ns::Int
    ml::Int
    lc::Int
    mr::Int
    VLL::Array{Float64,4}
    VRR::Array{Float64,4}
    WL::Vector{Array{Float64,4}}
    WR::Vector{Array{Float64,4}}
    WLswap::Vector{Array{Float64,4}}
    WRswap::Vector{Array{Float64,4}}
    VLR::Array{Float64,4}
    VRL::Array{Float64,4}
    center_slice::Vector{Int}
    center_local::Vector{Int}
    center_cols::Vector{Vector{Int}}
    Srho::Vector{Matrix{Float64}}
    Srho_up::Vector{Matrix{Float64}}
    Srho_dn::Vector{Matrix{Float64}}
    rho_slice::Vector{Matrix{Float64}}
    rho_slice_up::Vector{Matrix{Float64}}
    rho_slice_dn::Vector{Matrix{Float64}}
    K::Matrix{Float64}
    tmp_nj_n::Matrix{Float64}
end

function SlicedBasisBackendCached(layout::SliceLayout, V6::Array{Float64,6})
    ns = nslices(layout)
    nj = layout.dims[1]
    any(d -> d != nj, layout.dims) && error("layout must have fixed slice size")
    size(V6) == (nj, nj, nj, nj, ns, ns) || error("V6 must have size (nj,nj,nj,nj,ns,ns)")
    SlicedBasisBackendCached(layout, V6, nj, ns)
end

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

function _pair_map_mat!(U::Matrix{Float64}, P::Matrix{Float64})
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

function _pair_map_mat(P::Matrix{Float64})
    nj, m = size(P)
    U = zeros(Float64, nj * nj, m * m)
    _pair_map_mat!(U, P)
end

function _build_cached_state(ra::UnitRange{Int}, phi::Matrix{Float64},
    backend::SlicedBasisBackendCached)
    ns = backend.ns
    nj = backend.nj
    m = size(phi, 2)
    P = _build_slice_basis(backend.layout, ra, phi)

    Vijkl = zeros(Float64, m, m, m, m)
    for n = 1:ns, m2 = 1:ns
        Pn = P[n]
        Pm = P[m2]
        for i = 1:m, j = 1:m, k = 1:m, l = 1:m
            acc = 0.0
            for a = 1:nj, b = 1:nj, c = 1:nj, d = 1:nj
                acc += Pn[a, i] * Pn[b, k] * backend.V6[a, b, c, d, n, m2] *
                       Pm[c, j] * Pm[d, l]
            end
            Vijkl[i, j, k, l] += acc
        end
    end

    W = [zeros(Float64, m, m, nj, nj) for _ = 1:ns]
    Wswap = [zeros(Float64, m, m, nj, nj) for _ = 1:ns]
    for s = 1:ns
        Ws = W[s]
        Wt = Wswap[s]
        for i = 1:m, j = 1:m, a = 1:nj, b = 1:nj
            acc = 0.0
            acc_swap = 0.0
            for n = 1:ns
                Pn = P[n]
                for c = 1:nj, d = 1:nj
                    acc += Pn[c, i] * Pn[d, j] * backend.V6[c, d, a, b, n, s]
                    acc_swap += Pn[c, i] * Pn[d, j] * backend.V6[a, b, c, d, s, n]
                end
            end
            Ws[i, j, a, b] = acc
            Wt[i, j, a, b] = acc_swap
        end
    end

    SlicedCachedBlockState(ra, phi, Vijkl, W, Wswap, P)
end

function vee_init_block(side, ra, raV, phi, backend::SlicedBasisBackendCached)
    _check_side(side)
    _check_range(ra, "ra")
    _check_range(raV, "raV")
    size(phi, 1) == length(ra) || error("phi has wrong row count for ra")
    _build_cached_state(ra, Matrix{Float64}(phi), backend)
end

function vee_absorb_block(side, vee_old::SlicedCachedBlockState, cra, Phi_old, Phi_C,
    phi_new, raV_new, backend::SlicedBasisBackendCached)
    _check_side(side)
    _check_range(cra, "cra")
    _check_range(raV_new, "raV_new")
    oldra = vee_old.ra
    newra = side == :left ? (oldra[1]:cra[end]) : (cra[1]:oldra[end])
    size(phi_new, 1) == length(newra) || error("phi_new has wrong row count for new ra")
    _build_cached_state(newra, Matrix{Float64}(phi_new), backend)
end

function vee_window(Lvee::SlicedCachedBlockState, Rvee::SlicedCachedBlockState,
    Cra, backend::SlicedBasisBackendCached)
    _check_range(Cra, "Cra")
    Lra = Lvee.ra
    Rra = Rvee.ra
    Lra[end] < Cra[1] || error("Lra must end before Cra")
    Cra[end] < Rra[1] || error("Cra must end before Rra")
    N = backend.nj * backend.ns
    Rra[end] <= N || error("ranges exceed full basis size")

    ml = size(Lvee.phi, 2)
    mr = size(Rvee.phi, 2)
    lc = length(Cra)
    superdim = ml + lc + mr
    S = [zeros(Float64, backend.nj, superdim) for _ = 1:backend.ns]
    for s = 1:backend.ns
        S[s][:, 1:ml] = Lvee.P[s]
        S[s][:, ml + lc + 1:end] = Rvee.P[s]
    end
    for (idx, p) in enumerate(Cra)
        n, a = _slice_local(backend.layout, p)
        S[n][a, ml + idx] = 1.0
    end
    center_slice = zeros(Int, lc)
    center_local = zeros(Int, lc)
    center_cols = [fill(0, backend.layout.dims[s]) for s in 1:backend.ns]
    for (idx, p) in enumerate(Cra)
        n, a = _slice_local(backend.layout, p)
        center_slice[idx] = n
        center_local[idx] = a
        center_cols[n][a] = ml + idx
    end

    Srho = [zeros(Float64, backend.nj, superdim) for _ = 1:backend.ns]
    Srho_up = [zeros(Float64, backend.nj, superdim) for _ = 1:backend.ns]
    Srho_dn = [zeros(Float64, backend.nj, superdim) for _ = 1:backend.ns]
    rho_slice = [zeros(Float64, backend.nj, backend.nj) for _ = 1:backend.ns]
    rho_slice_up = [zeros(Float64, backend.nj, backend.nj) for _ = 1:backend.ns]
    rho_slice_dn = [zeros(Float64, backend.nj, backend.nj) for _ = 1:backend.ns]
    K = zeros(Float64, backend.nj, backend.nj)
    tmp_nj_n = zeros(Float64, backend.nj, superdim)

    VLR = zeros(Float64, ml, ml, mr, mr)
    VRL = zeros(Float64, mr, mr, ml, ml)
    VLR_mat = reshape(VLR, ml * ml, mr * mr)
    VRL_mat = reshape(VRL, mr * mr, ml * ml)
    Utmp_R = zeros(Float64, backend.nj * backend.nj, mr * mr)
    Utmp_L = zeros(Float64, backend.nj * backend.nj, ml * ml)
    for s = 1:backend.ns
        if any(!iszero, Rvee.P[s])
            _pair_map_mat!(Utmp_R, Rvee.P[s])
            WL_mat = reshape(Lvee.W[s], ml * ml, backend.nj * backend.nj)
            mul!(VLR_mat, WL_mat, Utmp_R, 1.0, 1.0)
        end
        if any(!iszero, Lvee.P[s])
            _pair_map_mat!(Utmp_L, Lvee.P[s])
            WR_mat = reshape(Rvee.W[s], mr * mr, backend.nj * backend.nj)
            mul!(VRL_mat, WR_mat, Utmp_L, 1.0, 1.0)
        end
    end

    SlicedBasisCachedWindow(S, backend.V6, backend.nj, backend.ns, ml, lc, mr,
        Lvee.Vijkl, Rvee.Vijkl, Lvee.W, Rvee.W, Lvee.Wswap, Rvee.Wswap, VLR, VRL,
        center_slice, center_local, center_cols, Srho, Srho_up, Srho_dn, rho_slice,
        rho_slice_up, rho_slice_dn, K, tmp_nj_n)
end

function vee_add_fock_r!(F, rho, win::SlicedBasisCachedWindow)
    S = win.S
    V6 = win.V6
    nj, ns = win.nj, win.ns
    ml, lc, mr = win.ml, win.lc, win.mr
    n = size(F, 1)
    size(F) == (n, n) || error("F has wrong size")
    size(rho) == size(F) || error("rho has wrong size")
    ml + lc + mr == n || error("window dimensions do not match F")

    Lrng = 1:ml
    Crng = ml + 1:ml + lc
    Rrng = ml + lc + 1:n

    # Direct (J) contributions.
    VLL = win.VLL
    for i = 1:ml, j = 1:ml
        acc = 0.0
        for k = 1:ml, l = 1:ml
            acc += rho[k, l] * VLL[i, j, k, l]
        end
        F[i, j] += 2.0 * acc
    end

    VRR = win.VRR
    for i = 1:mr, j = 1:mr
        acc = 0.0
        for k = 1:mr, l = 1:mr
            acc += rho[ml + lc + k, ml + lc + l] * VRR[i, j, k, l]
        end
        F[ml + lc + i, ml + lc + j] += 2.0 * acc
    end

    center_slice = win.center_slice
    center_local = win.center_local
    center_cols = win.center_cols
    for pidx = 1:lc
        np = center_slice[pidx]
        a = center_local[pidx]
        p = ml + pidx
        for qidx = 1:lc
            nq = center_slice[qidx]
            c = center_local[qidx]
            q = ml + qidx
            acc = 0.0
            cols_p = center_cols[np]
            cols_q = center_cols[nq]
            for b = 1:length(cols_p), d = 1:length(cols_q)
                rcol = cols_p[b]
                scol = cols_q[d]
                if rcol != 0 && scol != 0
                    acc += rho[rcol, scol] * V6[a, b, c, d, np, nq]
                end
            end
            F[p, q] += 2.0 * acc
        end
    end

    WL = win.WL
    WLswap = win.WLswap
    for s = 1:ns
        cols = center_cols[s]
        for a = 1:length(cols)
            col_a = cols[a]
            col_a == 0 && continue
            for i = 1:ml
                acc_lc = 0.0
                acc_cl = 0.0
                for b = 1:length(cols)
                    col_b = cols[b]
                    col_b == 0 && continue
                    for j = 1:ml
                        acc_lc += rho[j, col_b] * WL[s][i, j, a, b]
                        acc_cl += rho[col_b, j] * WLswap[s][i, j, a, b]
                    end
                end
                F[i, col_a] += 2.0 * acc_lc
                F[col_a, i] += 2.0 * acc_cl
            end
        end
    end

    WR = win.WR
    WRswap = win.WRswap
    for s = 1:ns
        cols = center_cols[s]
        for a = 1:length(cols)
            col_a = cols[a]
            col_a == 0 && continue
            for i = 1:mr
                acc_rc = 0.0
                acc_cr = 0.0
                for b = 1:length(cols)
                    col_b = cols[b]
                    col_b == 0 && continue
                    for j = 1:mr
                        acc_rc += rho[ml + lc + j, col_b] * WR[s][i, j, a, b]
                        acc_cr += rho[col_b, ml + lc + j] * WRswap[s][i, j, a, b]
                    end
                end
                F[ml + lc + i, col_a] += 2.0 * acc_rc
                F[col_a, ml + lc + i] += 2.0 * acc_cr
            end
        end
    end

    Srho = win.Srho
    tmp_nj_n = win.tmp_nj_n
    VLR = win.VLR
    VRL = win.VRL

    @views for n1 = 1:ns
        mul!(Srho[n1], S[n1], rho, 1.0, 0.0)
    end

    for i = 1:ml, k = 1:mr
        acc = 0.0
        for j = 1:ml, l = 1:mr
            acc += rho[j, ml + lc + l] * VLR[i, j, k, l]
        end
        F[i, ml + lc + k] += 2.0 * acc
    end
    for k = 1:mr, i = 1:ml
        acc = 0.0
        for l = 1:mr, j = 1:ml
            acc += rho[ml + lc + l, j] * VRL[k, l, i, j]
        end
        F[ml + lc + k, i] += 2.0 * acc
    end

    # Direct uses (p q| r s) and can couple slices in p,q.
    # Exchange uses (p r| q s) and is nonzero only when slice(p) == slice(q).
    rho_slice = win.rho_slice
    K = win.K
    @views for n1 = 1:ns
        mul!(rho_slice[n1], Srho[n1], S[n1]', 1.0, 0.0)
    end
    @views for n1 = 1:ns
        fill!(K, 0.0)
        for n2 = 1:ns
            rho_m = rho_slice[n2]
            for a = 1:nj, b = 1:nj
                acc = 0.0
                for c = 1:nj, d = 1:nj
                    acc += rho_m[c, d] * V6[a, b, c, d, n1, n2]
                end
                K[a, b] += acc
            end
        end
        mul!(tmp_nj_n, K, S[n1], 1.0, 0.0)
        mul!(F, S[n1]', tmp_nj_n, -1.0, 1.0)
    end
end

function vee_add_fock!(Fup, Fdn, rhoup, rhodn, win::SlicedBasisCachedWindow)
    S = win.S
    V6 = win.V6
    nj, ns = win.nj, win.ns
    ml, lc, mr = win.ml, win.lc, win.mr
    n = size(Fup, 1)
    size(Fup) == (n, n) || error("Fup has wrong size")
    size(Fdn) == (n, n) || error("Fdn has wrong size")
    size(rhoup) == size(Fup) || error("rhoup has wrong size")
    size(rhodn) == size(Fup) || error("rhodn has wrong size")
    ml + lc + mr == n || error("window dimensions do not match F")

    Lrng = 1:ml
    Crng = ml + 1:ml + lc
    Rrng = ml + lc + 1:n

    rtot = rhoup + rhodn

    Srho_up = win.Srho_up
    Srho_dn = win.Srho_dn
    tmp_nj_n = win.tmp_nj_n
    VLR = win.VLR
    VRL = win.VRL

    @views for n1 = 1:ns
        S1 = S[n1]
        mul!(Srho_up[n1], S1, rhoup, 1.0, 0.0)
        mul!(Srho_dn[n1], S1, rhodn, 1.0, 0.0)
    end

    VLL = win.VLL
    for i = 1:ml, j = 1:ml
        acc = 0.0
        for k = 1:ml, l = 1:ml
            acc += rtot[k, l] * VLL[i, j, k, l]
        end
        Fup[i, j] += acc
        Fdn[i, j] += acc
    end

    VRR = win.VRR
    for i = 1:mr, j = 1:mr
        acc = 0.0
        for k = 1:mr, l = 1:mr
            acc += rtot[ml + lc + k, ml + lc + l] * VRR[i, j, k, l]
        end
        Fup[ml + lc + i, ml + lc + j] += acc
        Fdn[ml + lc + i, ml + lc + j] += acc
    end

    center_slice = win.center_slice
    center_local = win.center_local
    center_cols = win.center_cols
    for pidx = 1:lc
        np = center_slice[pidx]
        a = center_local[pidx]
        p = ml + pidx
        for qidx = 1:lc
            nq = center_slice[qidx]
            c = center_local[qidx]
            q = ml + qidx
            acc = 0.0
            cols_p = center_cols[np]
            cols_q = center_cols[nq]
            for b = 1:length(cols_p), d = 1:length(cols_q)
                rcol = cols_p[b]
                scol = cols_q[d]
                if rcol != 0 && scol != 0
                    acc += rtot[rcol, scol] * V6[a, b, c, d, np, nq]
                end
            end
            Fup[p, q] += acc
            Fdn[p, q] += acc
        end
    end

    WL = win.WL
    WLswap = win.WLswap
    for s = 1:ns
        cols = center_cols[s]
        for a = 1:length(cols)
            col_a = cols[a]
            col_a == 0 && continue
            for i = 1:ml
                acc_lc = 0.0
                acc_cl = 0.0
                for b = 1:length(cols)
                    col_b = cols[b]
                    col_b == 0 && continue
                    for j = 1:ml
                        acc_lc += rtot[j, col_b] * WL[s][i, j, a, b]
                        acc_cl += rtot[col_b, j] * WLswap[s][i, j, a, b]
                    end
                end
                Fup[i, col_a] += acc_lc
                Fdn[i, col_a] += acc_lc
                Fup[col_a, i] += acc_cl
                Fdn[col_a, i] += acc_cl
            end
        end
    end

    WR = win.WR
    WRswap = win.WRswap
    for s = 1:ns
        cols = center_cols[s]
        for a = 1:length(cols)
            col_a = cols[a]
            col_a == 0 && continue
            for i = 1:mr
                acc_rc = 0.0
                acc_cr = 0.0
                for b = 1:length(cols)
                    col_b = cols[b]
                    col_b == 0 && continue
                    for j = 1:mr
                        acc_rc += rtot[ml + lc + j, col_b] * WR[s][i, j, a, b]
                        acc_cr += rtot[col_b, ml + lc + j] * WRswap[s][i, j, a, b]
                    end
                end
                Fup[ml + lc + i, col_a] += acc_rc
                Fdn[ml + lc + i, col_a] += acc_rc
                Fup[col_a, ml + lc + i] += acc_cr
                Fdn[col_a, ml + lc + i] += acc_cr
            end
        end
    end

    for i = 1:ml, k = 1:mr
        acc = 0.0
        for j = 1:ml, l = 1:mr
            acc += rtot[j, ml + lc + l] * VLR[i, j, k, l]
        end
        Fup[i, ml + lc + k] += acc
        Fdn[i, ml + lc + k] += acc
    end
    for k = 1:mr, i = 1:ml
        acc = 0.0
        for l = 1:mr, j = 1:ml
            acc += rtot[ml + lc + l, j] * VRL[k, l, i, j]
        end
        Fup[ml + lc + k, i] += acc
        Fdn[ml + lc + k, i] += acc
    end

    rho_slice_up = win.rho_slice_up
    rho_slice_dn = win.rho_slice_dn
    @views for n1 = 1:ns
        Sn = S[n1]
        mul!(rho_slice_up[n1], Srho_up[n1], Sn', 1.0, 0.0)
        mul!(rho_slice_dn[n1], Srho_dn[n1], Sn', 1.0, 0.0)
    end
    # Direct uses (p q| r s) and can couple slices in p,q.
    # Exchange uses (p r| q s) and is nonzero only when slice(p) == slice(q).
    K = win.K
    @views for n1 = 1:ns
        fill!(K, 0.0)
        for n2 = 1:ns
            rup = rho_slice_up[n2]
            for a = 1:nj, b = 1:nj
                acc = 0.0
                for c = 1:nj, d = 1:nj
                    acc += rup[c, d] * V6[a, b, c, d, n1, n2]
                end
                K[a, b] += acc
            end
        end
        mul!(tmp_nj_n, K, S[n1], 1.0, 0.0)
        mul!(Fup, S[n1]', tmp_nj_n, -1.0, 1.0)

        fill!(K, 0.0)
        for n2 = 1:ns
            rdn = rho_slice_dn[n2]
            for a = 1:nj, b = 1:nj
                acc = 0.0
                for c = 1:nj, d = 1:nj
                    acc += rdn[c, d] * V6[a, b, c, d, n1, n2]
                end
                K[a, b] += acc
            end
        end
        mul!(tmp_nj_n, K, S[n1], 1.0, 0.0)
        mul!(Fdn, S[n1]', tmp_nj_n, -1.0, 1.0)
    end
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
