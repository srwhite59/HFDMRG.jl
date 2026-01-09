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
    center_slice::Vector{Int}
    center_local::Vector{Int}
    center_cols::Vector{Vector{Int}}
end

function SlicedBasisBackendCached(layout::SliceLayout, V6::Array{Float64,6})
    ns = nslices(layout)
    nj = layout.dims[1]
    any(d -> d != nj, layout.dims) && error("layout must have fixed slice size")
    size(V6) == (nj, nj, nj, nj, ns, ns) || error("V6 must have size (nj,nj,nj,nj,ns,ns)")
    SlicedBasisBackendCached(layout, V6, nj, ns)
end

function _slice_local(layout::SliceLayout, p::Int)
    n = searchsortedlast(layout.offs, p - 1)
    1 <= n <= nslices(layout) || error("orbital index out of bounds")
    a = p - layout.offs[n]
    n, a
end

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
    for s = 1:ns
        Ws = W[s]
        for i = 1:m, j = 1:m, a = 1:nj, b = 1:nj
            acc = 0.0
            for n = 1:ns
                Pn = P[n]
                for c = 1:nj, d = 1:nj
                    acc += Pn[c, i] * Pn[d, j] * backend.V6[c, d, a, b, n, s]
                end
            end
            Ws[i, j, a, b] = acc
        end
    end

    SlicedCachedBlockState(ra, phi, Vijkl, W, P)
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

    SlicedBasisCachedWindow(S, backend.V6, backend.nj, backend.ns, ml, lc, mr,
        Lvee.Vijkl, Rvee.Vijkl, Lvee.W, Rvee.W, center_slice, center_local, center_cols)
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
    for s = 1:ns
        cols = center_cols[s]
        for a = 1:length(cols)
            col_a = cols[a]
            col_a == 0 && continue
            for i = 1:ml
                acc = 0.0
                for b = 1:length(cols)
                    col_b = cols[b]
                    col_b == 0 && continue
                    for j = 1:ml
                        acc += rho[j, col_b] * WL[s][i, j, a, b]
                    end
                end
                F[i, col_a] += 2.0 * acc
            end
        end
    end

    WR = win.WR
    for s = 1:ns
        cols = center_cols[s]
        for a = 1:length(cols)
            col_a = cols[a]
            col_a == 0 && continue
            for i = 1:mr
                acc = 0.0
                for b = 1:length(cols)
                    col_b = cols[b]
                    col_b == 0 && continue
                    for j = 1:mr
                        acc += rho[ml + lc + j, col_b] * WR[s][i, j, a, b]
                    end
                end
                F[ml + lc + i, col_a] += 2.0 * acc
            end
        end
    end

    J = [zeros(Float64, nj, nj) for _ = 1:ns, _ = 1:ns]
    for n1 = 1:ns
        S1 = S[n1]
        for n2 = 1:ns
            S2 = S[n2]
            rho_nm = S1 * rho * S2'
            Jnm = J[n1, n2]
            for a = 1:nj, c = 1:nj
                acc = 0.0
                for b = 1:nj, d = 1:nj
                    acc += V6[a, b, c, d, n1, n2] * rho_nm[b, d]
                end
                Jnm[a, c] = acc
            end
        end
    end

    for n1 = 1:ns
        S1 = S[n1]
        for n2 = 1:ns
            S2 = S[n2]
            Jnm = J[n1, n2]
            F[Crng, Lrng] .+= 2.0 * (S1[:, Crng]' * Jnm * S2[:, Lrng])
            F[Crng, Rrng] .+= 2.0 * (S1[:, Crng]' * Jnm * S2[:, Rrng])
            F[Lrng, Rrng] .+= 2.0 * (S1[:, Lrng]' * Jnm * S2[:, Rrng])
            F[Rrng, Lrng] .+= 2.0 * (S1[:, Rrng]' * Jnm * S2[:, Lrng])
        end
    end

    # Direct uses (p q| r s) and can couple slices in p,q.
    # Exchange uses (p r| q s) and is nonzero only when slice(p) == slice(q).
    rho_slice = Vector{Matrix{Float64}}(undef, ns)
    for n1 = 1:ns
        rho_slice[n1] = S[n1] * rho * S[n1]'
    end
    for n1 = 1:ns
        K = zeros(Float64, nj, nj)
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
        F .-= S[n1]' * K * S[n1]
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
    for s = 1:ns
        cols = center_cols[s]
        for a = 1:length(cols)
            col_a = cols[a]
            col_a == 0 && continue
            for i = 1:ml
                acc = 0.0
                for b = 1:length(cols)
                    col_b = cols[b]
                    col_b == 0 && continue
                    for j = 1:ml
                        acc += rtot[j, col_b] * WL[s][i, j, a, b]
                    end
                end
                Fup[i, col_a] += acc
                Fdn[i, col_a] += acc
            end
        end
    end

    WR = win.WR
    for s = 1:ns
        cols = center_cols[s]
        for a = 1:length(cols)
            col_a = cols[a]
            col_a == 0 && continue
            for i = 1:mr
                acc = 0.0
                for b = 1:length(cols)
                    col_b = cols[b]
                    col_b == 0 && continue
                    for j = 1:mr
                        acc += rtot[ml + lc + j, col_b] * WR[s][i, j, a, b]
                    end
                end
                Fup[ml + lc + i, col_a] += acc
                Fdn[ml + lc + i, col_a] += acc
            end
        end
    end

    J = [zeros(Float64, nj, nj) for _ = 1:ns, _ = 1:ns]
    for n1 = 1:ns
        S1 = S[n1]
        for n2 = 1:ns
            S2 = S[n2]
            rho_nm = S1 * rtot * S2'
            Jnm = J[n1, n2]
            for a = 1:nj, c = 1:nj
                acc = 0.0
                for b = 1:nj, d = 1:nj
                    acc += V6[a, b, c, d, n1, n2] * rho_nm[b, d]
                end
                Jnm[a, c] = acc
            end
        end
    end

    for n1 = 1:ns
        S1 = S[n1]
        for n2 = 1:ns
            S2 = S[n2]
            Jnm = J[n1, n2]
            block_CL = S1[:, Crng]' * Jnm * S2[:, Lrng]
            block_CR = S1[:, Crng]' * Jnm * S2[:, Rrng]
            block_LR = S1[:, Lrng]' * Jnm * S2[:, Rrng]
            block_RL = S1[:, Rrng]' * Jnm * S2[:, Lrng]
            Fup[Crng, Lrng] .+= block_CL
            Fdn[Crng, Lrng] .+= block_CL
            Fup[Crng, Rrng] .+= block_CR
            Fdn[Crng, Rrng] .+= block_CR
            Fup[Lrng, Rrng] .+= block_LR
            Fdn[Lrng, Rrng] .+= block_LR
            Fup[Rrng, Lrng] .+= block_RL
            Fdn[Rrng, Lrng] .+= block_RL
        end
    end

    rho_slice_up = Vector{Matrix{Float64}}(undef, ns)
    rho_slice_dn = Vector{Matrix{Float64}}(undef, ns)
    for n1 = 1:ns
        Sn = S[n1]
        rho_slice_up[n1] = Sn * rhoup * Sn'
        rho_slice_dn[n1] = Sn * rhodn * Sn'
    end
    # Direct uses (p q| r s) and can couple slices in p,q.
    # Exchange uses (p r| q s) and is nonzero only when slice(p) == slice(q).
    for n1 = 1:ns
        Kup = zeros(Float64, nj, nj)
        Kdn = zeros(Float64, nj, nj)
        for n2 = 1:ns
            rup = rho_slice_up[n2]
            rdn = rho_slice_dn[n2]
            for a = 1:nj, b = 1:nj
                acc_up = 0.0
                acc_dn = 0.0
                for c = 1:nj, d = 1:nj
                    v = V6[a, b, c, d, n1, n2]
                    acc_up += rup[c, d] * v
                    acc_dn += rdn[c, d] * v
                end
                Kup[a, b] += acc_up
                Kdn[a, b] += acc_dn
            end
        end
        Fup .-= S[n1]' * Kup * S[n1]
        Fdn .-= S[n1]' * Kdn * S[n1]
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
