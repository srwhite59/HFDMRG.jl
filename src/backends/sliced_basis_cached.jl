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
    SlicedBasisCachedWindow(S, backend.V6, backend.nj, backend.ns)
end

function vee_add_fock_r!(F, rho, win::SlicedBasisCachedWindow)
    S = win.S
    V6 = win.V6
    nj, ns = win.nj, win.ns
    n = size(F, 1)
    size(F) == (n, n) || error("F has wrong size")
    size(rho) == size(F) || error("rho has wrong size")

    for p = 1:n, q = 1:n
        acc = 0.0
        for r = 1:n, s = 1:n
            dir = 0.0
            ex = 0.0
            for n1 = 1:ns, n2 = 1:ns
                S1 = S[n1]
                S2 = S[n2]
                for a = 1:nj, b = 1:nj, c = 1:nj, d = 1:nj
                    v = V6[a, b, c, d, n1, n2]
                    dir += S1[a, p] * S1[b, r] * v * S2[c, q] * S2[d, s]
                    ex += S1[a, p] * S1[b, q] * v * S2[c, r] * S2[d, s]
                end
            end
            acc += rho[r, s] * (2.0 * dir - ex)
        end
        F[p, q] += acc
    end
end

function vee_add_fock!(Fup, Fdn, rhoup, rhodn, win::SlicedBasisCachedWindow)
    S = win.S
    V6 = win.V6
    nj, ns = win.nj, win.ns
    n = size(Fup, 1)
    size(Fup) == (n, n) || error("Fup has wrong size")
    size(Fdn) == (n, n) || error("Fdn has wrong size")
    size(rhoup) == size(Fup) || error("rhoup has wrong size")
    size(rhodn) == size(Fup) || error("rhodn has wrong size")

    rtot = rhoup + rhodn
    for p = 1:n, q = 1:n
        acc_up = 0.0
        acc_dn = 0.0
        for r = 1:n, s = 1:n
            dir = 0.0
            ex = 0.0
            for n1 = 1:ns, n2 = 1:ns
                S1 = S[n1]
                S2 = S[n2]
                for a = 1:nj, b = 1:nj, c = 1:nj, d = 1:nj
                    v = V6[a, b, c, d, n1, n2]
                    dir += S1[a, p] * S1[b, r] * v * S2[c, q] * S2[d, s]
                    ex += S1[a, p] * S1[b, q] * v * S2[c, r] * S2[d, s]
                end
            end
            acc_up += rtot[r, s] * dir - rhoup[r, s] * ex
            acc_dn += rtot[r, s] * dir - rhodn[r, s] * ex
        end
        Fup[p, q] += acc_up
        Fdn[p, q] += acc_dn
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
