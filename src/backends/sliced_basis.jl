"""
Sliced-basis backend.

V6[a,b,c,d,n,m] encodes two-electron integrals in a slice-preserving form:
p = (n,a), q = (m,c), r = (n,b), s = (m,d)
(pq|rs) = V6[a,b,c,d,n,m]
"""
struct SlicedBasisBackend
    layout::SliceLayout
    V6::Array{Float64,6}
    nj::Int
    ns::Int
end

struct SlicedBasisWindow
    V6::Array{Float64,6}
    nj::Int
    ns::Int
end

function SlicedBasisBackend(layout::SliceLayout, V6::Array{Float64,6})
    ns = nslices(layout)
    nj = layout.dims[1]
    any(d -> d != nj, layout.dims) && error("layout must have fixed slice size")
    size(V6) == (nj, nj, nj, nj, ns, ns) || error("V6 must have size (nj,nj,nj,nj,ns,ns)")
    SlicedBasisBackend(layout, V6, nj, ns)
end

SlicedBasisWindow(backend::SlicedBasisBackend) =
    SlicedBasisWindow(backend.V6, backend.nj, backend.ns)

function _check_side(side)
    (side === :left || side === :right) || error("side must be :left or :right")
end

function _check_range(ra, name)
    isa(ra, UnitRange) || error("$name must be UnitRange")
    eltype(ra) <: Integer || error("$name must be integer range")
    isempty(ra) && error("$name must be non-empty")
end

function _sliced_not_implemented()
    error("SlicedBasisBackend not implemented yet")
end

function vee_init_block(side, ra, raV, phi, backend::SlicedBasisBackend)
    _check_side(side)
    _check_range(ra, "ra")
    _check_range(raV, "raV")
    _sliced_not_implemented()
end

function vee_absorb_block(side, vee_old, cra, Phi_old, Phi_C, phi_new, raV_new,
    backend::SlicedBasisBackend)
    _check_side(side)
    _check_range(cra, "cra")
    _check_range(raV_new, "raV_new")
    _sliced_not_implemented()
end

function vee_window(Lvee, Rvee, Cra, backend::SlicedBasisBackend)
    _check_range(Cra, "Cra")
    _sliced_not_implemented()
end

function vee_window(::Nothing, ::Nothing, Cra, backend::SlicedBasisBackend)
    _check_range(Cra, "Cra")
    N = backend.nj * backend.ns
    Cra == 1:N || error("SlicedBasisBackend supports full-system center only")
    SlicedBasisWindow(backend)
end

function vee_add_fock!(Fup, Fdn, rhoup, rhodn, win::SlicedBasisWindow)
    error("SlicedBasisBackend UHF not implemented yet")
end

"""
Add RHF mean-field contribution using sliced integrals.

Flattened orbital index mapping:
p = (n,a) with n = 1:ns and a = 1:nj maps to p = (n-1)*nj + a.
V6[a,b,c,d,n,m] corresponds to (p,q|r,s) with p=(n,a), q=(m,c),
r=(n,b), s=(m,d).
"""
function sliced_add_fock_r!(F, rho, V6::Array{Float64,6}, nj::Int, ns::Int)
    N = nj * ns
    size(F, 1) == N || error("F has wrong size")
    size(F, 2) == N || error("F has wrong size")
    size(rho) == size(F) || error("rho has wrong size")
    size(V6) == (nj, nj, nj, nj, ns, ns) || error("V6 has wrong size")

    for np = 1:ns, nq = 1:ns
        for a = 1:nj, c = 1:nj
            p = (np - 1) * nj + a
            q = (nq - 1) * nj + c
            for b = 1:nj, d = 1:nj
                r = (np - 1) * nj + b
                s = (nq - 1) * nj + d
                F[p, q] += 2.0 * rho[r, s] * V6[a, b, c, d, np, nq]
            end
            if np == nq
                # Exchange term uses r,s within a common slice ms.
                for ms = 1:ns, b = 1:nj, d = 1:nj
                    r = (ms - 1) * nj + b
                    s = (ms - 1) * nj + d
                    F[p, q] -= rho[r, s] * V6[a, c, b, d, np, ms]
                end
            end
        end
    end
    F
end

function vee_add_fock_r!(F, rho, win::SlicedBasisWindow)
    sliced_add_fock_r!(F, rho, win.V6, win.nj, win.ns)
end
