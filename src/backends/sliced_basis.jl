"""
Sliced-basis backend.

Fixed-size form:
V6[a,b,c,d,n,m] encodes two-electron integrals in a slice-preserving form:
p = (n,a), q = (m,c), r = (n,b), s = (m,d)
(pq|rs) = V6[a,b,c,d,n,m]

Ragged form:
Vblocks[n][m] has size dims[n]×dims[n]×dims[m]×dims[m] with the same mapping.
"""
struct SlicedVeeRagged
    layout::SliceLayout
    V::Vector{Vector{Array{Float64,4}}}
    function SlicedVeeRagged(layout::SliceLayout, V::Vector{Vector{Array{Float64,4}}})
        ns = nslices(layout)
        length(V) == ns || error("Vblocks must have length ns")
        for n = 1:ns
            length(V[n]) == ns || error("Vblocks[$n] must have length ns")
        end
        dims = layout.dims
        for n = 1:ns, m = 1:ns
            size(V[n][m]) == (dims[n], dims[n], dims[m], dims[m]) ||
                error("Vblocks[$n][$m] has wrong size")
        end
        new(layout, V)
    end
end

struct SlicedBasisBackend{T}
    layout::SliceLayout
    V::T
end

struct SlicedBasisWindow{T}
    B::Matrix{Float64}
    V::T
    layout::SliceLayout
end

function SlicedBasisBackend(layout::SliceLayout, V6::Array{Float64,6})
    ns = nslices(layout)
    nj = layout.dims[1]
    any(d -> d != nj, layout.dims) && error("layout must have fixed slice size")
    size(V6) == (nj, nj, nj, nj, ns, ns) || error("V6 must have size (nj,nj,nj,nj,ns,ns)")
    SlicedBasisBackend{Array{Float64,6}}(layout, V6)
end

function SlicedBasisBackend(layout::SliceLayout, vee::SlicedVeeRagged)
    layout.dims == vee.layout.dims || error("layout does not match ragged interaction")
    SlicedBasisBackend{SlicedVeeRagged}(layout, vee)
end

SlicedBasisBackend(vee::SlicedVeeRagged) = SlicedBasisBackend(vee.layout, vee)

SlicedBasisBackend(layout::SliceLayout, Vblocks::Vector{Vector{Array{Float64,4}}}) =
    SlicedBasisBackend(SlicedVeeRagged(layout, Vblocks))

struct SlicedBlockState
    ra::UnitRange{Int}
    phi::Matrix{Float64}
end

SlicedBasisWindow(B::AbstractMatrix, backend::SlicedBasisBackend) =
    SlicedBasisWindow(Matrix{Float64}(B), backend.V, backend.layout)

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
    size(phi, 1) == length(ra) || error("phi has wrong row count for ra")
    SlicedBlockState(ra, Matrix{Float64}(phi))
end

function vee_absorb_block(side, vee_old, cra, Phi_old, Phi_C, phi_new, raV_new,
    backend::SlicedBasisBackend)
    _check_side(side)
    _check_range(cra, "cra")
    _check_range(raV_new, "raV_new")
    oldra = vee_old.ra
    newra = side == :left ? (oldra[1]:cra[end]) : (cra[1]:oldra[end])
    size(phi_new, 1) == length(newra) || error("phi_new has wrong row count for new ra")
    SlicedBlockState(newra, Matrix{Float64}(phi_new))
end

function vee_window(Lvee, Rvee, Cra, backend::SlicedBasisBackend)
    _check_range(Cra, "Cra")
    Lra = Lvee.ra
    Rra = Rvee.ra
    Lra[end] < Cra[1] || error("Lra must end before Cra")
    Cra[end] < Rra[1] || error("Cra must end before Rra")
    N = backend.layout.offs[end]
    Rra[end] <= N || error("ranges exceed full basis size")

    Lphi = Lvee.phi
    Rphi = Rvee.phi
    size(Lphi, 1) == length(Lra) || error("Lphi has wrong row count")
    size(Rphi, 1) == length(Rra) || error("Rphi has wrong row count")
    ml = size(Lphi, 2)
    mr = size(Rphi, 2)
    lc = length(Cra)
    B = zeros(Float64, N, ml + lc + mr)
    B[Lra, 1:ml] = Lphi
    for (i, idx) in enumerate(Cra)
        B[idx, ml + i] = 1.0
    end
    B[Rra, ml + lc + 1:end] = Rphi
    SlicedBasisWindow(B, backend)
end

function vee_window(::Nothing, ::Nothing, Cra, backend::SlicedBasisBackend)
    _check_range(Cra, "Cra")
    N = backend.layout.offs[end]
    Cra == 1:N || error("SlicedBasisBackend supports full-system center only")
    B = Matrix{Float64}(I, N, N)
    SlicedBasisWindow(B, backend)
end

function vee_add_fock!(Fup, Fdn, rhoup, rhodn, win::SlicedBasisWindow)
    B = win.B
    rhoup_full = B * rhoup * B'
    rhodn_full = B * rhodn * B'
    N = size(B, 1)
    Gup_full = zeros(Float64, N, N)
    Gdn_full = zeros(Float64, N, N)
    sliced_add_fock_uhf!(Gup_full, Gdn_full, rhoup_full, rhodn_full, win.V, win.layout)
    Fup .+= B' * Gup_full * B
    Fdn .+= B' * Gdn_full * B
end

"""
Add RHF mean-field contribution using sliced integrals.

Flattened orbital index mapping:
p = (n,a) with n = 1:ns and a = 1:nj maps to p = (n-1)*nj + a.
V6[a,b,c,d,n,m] corresponds to V[p,q,r,s] with p=(n,a), q=(m,c),
r=(n,b), s=(m,d), in the Hamiltonian convention
0.5 * sum(V[p,q,r,s] * cdag[p] * cdag[q] * c[s] * c[r]).
"""
function sliced_add_fock_r!(F, rho, V6::Array{Float64,6}, nj::Int, ns::Int)
    N = nj * ns
    size(F, 1) == N || error("F has wrong size")
    size(F, 2) == N || error("F has wrong size")
    size(rho) == size(F) || error("rho has wrong size")
    size(V6) == (nj, nj, nj, nj, ns, ns) || error("V6 has wrong size")

    for n = 1:ns, m = 1:ns, a = 1:nj, b = 1:nj, c = 1:nj, d = 1:nj
        p = (n - 1) * nj + a
        r = (n - 1) * nj + b
        q = (m - 1) * nj + c
        s = (m - 1) * nj + d
        v = V6[a, b, c, d, n, m]
        F[p, r] += 2.0 * rho[q, s] * v
        F[p, s] -= rho[r, q] * v
    end
    F
end

function sliced_add_fock_r!(F, rho, V6::Array{Float64,6}, layout::SliceLayout)
    nj = layout.dims[1]
    any(d -> d != nj, layout.dims) && error("layout must have fixed slice size")
    ns = nslices(layout)
    sliced_add_fock_r!(F, rho, V6, nj, ns)
end

function sliced_add_fock_r!(F, rho, vee::SlicedVeeRagged)
    layout = vee.layout
    dims = layout.dims
    offs = layout.offs
    ns = nslices(layout)
    N = offs[end]
    size(F, 1) == N || error("F has wrong size")
    size(F, 2) == N || error("F has wrong size")
    size(rho) == size(F) || error("rho has wrong size")

    Vblocks = vee.V
    for n = 1:ns, m = 1:ns
        dimn = dims[n]
        dimm = dims[m]
        offn = offs[n]
        offm = offs[m]
        Vnm = Vblocks[n][m]
        for a = 1:dimn, b = 1:dimn, c = 1:dimm, d = 1:dimm
            p = offn + a
            r = offn + b
            q = offm + c
            s = offm + d
            v = Vnm[a, b, c, d]
            F[p, r] += 2.0 * rho[q, s] * v
            F[p, s] -= rho[r, q] * v
        end
    end
    F
end

function sliced_add_fock_r!(F, rho, vee::SlicedVeeRagged, layout::SliceLayout)
    layout.dims == vee.layout.dims || error("layout does not match ragged interaction")
    sliced_add_fock_r!(F, rho, vee)
end

"""
Add UHF mean-field contributions using sliced integrals.

For V[p,q,r,s] in the documented Hamiltonian convention:
Fspin[p,r] += sum_{q,s} (rhoup + rhodn)[q,s] * V[p,q,r,s]
Fspin[p,s] -= sum_{q,r} rhospin[r,q] * V[p,q,r,s].
"""
function sliced_add_fock_uhf!(Fup, Fdn, rhoup, rhodn, V6::Array{Float64,6}, nj::Int, ns::Int)
    N = nj * ns
    size(Fup, 1) == N || error("Fup has wrong size")
    size(Fup, 2) == N || error("Fup has wrong size")
    size(Fdn) == size(Fup) || error("Fdn has wrong size")
    size(rhoup) == size(Fup) || error("rhoup has wrong size")
    size(rhodn) == size(Fup) || error("rhodn has wrong size")
    size(V6) == (nj, nj, nj, nj, ns, ns) || error("V6 has wrong size")

    rtot = rhoup + rhodn
    for n = 1:ns, m = 1:ns, a = 1:nj, b = 1:nj, c = 1:nj, d = 1:nj
        p = (n - 1) * nj + a
        r = (n - 1) * nj + b
        q = (m - 1) * nj + c
        s = (m - 1) * nj + d
        v = V6[a, b, c, d, n, m]
        Fup[p, r] += rtot[q, s] * v
        Fdn[p, r] += rtot[q, s] * v
        Fup[p, s] -= rhoup[r, q] * v
        Fdn[p, s] -= rhodn[r, q] * v
    end
    Fup, Fdn
end

function sliced_add_fock_uhf!(Fup, Fdn, rhoup, rhodn, V6::Array{Float64,6}, layout::SliceLayout)
    nj = layout.dims[1]
    any(d -> d != nj, layout.dims) && error("layout must have fixed slice size")
    ns = nslices(layout)
    sliced_add_fock_uhf!(Fup, Fdn, rhoup, rhodn, V6, nj, ns)
end

function sliced_add_fock_uhf!(Fup, Fdn, rhoup, rhodn, vee::SlicedVeeRagged)
    layout = vee.layout
    dims = layout.dims
    offs = layout.offs
    ns = nslices(layout)
    N = offs[end]
    size(Fup, 1) == N || error("Fup has wrong size")
    size(Fup, 2) == N || error("Fup has wrong size")
    size(Fdn) == size(Fup) || error("Fdn has wrong size")
    size(rhoup) == size(Fup) || error("rhoup has wrong size")
    size(rhodn) == size(Fup) || error("rhodn has wrong size")

    rtot = rhoup + rhodn
    Vblocks = vee.V
    for n = 1:ns, m = 1:ns
        dimn = dims[n]
        dimm = dims[m]
        offn = offs[n]
        offm = offs[m]
        Vnm = Vblocks[n][m]
        for a = 1:dimn, b = 1:dimn, c = 1:dimm, d = 1:dimm
            p = offn + a
            r = offn + b
            q = offm + c
            s = offm + d
            v = Vnm[a, b, c, d]
            Fup[p, r] += rtot[q, s] * v
            Fdn[p, r] += rtot[q, s] * v
            Fup[p, s] -= rhoup[r, q] * v
            Fdn[p, s] -= rhodn[r, q] * v
        end
    end
    Fup, Fdn
end

function sliced_add_fock_uhf!(Fup, Fdn, rhoup, rhodn, vee::SlicedVeeRagged,
    layout::SliceLayout)
    layout.dims == vee.layout.dims || error("layout does not match ragged interaction")
    sliced_add_fock_uhf!(Fup, Fdn, rhoup, rhodn, vee)
end

function vee_add_fock_r!(F, rho, win::SlicedBasisWindow)
    B = win.B
    rho_full = B * rho * B'
    N = size(B, 1)
    G_full = zeros(Float64, N, N)
    sliced_add_fock_r!(G_full, rho_full, win.V, win.layout)
    F .+= B' * G_full * B
end

"""
solve_hfdmrg(H, backend::SlicedBasisBackend, psiup0; kwargs...) -> (psiup, psidn, energy)

Restricted HF (RHF) sweep using a sliced-basis backend.
"""
function solve_hfdmrg(H, backend::SlicedBasisBackend, psiup0; kwargs...)
    solve_hfdmrg_core(H, backend, psiup0, psiup0; restricted = true, kwargs...)
end

"""
solve_hfdmrg(H, backend::SlicedBasisBackend, psiup0, psidn0; kwargs...) -> (psiup, psidn, energy)

Unrestricted HF (UHF) sweep using a sliced-basis backend.
"""
function solve_hfdmrg(H, backend::SlicedBasisBackend, psiup0, psidn0; kwargs...)
    solve_hfdmrg_core(H, backend, psiup0, psidn0; restricted = false, kwargs...)
end

"""
solve_hfdmrg(Hup, Hdn, backend::SlicedBasisBackend, psiup0, psidn0; kwargs...) -> (psiup, psidn, energy)

Unrestricted HF (UHF) sweep with spin-dependent one-body Hamiltonians and a
sliced-basis backend. If Hup and Hdn are the same matrix object, this routes to
the common-H fast path.
"""
function solve_hfdmrg(Hup, Hdn, backend::SlicedBasisBackend, psiup0, psidn0; kwargs...)
    Hup === Hdn && return solve_hfdmrg(Hup, backend, psiup0, psidn0; kwargs...)
    solve_hfdmrg_core_split(Hup, Hdn, backend, psiup0, psidn0; kwargs...)
end

"""
solve_hfdmrg(H, layout::SliceLayout, V, psiup0; kwargs...) -> (psiup, psidn, energy)

Restricted HF (RHF) sweep using a sliced-basis backend constructed from layout
and sliced two-electron integrals (fixed V6 or ragged Vblocks).
"""
function solve_hfdmrg(H, layout::SliceLayout, V, psiup0; kwargs...)
    backend = SlicedBasisBackend(layout, V)
    solve_hfdmrg(H, backend, psiup0; kwargs...)
end

"""
solve_hfdmrg(H, layout::SliceLayout, V, psiup0, psidn0; kwargs...) -> (psiup, psidn, energy)

Unrestricted HF (UHF) sweep using a sliced-basis backend constructed from layout
and sliced two-electron integrals (fixed V6 or ragged Vblocks).
"""
function solve_hfdmrg(H, layout::SliceLayout, V, psiup0, psidn0; kwargs...)
    backend = SlicedBasisBackend(layout, V)
    solve_hfdmrg(H, backend, psiup0, psidn0; kwargs...)
end

"""
solve_hfdmrg(Hup, Hdn, layout::SliceLayout, V, psiup0, psidn0; kwargs...) -> (psiup, psidn, energy)

Unrestricted HF (UHF) sweep with spin-dependent one-body Hamiltonians using a
sliced-basis backend constructed from layout and sliced two-electron integrals.
"""
function solve_hfdmrg(Hup, Hdn, layout::SliceLayout, V, psiup0, psidn0; kwargs...)
    Hup === Hdn && return solve_hfdmrg(Hup, layout, V, psiup0, psidn0; kwargs...)
    backend = SlicedBasisBackend(layout, V)
    solve_hfdmrg_core_split(Hup, Hdn, backend, psiup0, psidn0; kwargs...)
end
