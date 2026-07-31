using LinearAlgebra

mutable struct LRBlock{TV, T}
    ra::UnitRange{Int}
    raH1::UnitRange{Int}
    raV::UnitRange{Int}
    m::Int
    phi::Matrix{T}
    H1ij::Matrix{T}
    H1phi::Matrix{T}
    vee::TV
end

struct BlockH1Cache{T}
    H1ij::Matrix{T}
    H1phi::Matrix{T}
end

mutable struct SplitLRBlock{TV, TP, THU, THD}
    ra::UnitRange{Int}
    raH1::UnitRange{Int}
    raV::UnitRange{Int}
    m::Int
    phi::Matrix{TP}
    hup::THU
    hdn::THD
    vee::TV
end

"""
State passed to a per-sweep observer.

`psiup` and `psidn` are the full physical-basis orbitals that the solver will
return if it stops after this sweep. Observers must treat them as read-only and
copy them if they need mutable retained state.
"""
struct SweepInfo{E, U, D}
    sweep::Int
    energy::E
    psiup::U
    psidn::D
    converged::Bool
end

mutable struct _RunDiagnostics
    mode::Symbol
    spectra::Vector{NamedTuple}
    windows::Vector{NamedTuple}
    eigsym_seconds::Float64
    eigsym_calls::Int
    aufbau_seconds::Float64
    aufbau_bytes::Int
    aufbau_calls::Int
    frozen_rebuilds::Int
    frozen_max_generations::Int
    frozen_live_generations::Int
    frozen_recurrence_error::Float64
    frozen_zero_up_windows::Int
    frozen_zero_dn_windows::Int
    frozen_zero_both_windows::Int
end

function _RunDiagnostics(mode::Symbol = :detailed)
    mode in (:detailed, :timing) ||
        error("_RunDiagnostics mode must be :detailed or :timing")
    _RunDiagnostics(mode, NamedTuple[], NamedTuple[], 0.0, 0, 0.0, 0, 0,
        0, 0, 0, 0.0, 0, 0, 0)
end

_detailed(diagnostics) =
    diagnostics !== nothing && diagnostics.mode === :detailed

function _record_fixed_edge!(diagnostics, side, range)
    _detailed(diagnostics) || return
    push!(diagnostics.spectra, (; sweep = 0, direction = :initialization,
        window_ordinal = 0, side, physical_range = range,
        retained_rank = length(range), status = :not_applicable,
        singular_values = nothing))
end

function _record_window!(diagnostics, sweep, direction, window_ordinal,
        left, center_range, right, local_iterations, locally_converged,
        damping_before, damping_after, rise_updates, local_energies, final_energy,
        eigsym_before)
    _detailed(diagnostics) || return
    push!(diagnostics.windows, (; sweep, direction, window_ordinal,
        left_range = left.ra, center_range, right_range = right.ra,
        left_rank = left.m, right_rank = right.m, local_iterations,
        reached_iteration_four = local_iterations == 4, locally_converged,
        cap_exhausted = local_iterations == 4 && !locally_converged,
        damping_before, damping_after, rise_updates, local_energies, final_energy,
        eigsym_seconds = diagnostics.eigsym_seconds - eigsym_before))
end

function _notify_observer(observer, info::SweepInfo)
    observer === nothing && return false
    stop = observer(info)
    (stop === nothing || stop isa Bool) || error("observer must return Bool or nothing")
    stop === true
end

# Input: ranges for A and B in absolute sites
# Output: intersection in abs. sites, as part of A range, then B, then isempty?
function intersectrange(Ara, Bra)
    firstidx = max(first(Ara), first(Bra))
    lastidx = min(last(Ara), last(Bra))
    absra = firstidx:lastidx
    ara = (absra) .- (first(Ara) - 1)
    bra = (absra) .- (first(Bra) - 1)
    absra, ara, bra, (lastidx >= firstidx)
end

function contract(A, B, transpA = false, transpB = false)
    As, Bs = size(A), size(B)
    AA = reshape(A, prod(As[1:end-1]), As[end])
    transpA && (AA = reshape(A, As[1], prod(As[2:end])))
    BB = reshape(B, Bs[1], prod(Bs[2:end]))
    transpB && (BB = reshape(B, prod(Bs[1:end-1]), Bs[end]))
    if !transpA && !transpB
        reshape(AA * BB, As[1:end-1]..., Bs[2:end]...)
    elseif transpA && !transpB
        reshape(AA' * BB, As[2:end]..., Bs[2:end]...)
    elseif !transpA && transpB
        reshape(AA * BB', As[1:end-1]..., Bs[1:end-1]...)
    else
        reshape(AA' * BB', As[2:end]..., Bs[1:end-1]...)
    end
end

function trans2(V, U)
    n = length(size(V))
    swaplast = [1:n-2..., n, n-1]
    V2 = permutedims(contract(V, U), swaplast)
    permutedims(contract(V2, U), swaplast)
end

function trans4(V, U)
    V2 = contract(U, contract(V, U), true, false)
    V2 = permutedims(V2, [2, 1, 4, 3])
    V4 = contract(U, contract(V2, U), true, false)
    permutedims(V4, [2, 1, 4, 3])
end

function transrangeleft(A, ra::UnitRange, O)
    n = size(A, 1)
    oa = O * A[ra, :]
    ra[1] == 1 && return vcat(oa, A[ra[end] + 1:n, :])
    ra[end] == n && return vcat(A[1:ra[1] - 1, :], oa)
    vcat(A[1:ra[1] - 1, :], oa, A[ra[end] + 1:n, :])
end

function transrangeleft(A, l::Integer, O)
    transrangeleft(A, l:l + size(O, 2) - 1, O)
end

function getraH1(ra, firstindsH, finalindsH)
    if ra[1] == 1
        ra[end] + 1:maximum(finalindsH[ra])
    else
        minimum(firstindsH[ra]):ra[1] - 1
    end
end

function getraV(ra, N)
    ra[1] != 1 && return 1:ra[1] - 1
    ra[end] + 1:N
end

function _check_environment_cutoff(environment_cutoff)
    environment_cutoff isa Real && isfinite(environment_cutoff) &&
        environment_cutoff >= 0 ||
        error("environment_cutoff must be a finite, nonnegative real value")
    environment_cutoff
end

function getphi(psira, diagnostics = nothing, sweep = 0, direction = :none,
        window_ordinal = 0, side = :none, physical_range = 1:0,
        environment_cutoff = 1e-10)
    u, d, v = svd(psira)
    i = findlast(x -> x > environment_cutoff, d)
    i === nothing && (i = 1)
    mkeep = max(1, i)
    if _detailed(diagnostics)
        push!(diagnostics.spectra, (; sweep, direction, window_ordinal, side,
            physical_range, retained_rank = mkeep, status = :measured,
            singular_values = copy(d)))
    end
    u[:, 1:mkeep], mkeep
end

function _makeblock(side, range, phi, H, Vee, firstindsH, finalindsH;
    dofull = true)
    N = size(H, 1)
    raH1 = getraH1(range, firstindsH, finalindsH)
    raV = getraV(range, N)
    m = size(phi, 2)
    if dofull
        H1ij = phi' * H[range, range] * phi
        H1phi = H[raH1, range] * phi
    else
        T = eltype(phi)
        H1ij = zeros(T, m, m)
        H1phi = zeros(T, length(raH1), m)
    end
    vee_state = vee_init_block(side, range, raV, phi, Vee)
    LRBlock(range, raH1, raV, m, phi, H1ij, H1phi, vee_state)
end

function getblock(side, range, psi, H, Vee, firstindsH, finalindsH; dofull = true)
    phi, _ = getphi(psi[range, :])
    _makeblock(side, range, phi, H, Vee, firstindsH, finalindsH; dofull)
end

function _h1cache(range, phi, m, H, raH1; dofull = true)
    if dofull
        return BlockH1Cache(phi' * H[range, range] * phi, H[raH1, range] * phi)
    end
    T = promote_type(eltype(phi), eltype(H))
    BlockH1Cache(zeros(T, m, m), zeros(T, length(raH1), m))
end

function _makeblock_split(side, range, phi, Hup, Hdn, Vee, firstindsH, finalindsH;
    dofull = true)
    N = size(Hup, 1)
    raH1 = getraH1(range, firstindsH, finalindsH)
    raV = getraV(range, N)
    m = size(phi, 2)
    hup = _h1cache(range, phi, m, Hup, raH1; dofull)
    hdn = _h1cache(range, phi, m, Hdn, raH1; dofull)
    vee_state = vee_init_block(side, range, raV, phi, Vee)
    SplitLRBlock(range, raH1, raV, m, phi, hup, hdn, vee_state)
end

function getblock_split(side, range, psi, Hup, Hdn, Vee, firstindsH, finalindsH;
    dofull = true)
    phi, _ = getphi(psi[range, :])
    _makeblock_split(
        side, range, phi, Hup, Hdn, Vee, firstindsH, finalindsH; dofull)
end

function _addblockright_h1(cra, PhiC, PhiR, oldcache, oldraH1, H, raH1)
    H1ij = PhiR' * oldcache.H1ij * PhiR + PhiC' * H[cra, cra] * PhiC

    if !isempty(cra) && !isempty(oldraH1)
        interabs, cr, rr, nonnull = intersectrange(cra, oldraH1)
        if nonnull
            offterm = PhiC[cr, :]' * oldcache.H1phi[rr, :] * PhiR
            H1ij += offterm + offterm'
        end
    end

    H1phi = H[raH1, cra] * PhiC
    if !isempty(raH1) && !isempty(oldraH1)
        interabs, cr, rr, nonnull = intersectrange(raH1, oldraH1)
        nonnull && (H1phi[cr, :] += oldcache.H1phi[rr, :] * PhiR)
    end

    BlockH1Cache(H1ij, H1phi)
end

function _addblockleft_h1(cra, PhiL, PhiC, oldcache, oldraH1, H, raH1)
    H1ij = PhiL' * oldcache.H1ij * PhiL + PhiC' * H[cra, cra] * PhiC

    if !isempty(cra) && !isempty(oldraH1)
        interabs, cr, lr, nonnull = intersectrange(cra, oldraH1)
        if nonnull
            offterm = PhiC[cr, :]' * oldcache.H1phi[lr, :] * PhiL
            H1ij += offterm + offterm'
        end
    end

    H1phi = H[raH1, cra] * PhiC
    if !isempty(raH1) && !isempty(oldraH1)
        interabs, cr, lr, nonnull = intersectrange(raH1, oldraH1)
        nonnull && (H1phi[cr, :] += oldcache.H1phi[lr, :] * PhiL)
    end

    BlockH1Cache(H1ij, H1phi)
end

function addblockright(cra, O, rblock, N, H, Vee, firstindsH)
    ra = cra[1]:rblock.ra[end]
    mold = rblock.m
    m = size(O, 2)
    raH1 = minimum(firstindsH[ra]):cra[1] - 1
    raV = 1:cra[1] - 1
    lc = length(cra)
    PhiC, PhiR = O[1:lc, :], O[lc + 1:end, :]
    phi = vcat(PhiC, rblock.phi * PhiR)
    Ophi = phi' * phi
    R = inv(sqrt(Symmetric(Ophi)))
    phi = phi * R
    PhiC, PhiR = PhiC * R, PhiR * R
    H1ij = PhiR' * rblock.H1ij * PhiR + PhiC' * H[cra, cra] * PhiC

    interabs, cr, rr, nonnull = intersectrange(cra, rblock.raH1)
    offterm = PhiC[cr, :]' * rblock.H1phi[rr, :] * PhiR
    H1ij += offterm + offterm'

    H1phi = H[raH1, cra] * PhiC
    interabs, cr, rr, nonnull = intersectrange(raH1, rblock.raH1)
    nonnull && (H1phi[cr, :] += rblock.H1phi[rr, :] * PhiR)

    vee_new = vee_absorb_block(:right, rblock.vee, cra, PhiR, PhiC, phi, raV, Vee)
    LRBlock(ra, raH1, raV, m, phi, H1ij, H1phi, vee_new)
end

function addblockleft(cra, O, lblock, N, H, Vee, finalindsH)
    mold = lblock.m
    m = size(O, 2)
    ra = lblock.ra[1]:cra[end]
    raH1 = cra[end] + 1:maximum(finalindsH[ra])
    raV = cra[end] + 1:N
    lc = length(cra)
    PhiC, PhiL = O[mold + 1:end, :], O[1:mold, :]
    phi = vcat(lblock.phi * PhiL, PhiC)
    Ophi = phi' * phi
    R = inv(sqrt(Symmetric(Ophi)))
    phi = phi * R
    PhiC, PhiL = PhiC * R, PhiL * R
    H1ij = PhiL' * lblock.H1ij * PhiL + PhiC' * H[cra, cra] * PhiC

    interabs, cr, lr, nonnull = intersectrange(cra, lblock.raH1)
    offterm = PhiC[cr, :]' * lblock.H1phi[lr, :] * PhiL
    H1ij += offterm + offterm'

    H1phi = H[raH1, cra] * PhiC
    interabs, cr, lr, nonnull = intersectrange(raH1, lblock.raH1)
    nonnull && (H1phi[cr, :] += lblock.H1phi[lr, :] * PhiL)

    vee_new = vee_absorb_block(:left, lblock.vee, cra, PhiL, PhiC, phi, raV, Vee)
    LRBlock(ra, raH1, raV, m, phi, H1ij, H1phi, vee_new)
end

function addblockright_split(cra, O, rblock, N, Hup, Hdn, Vee, firstindsH)
    ra = cra[1]:rblock.ra[end]
    m = size(O, 2)
    raH1 = minimum(firstindsH[ra]):cra[1] - 1
    raV = 1:cra[1] - 1
    lc = length(cra)
    PhiC, PhiR = O[1:lc, :], O[lc + 1:end, :]
    phi = vcat(PhiC, rblock.phi * PhiR)
    Ophi = phi' * phi
    R = inv(sqrt(Symmetric(Ophi)))
    phi = phi * R
    PhiC, PhiR = PhiC * R, PhiR * R

    hup = _addblockright_h1(cra, PhiC, PhiR, rblock.hup, rblock.raH1, Hup, raH1)
    hdn = _addblockright_h1(cra, PhiC, PhiR, rblock.hdn, rblock.raH1, Hdn, raH1)
    vee_new = vee_absorb_block(:right, rblock.vee, cra, PhiR, PhiC, phi, raV, Vee)
    SplitLRBlock(ra, raH1, raV, m, phi, hup, hdn, vee_new)
end

function addblockleft_split(cra, O, lblock, N, Hup, Hdn, Vee, finalindsH)
    mold = lblock.m
    m = size(O, 2)
    ra = lblock.ra[1]:cra[end]
    raH1 = cra[end] + 1:maximum(finalindsH[ra])
    raV = cra[end] + 1:N
    lc = length(cra)
    PhiC, PhiL = O[mold + 1:end, :], O[1:mold, :]
    phi = vcat(lblock.phi * PhiL, PhiC)
    Ophi = phi' * phi
    R = inv(sqrt(Symmetric(Ophi)))
    phi = phi * R
    PhiC, PhiL = PhiC * R, PhiL * R

    hup = _addblockleft_h1(cra, PhiL, PhiC, lblock.hup, lblock.raH1, Hup, raH1)
    hdn = _addblockleft_h1(cra, PhiL, PhiC, lblock.hdn, lblock.raH1, Hdn, raH1)
    vee_new = vee_absorb_block(:left, lblock.vee, cra, PhiL, PhiC, phi, raV, Vee)
    SplitLRBlock(ra, raH1, raV, m, phi, hup, hdn, vee_new)
end

function getblocksizes(N, m, blocksize; verbose = false)
    nblocks = div(N - 2 * m, blocksize) + 2
    while nblocks <= 4
        blocksize -= 1
        if blocksize == 0
            @show N, m, blocksize
        end
        nblocks = div(N - 2 * m, blocksize) + 2
    end
    aveblocksize = div(N - 2 * m, nblocks - 2)
    blocksizes = [aveblocksize for _ = 1:nblocks]
    blocksizes[1] = m
    blocksizes[end] = m
    missing = N - sum(blocksizes)
    blocksizes[2] += missing
    Cranges = [1:m for _ = 1:nblocks]
    for b = 2:nblocks
        Cranges[b] = Cranges[b - 1][end] + 1:Cranges[b - 1][end] + blocksizes[b]
    end
    verbose && @show nblocks, blocksizes
    nblocks, blocksizes, Cranges
end

function getblocksizes(N, m, blocksize, nblockcenter, block_partition, Vee;
    verbose = false)
    block_partition === nothing && return getblocksizes(N, m, blocksize; verbose)
    block_partition isa SliceLayout || error("block_partition must be a SliceLayout")
    nblockcenter >= 1 || error("block_partition requires nblockcenter >= 1")
    dims, offs = block_partition.dims, block_partition.offs
    valid = !isempty(dims) && all(>(0), dims) && offs == [0; cumsum(dims)] &&
            offs[end] == N
    valid || error("block_partition must cover exactly 1:N with nonempty contiguous ranges")
    nblocks = length(dims)
    nblocks >= 5 || error("block_partition must contain at least five blocks")
    nblocks >= nblockcenter + 4 || error("block_partition has too few blocks for nblockcenter")
    if hasproperty(Vee, :layout) && getproperty(Vee, :layout) isa SliceLayout
        layout = getproperty(Vee, :layout)
        (layout.dims == dims && layout.offs == offs) ||
            error("block_partition must exactly match the backend SliceLayout")
    end
    blocksizes = copy(dims)
    Cranges = [orb_range(block_partition, b) for b = 1:nblocks]
    verbose && @show nblocks, blocksizes
    nblocks, blocksizes, Cranges
end

function getfinal(H, N)
    function fl(x, r)
        res = findlast(x, r)
        res === nothing && (res = 1)
        res
    end
    function ff(x, r)
        res = findfirst(x, r)
        res === nothing && (res = N)
        res
    end
    finalindsH = [fl(i -> abs(H[j, i]) > 1e-13, 1:N) for j = 1:N]
    firstindsH = [ff(i -> abs(H[j, i]) > 1e-13, 1:N) for j = 1:N]
    finalindsH, firstindsH
end

function getfinal(Hup, Hdn, N)
    finalup, firstup = getfinal(Hup, N)
    finaldn, firstdn = getfinal(Hdn, N)
    max.(finalup, finaldn), min.(firstup, firstdn)
end

function getinitialblocks(nblocks, blocksizes, Cranges, psi, H, Vee, firstindsH,
        finalindsH; verbose = false, _diagnostics = nothing,
        environment_cutoff = 1e-10, _frozen = nothing)
    N = size(H, 1)
    left = 1:blocksizes[1]
    right = N - blocksizes[nblocks] + 1:N
    # Fixed edge blocks never become centers, so retain their complete physical spans.
    block1 = _makeblock(:left, left, Matrix{eltype(psi)}(I, length(left), length(left)),
        H, Vee, firstindsH, finalindsH)
    blockn = _makeblock(:right, right,
        Matrix{eltype(psi)}(I, length(right), length(right)),
        H, Vee, firstindsH, finalindsH)
    block = Vector{typeof(block1)}(undef, nblocks)
    block[1] = block1
    block[nblocks] = blockn
    _frozen === nothing ||
        _frozen_begin_blocks!(_frozen, nblocks, block1.m, blockn.m)
    _record_fixed_edge!(_diagnostics, :left, left)
    _record_fixed_edge!(_diagnostics, :right, right)
    verbose && @show block[1].ra
    for b = nblocks - 1:-1:3
        cra = Cranges[b]
        lc = length(cra)
        crra = cra[1]:N
        if _frozen === nothing
            phi, m = getphi(psi[crra, :], _diagnostics, 0, :initialization,
                nblocks - b, :right, crra, environment_cutoff)
            phi1 = block[b + 1].phi
            O = vcat(phi[1:lc, :], phi1' * phi[lc + 1:end, :])
        else
            growth = _frozen_growth!(
                _frozen, :right, block[b + 1], cra, b + 1, environment_cutoff)
            O = growth.O
        end
        block[b] = addblockright(cra, O, block[b + 1], N, H, Vee, firstindsH)
        _frozen === nothing ||
            _frozen_store_growth!(_frozen, b, b + 1, :right,
                block[b + 1], block[b], growth, _diagnostics)
    end
    _frozen === nothing || (_frozen.initializing = false)
    block
end

function getinitialblocks_split(nblocks, blocksizes, Cranges, psi, Hup, Hdn, Vee,
    firstindsH, finalindsH; verbose = false, _diagnostics = nothing,
    environment_cutoff = 1e-10, _frozen = nothing)
    N = size(Hup, 1)
    left = 1:blocksizes[1]
    right = N - blocksizes[nblocks] + 1:N
    block1 = _makeblock_split(:left, left,
        Matrix{eltype(psi)}(I, length(left), length(left)),
        Hup, Hdn, Vee, firstindsH, finalindsH)
    blockn = _makeblock_split(:right, right,
        Matrix{eltype(psi)}(I, length(right), length(right)),
        Hup, Hdn, Vee, firstindsH, finalindsH)
    block = Vector{typeof(block1)}(undef, nblocks)
    block[1] = block1
    block[nblocks] = blockn
    _frozen === nothing ||
        _frozen_begin_blocks!(_frozen, nblocks, block1.m, blockn.m)
    _record_fixed_edge!(_diagnostics, :left, left)
    _record_fixed_edge!(_diagnostics, :right, right)
    verbose && @show block[1].ra
    for b = nblocks - 1:-1:3
        cra = Cranges[b]
        lc = length(cra)
        crra = cra[1]:N
        if _frozen === nothing
            phi, m = getphi(psi[crra, :], _diagnostics, 0, :initialization,
                nblocks - b, :right, crra, environment_cutoff)
            phi1 = block[b + 1].phi
            O = vcat(phi[1:lc, :], phi1' * phi[lc + 1:end, :])
        else
            growth = _frozen_growth!(
                _frozen, :right, block[b + 1], cra, b + 1, environment_cutoff)
            O = growth.O
        end
        block[b] = addblockright_split(cra, O, block[b + 1], N, Hup, Hdn, Vee,
            firstindsH)
        _frozen === nothing ||
            _frozen_store_growth!(_frozen, b, b + 1, :right,
                block[b + 1], block[b], growth, _diagnostics)
    end
    _frozen === nothing || (_frozen.initializing = false)
    block
end

function getH1(LB::LRBlock, RB::LRBlock, H)
    Cra = LB.ra[end] + 1:RB.ra[1] - 1
    Rra = RB.ra
    lc = length(Cra)
    ml, mr = LB.m, RB.m
    n = ml + lc + mr
    lra = 1:ml
    cra = ml + 1:ml + lc
    rra = ml + lc + 1:n

    Hb = zeros(eltype(H), n, n)
    Hb[lra, lra] = LB.H1ij
    Hb[cra, cra] = H[Cra, Cra]
    Hb[rra, rra] = RB.H1ij

    interabs, lr, cr, nonnull = intersectrange(LB.raH1, Cra)
    if nonnull
        hcr = ml .+ cr
        Hb[lra, hcr] = LB.H1phi[lr, :]'
        Hb[hcr, lra] = Hb[lra, hcr]'
    end

    interabs, lr, rr, nonnull = intersectrange(LB.raH1, Rra)
    if nonnull
        Hb[lra, rra] = LB.H1phi[lr, :]' * RB.phi[rr, :]
        Hb[rra, lra] = Hb[lra, rra]'
    end

    interabs, cr, rr, nonnull = intersectrange(Cra, RB.raH1)
    if nonnull
        hcr = ml .+ cr
        Hb[hcr, rra] = RB.H1phi[rr, :]
        Hb[rra, hcr] = Hb[hcr, rra]'
    end
    Hb
end

_split_h1cache(block::SplitLRBlock, ::Val{:up}) = block.hup
_split_h1cache(block::SplitLRBlock, ::Val{:dn}) = block.hdn

function getH1(LB::SplitLRBlock, RB::SplitLRBlock, H, spin::Val)
    Lh = _split_h1cache(LB, spin)
    Rh = _split_h1cache(RB, spin)
    Cra = LB.ra[end] + 1:RB.ra[1] - 1
    Rra = RB.ra
    lc = length(Cra)
    ml, mr = LB.m, RB.m
    n = ml + lc + mr
    lra = 1:ml
    cra = ml + 1:ml + lc
    rra = ml + lc + 1:n

    T = promote_type(eltype(H), eltype(Lh.H1ij), eltype(Rh.H1ij))
    Hb = zeros(T, n, n)
    Hb[lra, lra] = Lh.H1ij
    Hb[cra, cra] = H[Cra, Cra]
    Hb[rra, rra] = Rh.H1ij

    if !isempty(LB.raH1) && !isempty(Cra)
        interabs, lr, cr, nonnull = intersectrange(LB.raH1, Cra)
        if nonnull
            hcr = ml .+ cr
            Hb[lra, hcr] = Lh.H1phi[lr, :]'
            Hb[hcr, lra] = Hb[lra, hcr]'
        end
    end

    if !isempty(LB.raH1) && !isempty(Rra)
        interabs, lr, rr, nonnull = intersectrange(LB.raH1, Rra)
        if nonnull
            Hb[lra, rra] = Lh.H1phi[lr, :]' * RB.phi[rr, :]
            Hb[rra, lra] = Hb[lra, rra]'
        end
    end

    if !isempty(Cra) && !isempty(RB.raH1)
        interabs, cr, rr, nonnull = intersectrange(Cra, RB.raH1)
        if nonnull
            hcr = ml .+ cr
            Hb[hcr, rra] = Rh.H1phi[rr, :]
            Hb[rra, hcr] = Hb[hcr, rra]'
        end
    end
    Hb
end

function getpsiexpanded(psiup, Lra, Lphi, Cra, Rra, Rphi)
    ml = size(Lphi, 2)
    mr = size(Rphi, 2)
    lra = 1:ml
    lc = length(Cra)
    cra = ml + 1:ml + lc
    rra = ml + lc + 1:ml + lc + mr
    psiL = Lphi * psiup[lra, :]
    psiC = psiup[cra, :]
    psiR = Rphi * psiup[rra, :]
    vcat(psiL, psiC, psiR)
end

function getpsireduced(psiup, Lra, Lphi, Cra, Rra, Rphi)
    psiL = Lphi' * psiup[Lra, :]
    psiC = psiup[Cra, :]
    psiR = Rphi' * psiup[Rra, :]
    vcat(psiL, psiC, psiR)
end

function eigsym(A, diagnostics = nothing)
    start_ns = diagnostics === nothing ? 0 : time_ns()
    E = eigen(Symmetric(A))
    diagnostics === nothing || (diagnostics.eigsym_calls += 1)
    diagnostics === nothing ||
        (diagnostics.eigsym_seconds += (time_ns() - start_ns) * 1.0e-9)
    E.values, E.vectors
end

function _rel_converged(energy_old, energy_new, cutoff)
    scale = max(1.0, abs(energy_new))
    abs(energy_old - energy_new) < cutoff * scale
end

function solve_hfdmrg_core(H, Vee, psiup0, psidn0; frozen_occupied = :off, kwargs...)
    frozen_occupied in (:off, :roundoff_exact) ||
        error("frozen_occupied must be :off or :roundoff_exact")
    _solve_hfdmrg_core(Val(frozen_occupied), H, Vee, psiup0, psidn0; kwargs...)
end

function _solve_hfdmrg_core(::Val{frozen_mode}, H, Vee, psiup0, psidn0;
    restricted = false,
    nblockcenter = 1,
    blocksize = 200,
    block_partition = nothing,
    maxiter = 1000,
    cutoff = 1e-11,
    scf_cutoff = nothing,
    environment_cutoff = 1e-10,
    observer = nothing,
    _diagnostics = nothing,
    verbose = false) where {frozen_mode}

    environment_cutoff = _check_environment_cutoff(environment_cutoff)
    Nup, Ndn, N = size(psiup0, 2), size(psidn0, 2), size(H, 1)
    frozen = frozen_mode === :off ? nothing :
        _frozen_policy(H, H, Vee, psiup0, psidn0, restricted)
    finalindsH, firstindsH = getfinal(H, N)
    m = restricted ? Nup : Nup + Ndn
    psiall = restricted ? psiup0 : hcat(psiup0, psidn0)
    nblocks, blocksizes, Cranges = getblocksizes(N, m, blocksize, nblockcenter,
        block_partition, Vee; verbose)
    block = getinitialblocks(nblocks, blocksizes, Cranges, psiall, H, Vee,
        firstindsH, finalindsH; verbose, _diagnostics, environment_cutoff,
        _frozen = frozen)

    Lra = block[1].ra
    Rra = block[2 + nblockcenter].ra
    Cra = Lra[end] + 1:Rra[1] - 1
    psiup = getpsireduced(psiup0, Lra, block[1].phi, Cra, Rra, block[2 + nblockcenter].phi)
    psidn = restricted ? psiup : getpsireduced(psidn0, Lra, block[1].phi, Cra, Rra,
        block[2 + nblockcenter].phi)

    left2right = [(b, 1, 1) for b = 1:nblocks - 1 - nblockcenter]
    left2right[end] = (nblocks - 1 - nblockcenter, 1, -1)
    right2left = [(b, -1, -1) for b = nblocks - 2 - nblockcenter:-1:2]
    fullsweep = vcat(left2right, right2left)

    psiallup = psiup0
    psialldn = restricted ? psiup0 : psidn0
    energyiter = 10000.0
    lambda = [1.0 for _ = 1:nblocks]
    scf_cutoff === nothing && (scf_cutoff = cutoff)

    for iter = 1:maxiter
        energy = 0.0
        restart_frozen = false
        verbose && println()
        verbose && @show iter
        for (b, dir, transdir) in fullsweep
            detailed = _detailed(_diagnostics)
            direction = detailed ? (dir == 1 ? :left_to_right : :right_to_left) : :none
            window_ordinal = detailed ?
                (dir == 1 ? b : nblocks - 1 - nblockcenter - b) : 0
            H1B = getH1(block[b], block[b + 1 + nblockcenter], H)
            Cra = block[b].ra[end] + 1:block[b + 1 + nblockcenter].ra[1] - 1
            win = vee_window(block[b].vee, block[b + 1 + nblockcenter].vee, Cra, Vee)
            energy = energylast = 1e10
            local_energies = detailed ? Float64[] : nothing
            rise_updates = detailed ? Int[] : nothing
            damping_before = detailed ? lambda[b] : 0.0
            eigsym_before = _diagnostics === nothing ? 0.0 : _diagnostics.eigsym_seconds
            local_iterations, locally_converged = 0, false
            if frozen === nothing
            psiup_occ = @view psiup[:, 1:Nup]
            rhoup = psiup_occ * psiup_occ'
            if !restricted
                psidn_occ = @view psidn[:, 1:Ndn]
                rhodn = psidn_occ * psidn_occ'
            end
            Fup = copy(H1B)
            if restricted
                vee_add_fock_r!(Fup, rhoup, win)
            else
                Fdn = copy(H1B)
                vee_add_fock!(Fup, Fdn, rhoup, rhodn, win)
            end
            for s = 1:4
                if restricted
                    evals, evecs = eigsym(Fup, _diagnostics)
                    psiup = evecs[:, 1:Nup]
                    rhoup = (1 - lambda[b]) * rhoup + lambda[b] * psiup * psiup'
                    Fup = copy(H1B)
                    vee_add_fock_r!(Fup, rhoup, win)
                    energy = tr(rhoup * (Fup + H1B))
                else
                    evals, evecs = eigsym(Fup, _diagnostics)
                    psiup = evecs[:, 1:Nup]
                    rhoup = (1 - lambda[b]) * rhoup + lambda[b] * psiup * psiup'
                    evals, evecs = eigsym(Fdn, _diagnostics)
                    psidn = evecs[:, 1:Ndn]
                    rhodn = (1 - lambda[b]) * rhodn + lambda[b] * psidn * psidn'
                    Fup, Fdn = copy(H1B), copy(H1B)
                    vee_add_fock!(Fup, Fdn, rhoup, rhodn, win)
                    energy = 0.5 * tr(rhoup * (Fup + H1B)) + 0.5 * tr(rhodn * (Fdn + H1B))
                end
                local_iterations = s
                detailed && push!(local_energies, energy)
                if energy > energylast
                    detailed && push!(rise_updates, s)
                    lambda[b] *= 0.5
                end
                locally_converged = abs(energylast - energy) < scf_cutoff
                (locally_converged || s == 4) && break
                energylast = energy
            end
            else
                result = _frozen_local!(frozen, b, b + 1 + nblockcenter, block[b],
                    Cra, block[b + 1 + nblockcenter], H1B, H1B, win, lambda[b],
                    scf_cutoff, false, _diagnostics)
                if result.stale
                    _protect_stale!(frozen, result.up, result.dn, _diagnostics)
                    restart_frozen = true
                    break
                end
                psiup, psidn, energy = result.psiup, result.psidn, result.energy
                local_iterations, locally_converged = result.nit, result.converged
                lambda[b] = result.lambda
                detailed && (local_energies = result.energies; rise_updates = result.rises)
            end

            _record_window!(_diagnostics, iter, direction, window_ordinal, block[b],
                Cra, block[b + 1 + nblockcenter], local_iterations,
                locally_converged, damping_before, lambda[b], rise_updates,
                local_energies, energy, eigsym_before)

            rbl = block[b + 1 + nblockcenter]
            if frozen === nothing
                psiallup = getpsiexpanded(psiup, block[b].ra, block[b].phi,
                    block[b].ra[end] + 1:rbl.ra[1] - 1, rbl.ra, rbl.phi)
                !restricted && (psialldn = getpsiexpanded(psidn, block[b].ra,
                    block[b].phi, block[b].ra[end] + 1:rbl.ra[1] - 1,
                    rbl.ra, rbl.phi))
                psi = restricted ? psiup : hcat(psiup, psidn)
            else
                psiallup, psialldn = frozen.psiup, frozen.psidn
            end
            if frozen !== nothing
                if transdir == 1
                    growth = _frozen_growth!(frozen, :left, block[b],
                        Cranges[b + 1], b, environment_cutoff)
                    block[b + 1] = addblockleft(
                        Cranges[b + 1], growth.O, block[b], N, H, Vee, finalindsH)
                    _frozen_store_growth!(frozen, b + 1, b, :left,
                        block[b], block[b + 1], growth, _diagnostics)
                else
                    d = nblockcenter
                    growth = _frozen_growth!(frozen, :right, block[b + 1 + d],
                        Cranges[b + d], b + 1 + d, environment_cutoff)
                    block[b + d] = addblockright(Cranges[b + d], growth.O,
                        block[b + 1 + d], N, H, Vee, firstindsH)
                    _frozen_store_growth!(frozen, b + d, b + 1 + d, :right,
                        block[b + 1 + d], block[b + d], growth, _diagnostics)
                end
            elseif nblockcenter < 1
                Cra = block[b].ra[end] + 1:block[b + 1 + nblockcenter].ra[1] - 1
                oldlc = length(Cra)
                if transdir == 1
                    oldm = block[b].m
                    O, m = getphi(psi[1:oldm + oldlc, :], _diagnostics, iter,
                        direction, window_ordinal, :left,
                        block[b].ra[1]:last(Cra), environment_cutoff)
                    block[b + 1] = addblockleft(Cra, O, block[b], N, H, Vee, finalindsH)
                    lc = length(Cranges[b + 2])
                    psibb1 = transrangeleft(psi, 1, O')
                    psiR = transrangeleft(psibb1, m + 1, block[b + 2].phi)
                    psi = transrangeleft(psiR, m + lc + 1, block[b + 3].phi')
                else
                    oldm = block[b].m
                    O, m = getphi(psi[oldm + 1:end, :], _diagnostics, iter,
                        direction, window_ordinal, :right,
                        first(Cra):block[b + 2].ra[end], environment_cutoff)
                    block[b + 1] = addblockright(Cra, O, block[b + 2], N, H, Vee, firstindsH)
                    psiL = transrangeleft(psi, m + 1, O')
                    psiLL = transrangeleft(psiL, 1:oldm, block[b].phi)
                    psi = transrangeleft(psiLL, 1, block[b - 1].phi')
                end
            else
                d = nblockcenter
                if transdir == 1
                    oldm = block[b].m
                    O, m = getphi(psi[1:oldm + blocksizes[b + 1], :],
                        _diagnostics, iter, direction, window_ordinal, :left,
                        block[b].ra[1]:Cranges[b + 1][end], environment_cutoff)
                    block[b + 1] = addblockleft(Cranges[b + 1], O, block[b], N, H, Vee, finalindsH)
                    Cra = Cranges[b + 2][1]:Cranges[b + 1 + d][end]
                    lc = length(Cra)
                    psibb1 = transrangeleft(psiup, 1, O')
                    posright = m + sum(@view blocksizes[b + 2:b + d]) + 1
                    psiR = transrangeleft(psibb1, posright, block[b + 1 + d].phi)
                    psiup = transrangeleft(psiR, m + lc + 1, block[b + 2 + d].phi')
                    if !restricted
                        psibb1 = transrangeleft(psidn, 1, O')
                        psiR = transrangeleft(psibb1, posright, block[b + 1 + d].phi)
                        psidn = transrangeleft(psiR, m + lc + 1, block[b + 2 + d].phi')
                    end
                else
                    startO = size(psi, 1) - block[b + 1 + d].m - blocksizes[b + d] + 1
                    oldm = block[b].m
                    O, m = getphi(psi[startO:end, :], _diagnostics, iter,
                        direction, window_ordinal, :right,
                        Cranges[b + d][1]:block[b + 1 + d].ra[end],
                        environment_cutoff)
                    block[b + d] = addblockright(Cranges[b + d], O, block[b + 1 + d], N, H, Vee, firstindsH)
                    psiL = transrangeleft(psiup, startO, O')
                    psiLL = transrangeleft(psiL, 1:oldm, block[b].phi)
                    psiup = transrangeleft(psiLL, 1, block[b - 1].phi')
                    if !restricted
                        psiL = transrangeleft(psidn, startO, O')
                        psiLL = transrangeleft(psiL, 1:oldm, block[b].phi)
                        psidn = transrangeleft(psiLL, 1, block[b - 1].phi')
                    end
                end
            end
        end
        if restart_frozen
            energyiter = _frozen_full_energy(frozen)
            block = getinitialblocks(nblocks, blocksizes, Cranges,
                restricted ? frozen.psiup : hcat(frozen.psiup, frozen.psidn),
                H, Vee, firstindsH, finalindsH; verbose, _diagnostics,
                environment_cutoff, _frozen = frozen)
            psiallup, psialldn = frozen.psiup, frozen.psidn
            stop_requested = _notify_observer(observer, SweepInfo(iter, energyiter,
                psiallup, restricted ? psiallup : psialldn, false))
            stop_requested && break
            continue
        end
        verbose && @show iter, energy, energyiter
        converged = abs(energyiter - energy) < cutoff
        energyiter = energy
        if frozen !== nothing
            rebuild, converged = _frozen_audit!(frozen, converged, _diagnostics)
            if rebuild
                block = getinitialblocks(nblocks, blocksizes, Cranges,
                    restricted ? frozen.psiup : hcat(frozen.psiup, frozen.psidn),
                    H, Vee, firstindsH, finalindsH; verbose, _diagnostics,
                    environment_cutoff, _frozen = frozen)
            end
            psiallup, psialldn = frozen.psiup, frozen.psidn
        end
        info = SweepInfo(iter, energyiter, psiallup,
            restricted ? psiallup : psialldn, converged)
        stop_requested = _notify_observer(observer, info)
        (converged || stop_requested) && break
    end

    restricted && (psialldn = psiallup)
    psiallup, psialldn, energyiter
end

function solve_hfdmrg_core_split(Hup, Hdn, Vee, psiup0, psidn0;
        frozen_occupied = :off, kwargs...)
    frozen_occupied in (:off, :roundoff_exact) ||
        error("frozen_occupied must be :off or :roundoff_exact")
    _solve_hfdmrg_core_split(
        Val(frozen_occupied), Hup, Hdn, Vee, psiup0, psidn0; kwargs...)
end

function _solve_hfdmrg_core_split(::Val{frozen_mode}, Hup, Hdn, Vee, psiup0, psidn0;
    nblockcenter = 1,
    blocksize = 200,
    block_partition = nothing,
    maxiter = 1000,
    cutoff = 1e-11,
    scf_cutoff = nothing,
    environment_cutoff = 1e-10,
    observer = nothing,
    _diagnostics = nothing,
    verbose = false) where {frozen_mode}

    environment_cutoff = _check_environment_cutoff(environment_cutoff)
    Nup, Ndn, N = size(psiup0, 2), size(psidn0, 2), size(Hup, 1)
    frozen = frozen_mode === :off ? nothing :
        _frozen_policy(Hup, Hdn, Vee, psiup0, psidn0, false)
    size(Hup, 2) == N || error("Hup must be square")
    size(Hdn) == (N, N) || error("Hdn must have the same size as Hup")
    size(psiup0, 1) == N || error("psiup0 has wrong row dimension")
    size(psidn0, 1) == N || error("psidn0 has wrong row dimension")

    finalindsH, firstindsH = getfinal(Hup, Hdn, N)
    m = Nup + Ndn
    psiall = hcat(psiup0, psidn0)
    nblocks, blocksizes, Cranges = getblocksizes(N, m, blocksize, nblockcenter,
        block_partition, Vee; verbose)
    block = getinitialblocks_split(nblocks, blocksizes, Cranges, psiall, Hup, Hdn,
        Vee, firstindsH, finalindsH; verbose, _diagnostics, environment_cutoff,
        _frozen = frozen)

    Lra = block[1].ra
    Rra = block[2 + nblockcenter].ra
    Cra = Lra[end] + 1:Rra[1] - 1
    psiup = getpsireduced(psiup0, Lra, block[1].phi, Cra, Rra,
        block[2 + nblockcenter].phi)
    psidn = getpsireduced(psidn0, Lra, block[1].phi, Cra, Rra,
        block[2 + nblockcenter].phi)

    left2right = [(b, 1, 1) for b = 1:nblocks - 1 - nblockcenter]
    left2right[end] = (nblocks - 1 - nblockcenter, 1, -1)
    right2left = [(b, -1, -1) for b = nblocks - 2 - nblockcenter:-1:2]
    fullsweep = vcat(left2right, right2left)

    psiallup = psiup0
    psialldn = psidn0
    energyiter = 10000.0
    lambda = [1.0 for _ = 1:nblocks]
    scf_cutoff === nothing && (scf_cutoff = cutoff / 10)

    for iter = 1:maxiter
        energy = 0.0
        restart_frozen = false
        verbose && println()
        verbose && @show iter
        verbose && flush(stdout)
        for (b, dir, transdir) in fullsweep
            detailed = _detailed(_diagnostics)
            direction = detailed ? (dir == 1 ? :left_to_right : :right_to_left) : :none
            window_ordinal = detailed ?
                (dir == 1 ? b : nblocks - 1 - nblockcenter - b) : 0
            H1Bup = getH1(block[b], block[b + 1 + nblockcenter], Hup, Val(:up))
            H1Bdn = getH1(block[b], block[b + 1 + nblockcenter], Hdn, Val(:dn))
            Cra = block[b].ra[end] + 1:block[b + 1 + nblockcenter].ra[1] - 1
            win = vee_window(block[b].vee, block[b + 1 + nblockcenter].vee, Cra, Vee)
            energy = energylast = 1e10
            local_energies = detailed ? Float64[] : nothing
            rise_updates = detailed ? Int[] : nothing
            damping_before = detailed ? lambda[b] : 0.0
            eigsym_before = _diagnostics === nothing ? 0.0 : _diagnostics.eigsym_seconds
            local_iterations, locally_converged = 0, false
            if frozen === nothing
            psiup_occ = @view psiup[:, 1:Nup]
            psidn_occ = @view psidn[:, 1:Ndn]
            rhoup = psiup_occ * psiup_occ'
            rhodn = psidn_occ * psidn_occ'
            Fup, Fdn = copy(H1Bup), copy(H1Bdn)
            vee_add_fock!(Fup, Fdn, rhoup, rhodn, win)
            for s = 1:4
                evals, evecs = eigsym(Fup, _diagnostics)
                psiup = evecs[:, 1:Nup]
                rhoup = (1 - lambda[b]) * rhoup + lambda[b] * psiup * psiup'
                evals, evecs = eigsym(Fdn, _diagnostics)
                psidn = evecs[:, 1:Ndn]
                rhodn = (1 - lambda[b]) * rhodn + lambda[b] * psidn * psidn'
                Fup, Fdn = copy(H1Bup), copy(H1Bdn)
                vee_add_fock!(Fup, Fdn, rhoup, rhodn, win)
                energy = 0.5 * tr(rhoup * (Fup + H1Bup)) +
                         0.5 * tr(rhodn * (Fdn + H1Bdn))

                local_iterations = s
                detailed && push!(local_energies, energy)
                if energy > energylast
                    detailed && push!(rise_updates, s)
                    lambda[b] *= 0.5
                end
                locally_converged = _rel_converged(energylast, energy, scf_cutoff)
                (locally_converged || s == 4) && break
                energylast = energy
            end
            else
                result = _frozen_local!(frozen, b, b + 1 + nblockcenter, block[b],
                    Cra, block[b + 1 + nblockcenter], H1Bup, H1Bdn, win, lambda[b],
                    scf_cutoff, true, _diagnostics)
                if result.stale
                    _protect_stale!(frozen, result.up, result.dn, _diagnostics)
                    restart_frozen = true
                    break
                end
                psiup, psidn, energy = result.psiup, result.psidn, result.energy
                local_iterations, locally_converged = result.nit, result.converged
                lambda[b] = result.lambda
                detailed && (local_energies = result.energies; rise_updates = result.rises)
            end

            _record_window!(_diagnostics, iter, direction, window_ordinal, block[b],
                Cra, block[b + 1 + nblockcenter], local_iterations,
                locally_converged, damping_before, lambda[b], rise_updates,
                local_energies, energy, eigsym_before)

            rbl = block[b + 1 + nblockcenter]
            if frozen === nothing
                psiallup = getpsiexpanded(psiup, block[b].ra, block[b].phi,
                    block[b].ra[end] + 1:rbl.ra[1] - 1, rbl.ra, rbl.phi)
                psialldn = getpsiexpanded(psidn, block[b].ra, block[b].phi,
                    block[b].ra[end] + 1:rbl.ra[1] - 1, rbl.ra, rbl.phi)
                psi = hcat(psiup, psidn)
            else
                psiallup, psialldn = frozen.psiup, frozen.psidn
            end
            if frozen !== nothing
                if transdir == 1
                    growth = _frozen_growth!(frozen, :left, block[b],
                        Cranges[b + 1], b, environment_cutoff)
                    block[b + 1] = addblockleft_split(Cranges[b + 1], growth.O,
                        block[b], N, Hup, Hdn, Vee, finalindsH)
                    _frozen_store_growth!(frozen, b + 1, b, :left,
                        block[b], block[b + 1], growth, _diagnostics)
                else
                    d = nblockcenter
                    growth = _frozen_growth!(frozen, :right, block[b + 1 + d],
                        Cranges[b + d], b + 1 + d, environment_cutoff)
                    block[b + d] = addblockright_split(Cranges[b + d], growth.O,
                        block[b + 1 + d], N, Hup, Hdn, Vee, firstindsH)
                    _frozen_store_growth!(frozen, b + d, b + 1 + d, :right,
                        block[b + 1 + d], block[b + d], growth, _diagnostics)
                end
            elseif nblockcenter < 1
                Cra = block[b].ra[end] + 1:block[b + 1 + nblockcenter].ra[1] - 1
                oldlc = length(Cra)
                if transdir == 1
                    oldm = block[b].m
                    O, m = getphi(psi[1:oldm + oldlc, :], _diagnostics, iter,
                        direction, window_ordinal, :left,
                        block[b].ra[1]:last(Cra), environment_cutoff)
                    block[b + 1] = addblockleft_split(Cra, O, block[b], N, Hup, Hdn,
                        Vee, finalindsH)
                    lc = length(Cranges[b + 2])
                    psibb1 = transrangeleft(psi, 1, O')
                    psiR = transrangeleft(psibb1, m + 1, block[b + 2].phi)
                    psi = transrangeleft(psiR, m + lc + 1, block[b + 3].phi')
                else
                    oldm = block[b].m
                    O, m = getphi(psi[oldm + 1:end, :], _diagnostics, iter,
                        direction, window_ordinal, :right,
                        first(Cra):block[b + 2].ra[end], environment_cutoff)
                    block[b + 1] = addblockright_split(Cra, O, block[b + 2], N, Hup,
                        Hdn, Vee, firstindsH)
                    psiL = transrangeleft(psi, m + 1, O')
                    psiLL = transrangeleft(psiL, 1:oldm, block[b].phi)
                    psi = transrangeleft(psiLL, 1, block[b - 1].phi')
                end
                psiup = psi[:, 1:Nup]
                psidn = psi[:, Nup + 1:Nup + Ndn]
            else
                d = nblockcenter
                if transdir == 1
                    oldm = block[b].m
                    O, m = getphi(psi[1:oldm + blocksizes[b + 1], :],
                        _diagnostics, iter, direction, window_ordinal, :left,
                        block[b].ra[1]:Cranges[b + 1][end], environment_cutoff)
                    block[b + 1] = addblockleft_split(Cranges[b + 1], O, block[b], N,
                        Hup, Hdn, Vee, finalindsH)
                    Cra = Cranges[b + 2][1]:Cranges[b + 1 + d][end]
                    lc = length(Cra)
                    psibb1 = transrangeleft(psiup, 1, O')
                    posright = m + sum(@view blocksizes[b + 2:b + d]) + 1
                    psiR = transrangeleft(psibb1, posright, block[b + 1 + d].phi)
                    psiup = transrangeleft(psiR, m + lc + 1, block[b + 2 + d].phi')
                    psibb1 = transrangeleft(psidn, 1, O')
                    psiR = transrangeleft(psibb1, posright, block[b + 1 + d].phi)
                    psidn = transrangeleft(psiR, m + lc + 1, block[b + 2 + d].phi')
                else
                    startO = size(psi, 1) - block[b + 1 + d].m - blocksizes[b + d] + 1
                    oldm = block[b].m
                    O, m = getphi(psi[startO:end, :], _diagnostics, iter,
                        direction, window_ordinal, :right,
                        Cranges[b + d][1]:block[b + 1 + d].ra[end],
                        environment_cutoff)
                    block[b + d] = addblockright_split(Cranges[b + d], O,
                        block[b + 1 + d], N, Hup, Hdn, Vee, firstindsH)
                    psiL = transrangeleft(psiup, startO, O')
                    psiLL = transrangeleft(psiL, 1:oldm, block[b].phi)
                    psiup = transrangeleft(psiLL, 1, block[b - 1].phi')
                    psiL = transrangeleft(psidn, startO, O')
                    psiLL = transrangeleft(psiL, 1:oldm, block[b].phi)
                    psidn = transrangeleft(psiLL, 1, block[b - 1].phi')
                end
            end
        end
        if restart_frozen
            energyiter = _frozen_full_energy(frozen)
            block = getinitialblocks_split(nblocks, blocksizes, Cranges,
                hcat(frozen.psiup, frozen.psidn), Hup, Hdn, Vee,
                firstindsH, finalindsH; verbose, _diagnostics,
                environment_cutoff, _frozen = frozen)
            psiallup, psialldn = frozen.psiup, frozen.psidn
            stop_requested = _notify_observer(observer,
                SweepInfo(iter, energyiter, psiallup, psialldn, false))
            stop_requested && break
            continue
        end
        verbose && @show iter, energy, energyiter
        verbose && flush(stdout)
        converged = _rel_converged(energyiter, energy, cutoff)
        energyiter = energy
        if frozen !== nothing
            rebuild, converged = _frozen_audit!(frozen, converged, _diagnostics)
            if rebuild
                block = getinitialblocks_split(nblocks, blocksizes, Cranges,
                    hcat(frozen.psiup, frozen.psidn), Hup, Hdn, Vee,
                    firstindsH, finalindsH; verbose, _diagnostics,
                    environment_cutoff, _frozen = frozen)
            end
            psiallup, psialldn = frozen.psiup, frozen.psidn
        end
        info = SweepInfo(iter, energyiter, psiallup, psialldn, converged)
        stop_requested = _notify_observer(observer, info)
        (converged || stop_requested) && break
    end

    psiallup, psialldn, energyiter
end

"""
solve_hfdmrg(H, V, psiup0; kwargs...) -> (psiup, psidn, energy)

Restricted HF (RHF) sweep using the density-density backend. H and V are N x N
matrices. psiup0 is N x Nup with orthonormal columns. Returns the optimized
orbitals in the physical site basis for the last sweep position and the final
energy. For RHF, psidn == psiup.
"""
function solve_hfdmrg(H, V, psiup0; kwargs...)
    backend = DensityDensityBackend(V)
    haskey(kwargs, :frozen_occupied) &&
        return solve_hfdmrg_core(H, backend, psiup0, psiup0; restricted = true, kwargs...)
    _solve_hfdmrg_core(Val(:off), H, backend, psiup0, psiup0; restricted = true, kwargs...)
end

"""
solve_hfdmrg(H, V, psiup0, psidn0; kwargs...) -> (psiup, psidn, energy)

Unrestricted HF (UHF) sweep using the density-density backend. H and V are N x N
matrices. psiup0/psidn0 are N x Nup/Ndn with orthonormal columns. Returns the
optimized orbitals in the physical site basis for the last sweep position and
the final energy.
"""
function solve_hfdmrg(H, V, psiup0, psidn0; kwargs...)
    backend = DensityDensityBackend(V)
    haskey(kwargs, :frozen_occupied) &&
        return solve_hfdmrg_core(H, backend, psiup0, psidn0; restricted = false, kwargs...)
    _solve_hfdmrg_core(Val(:off), H, backend, psiup0, psidn0; restricted = false, kwargs...)
end

"""
solve_hfdmrg(Hup, Hdn, V, psiup0, psidn0; kwargs...) -> (psiup, psidn, energy)

Unrestricted HF (UHF) sweep with spin-dependent one-body Hamiltonians and the
density-density backend. The two-body interaction V is unchanged. If Hup and
Hdn are the same matrix object, this routes to the common-H fast path.
"""
function solve_hfdmrg(Hup, Hdn, V, psiup0, psidn0; kwargs...)
    Hup === Hdn && return solve_hfdmrg(Hup, V, psiup0, psidn0; kwargs...)
    backend = DensityDensityBackend(V)
    haskey(kwargs, :frozen_occupied) &&
        return solve_hfdmrg_core_split(Hup, Hdn, backend, psiup0, psidn0; kwargs...)
    _solve_hfdmrg_core_split(Val(:off), Hup, Hdn, backend, psiup0, psidn0; kwargs...)
end
