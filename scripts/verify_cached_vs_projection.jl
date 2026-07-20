using LinearAlgebra
using Random
using HFDMRG

function orthonormal_cols(rng, n, k)
    q = qr(randn(rng, n, k)).Q
    Matrix(q[:, 1:k])
end

function build_windows(layout, V, Lra, Cra, Rra, Lphi, Rphi)
    backend_proj = HFDMRG.SlicedBasisBackend(layout, V)
    backend_cached = HFDMRG.SlicedBasisBackendCached(layout, V)
    N = layout.offs[end]

    Lvee_proj = HFDMRG.vee_init_block(:left, Lra, (Lra[end] + 1):N, Lphi, backend_proj)
    Rvee_proj = HFDMRG.vee_init_block(:right, Rra, 1:(Rra[1] - 1), Rphi, backend_proj)
    win_proj = HFDMRG.vee_window(Lvee_proj, Rvee_proj, Cra, backend_proj)

    Lvee_cached = HFDMRG.vee_init_block(:left, Lra, (Lra[end] + 1):N, Lphi, backend_cached)
    Rvee_cached = HFDMRG.vee_init_block(:right, Rra, 1:(Rra[1] - 1), Rphi, backend_cached)
    win_cached = HFDMRG.vee_window(Lvee_cached, Rvee_cached, Cra, backend_cached)

    win_proj, win_cached
end

function compare_window_rhf(layout, V, Lra, Cra, Rra, Lphi, Rphi, rho)
    win_proj, win_cached = build_windows(layout, V, Lra, Cra, Rra, Lphi, Rphi)
    F_proj = zeros(size(rho))
    F_cached = zeros(size(rho))
    HFDMRG.vee_add_fock_r!(F_proj, rho, win_proj)
    HFDMRG.vee_add_fock_r!(F_cached, rho, win_cached)
    maximum(abs.(F_proj .- F_cached))
end

function compare_window_uhf(layout, V, Lra, Cra, Rra, Lphi, Rphi, rhoup, rhodn)
    win_proj, win_cached = build_windows(layout, V, Lra, Cra, Rra, Lphi, Rphi)
    Fup_proj = zeros(size(rhoup))
    Fdn_proj = zeros(size(rhoup))
    Fup_cached = zeros(size(rhoup))
    Fdn_cached = zeros(size(rhoup))
    HFDMRG.vee_add_fock!(Fup_proj, Fdn_proj, rhoup, rhodn, win_proj)
    HFDMRG.vee_add_fock!(Fup_cached, Fdn_cached, rhoup, rhodn, win_cached)
    (maxdiff_up = maximum(abs.(Fup_proj .- Fup_cached)),
     maxdiff_dn = maximum(abs.(Fdn_proj .- Fdn_cached)))
end

function compare_sweep_rhf(layout, V, H, psiup0; kwargs...)
    backend_proj = HFDMRG.SlicedBasisBackend(layout, V)
    backend_cached = HFDMRG.SlicedBasisBackendCached(layout, V)
    _, _, e_proj = HFDMRG.solve_hfdmrg(H, backend_proj, psiup0; kwargs...)
    _, _, e_cached = HFDMRG.solve_hfdmrg(H, backend_cached, psiup0; kwargs...)
    (energy_proj = e_proj, energy_cached = e_cached, energy_diff = abs(e_proj - e_cached))
end

rng = MersenneTwister(1234)

println("Fixed-size (nj=2, ns=5)")
ns = 5
nj = 2
layout = HFDMRG.SliceLayout(fill(nj, ns))
N = layout.offs[end]
V6 = randn(rng, nj, nj, nj, nj, ns, ns)
H = randn(rng, N, N)
H = (H + H') / 2

Lra = 1:2
Cra = 3:6
Rra = 7:10
Lphi = orthonormal_cols(rng, length(Lra), 1)
Rphi = orthonormal_cols(rng, length(Rra), 1)

superdim = size(Lphi, 2) + length(Cra) + size(Rphi, 2)
rho = randn(rng, superdim, superdim)
rho = (rho + rho') / 2
rhoup = randn(rng, superdim, superdim)
rhoup = (rhoup + rhoup') / 2
rhodn = randn(rng, superdim, superdim)
rhodn = (rhodn + rhodn') / 2

res_rhf = compare_window_rhf(layout, V6, Lra, Cra, Rra, Lphi, Rphi, rho)
println("  window RHF max |ΔF| = ", res_rhf)
res_uhf = compare_window_uhf(layout, V6, Lra, Cra, Rra, Lphi, Rphi, rhoup, rhodn)
println("  window UHF max |ΔFup| = ", res_uhf.maxdiff_up,
    ", max |ΔFdn| = ", res_uhf.maxdiff_dn)

psiup0 = orthonormal_cols(rng, N, 2)
res_sweep = compare_sweep_rhf(layout, V6, H, psiup0;
    maxiter = 2, blocksize = 2, block_partition = layout, cutoff = 1e-8, verbose = false)
println("  sweep RHF energy proj = ", res_sweep.energy_proj,
    ", cached = ", res_sweep.energy_cached,
    ", |ΔE| = ", res_sweep.energy_diff)

println("Ragged dims [1,2,3,2,2]")
dims = [1, 2, 3, 2, 2]
layout_r = HFDMRG.SliceLayout(dims)
ns_r = length(dims)
N_r = layout_r.offs[end]
Vblocks = [[randn(rng, dims[n], dims[n], dims[m], dims[m]) for m in 1:ns_r]
           for n in 1:ns_r]
H_r = randn(rng, N_r, N_r)
H_r = (H_r + H_r') / 2

Lra = 1:3
Cra = 4:6
Rra = 7:10
Lphi = orthonormal_cols(rng, length(Lra), 1)
Rphi = orthonormal_cols(rng, length(Rra), 1)

superdim = size(Lphi, 2) + length(Cra) + size(Rphi, 2)
rho = randn(rng, superdim, superdim)
rho = (rho + rho') / 2
rhoup = randn(rng, superdim, superdim)
rhoup = (rhoup + rhoup') / 2
rhodn = randn(rng, superdim, superdim)
rhodn = (rhodn + rhodn') / 2

res_rhf = compare_window_rhf(layout_r, Vblocks, Lra, Cra, Rra, Lphi, Rphi, rho)
println("  window RHF max |ΔF| = ", res_rhf)
res_uhf = compare_window_uhf(layout_r, Vblocks, Lra, Cra, Rra, Lphi, Rphi, rhoup, rhodn)
println("  window UHF max |ΔFup| = ", res_uhf.maxdiff_up,
    ", max |ΔFdn| = ", res_uhf.maxdiff_dn)

psiup0 = orthonormal_cols(rng, N_r, 2)
res_sweep = compare_sweep_rhf(layout_r, Vblocks, H_r, psiup0;
    maxiter = 2, blocksize = 2, block_partition = layout_r, cutoff = 1e-8, verbose = false)
println("  sweep RHF energy proj = ", res_sweep.energy_proj,
    ", cached = ", res_sweep.energy_cached,
    ", |ΔE| = ", res_sweep.energy_diff)
