using LinearAlgebra
using Random

using HFDMRG

function orthonormal_cols(rng, n, m)
    Q = Matrix(qr(randn(rng, n, m)).Q)
    Q[:, 1:m]
end

function build_window(backend, Lra, Cra, Rra, Lphi, Rphi, N)
    Lvee = HFDMRG.vee_init_block(:left, Lra, (Lra[end] + 1):N, Lphi, backend)
    Rvee = HFDMRG.vee_init_block(:right, Rra, 1:(Rra[1] - 1), Rphi, backend)
    HFDMRG.vee_window(Lvee, Rvee, Cra, backend)
end

function bench_backend(name, backend, Lra, Cra, Rra, Lphi, Rphi, N, rho, rhoup, rhodn)
    superdim = size(rho, 1)
    println("== ", name, " ==")

    # Warmup
    win = build_window(backend, Lra, Cra, Rra, Lphi, Rphi, N)
    F = zeros(superdim, superdim)
    HFDMRG.vee_add_fock_r!(F, rho, win)
    Fup = zeros(superdim, superdim)
    Fdn = zeros(superdim, superdim)
    HFDMRG.vee_add_fock!(Fup, Fdn, rhoup, rhodn, win)

    println("RHF (window + Fock)")
    t = @elapsed begin
        win = build_window(backend, Lra, Cra, Rra, Lphi, Rphi, N)
        F = zeros(superdim, superdim)
        HFDMRG.vee_add_fock_r!(F, rho, win)
    end
    b = @allocated begin
        win = build_window(backend, Lra, Cra, Rra, Lphi, Rphi, N)
        F = zeros(superdim, superdim)
        HFDMRG.vee_add_fock_r!(F, rho, win)
    end
    println("  elapsed: ", round(t, digits = 3), " s")
    println("  allocated: ", round(b / 1e6, digits = 2), " MB")

    println("UHF (window + Fock)")
    t = @elapsed begin
        win = build_window(backend, Lra, Cra, Rra, Lphi, Rphi, N)
        Fup = zeros(superdim, superdim)
        Fdn = zeros(superdim, superdim)
        HFDMRG.vee_add_fock!(Fup, Fdn, rhoup, rhodn, win)
    end
    b = @allocated begin
        win = build_window(backend, Lra, Cra, Rra, Lphi, Rphi, N)
        Fup = zeros(superdim, superdim)
        Fdn = zeros(superdim, superdim)
        HFDMRG.vee_add_fock!(Fup, Fdn, rhoup, rhodn, win)
    end
    println("  elapsed: ", round(t, digits = 3), " s")
    println("  allocated: ", round(b / 1e6, digits = 2), " MB")
end

function main()
    rng = MersenneTwister(1234)
    ns = 30
    nj = 4
    N = ns * nj

    A = randn(rng, N, N)
    H = (A + A') / 2
    V6 = randn(rng, nj, nj, nj, nj, ns, ns)

    Lra = 1:40
    Cra = 41:80
    Rra = 81:120
    ml = 10
    mr = 10
    Lphi = orthonormal_cols(rng, length(Lra), ml)
    Rphi = orthonormal_cols(rng, length(Rra), mr)
    superdim = ml + length(Cra) + mr

    rho = randn(rng, superdim, superdim)
    rho = (rho + rho') / 2
    rhoup = randn(rng, superdim, superdim)
    rhoup = (rhoup + rhoup') / 2
    rhodn = randn(rng, superdim, superdim)
    rhodn = (rhodn + rhodn') / 2

    layout = HFDMRG.SliceLayout(fill(nj, ns))
    backend_proj = HFDMRG.SlicedBasisBackend(layout, V6)
    backend_cached = HFDMRG.SlicedBasisBackendCached(layout, V6)

    println("N = ", N, ", ml = ", ml, ", lc = ", length(Cra), ", mr = ", mr)
    println("H norm = ", norm(H))
    bench_backend("Projection backend", backend_proj, Lra, Cra, Rra, Lphi, Rphi, N, rho,
        rhoup, rhodn)
    bench_backend("Cached backend", backend_cached, Lra, Cra, Rra, Lphi, Rphi, N, rho,
        rhoup, rhodn)
end

main()
