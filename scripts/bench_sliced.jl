# Benchmark snapshot (make bench):
# N = 120, ns = 30, nj = 4, ml = 10, lc = 40, mr = 10, nrep = 200
# H norm = 85.22616536592902
# == Projection backend ==
#   split slices: n/a
#   RHF: 1.233 ms/call, 0.423 MB/call
#   UHF: 1.459 ms/call, 0.977 MB/call
# == Cached backend ==
#   split slices: none
#   RHF: 0.269 ms/call, 0.051 MB/call
#   RHF timing breakdown (ms/call):
#   direct_LL_RR: 0.011
#   direct_CL_CR: 0.099
#   direct_LR: 0.011
#   exchange: 0.111
#   UHF: 0.435 ms/call, 0.109 MB/call
#   UHF timing breakdown (ms/call):
#   direct_LL_RR: 0.012
#   direct_CL_CR: 0.106
#   direct_LR: 0.012
#   exchange: 0.224
ENV["HFDMRG_BENCH_TIMING"] = "1"

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

function split_label(win)
    if hasproperty(win, :split_slices)
        isempty(win.split_slices) ? "none" : join(win.split_slices, ",")
    else
        "n/a"
    end
end

function bench_rhf!(win, rho; nrep)
    F = zeros(size(rho))
    HFDMRG.vee_add_fock_r!(F, rho, win)

    if win isa HFDMRG.SlicedBasisCachedWindow
        HFDMRG._bench_timing_reset!()
    end

    t0 = time_ns()
    for _ = 1:nrep
        fill!(F, 0.0)
        HFDMRG.vee_add_fock_r!(F, rho, win)
    end
    dt = (time_ns() - t0) * 1e-9

    timing = win isa HFDMRG.SlicedBasisCachedWindow ? HFDMRG._bench_timing_snapshot() : Dict{Symbol, Float64}()

    env_flag = get(ENV, "HFDMRG_BENCH_TIMING", "")
    if win isa HFDMRG.SlicedBasisCachedWindow
        ENV["HFDMRG_BENCH_TIMING"] = ""
    end
    alloc = @allocated begin
        for _ = 1:nrep
            fill!(F, 0.0)
            HFDMRG.vee_add_fock_r!(F, rho, win)
        end
    end
    if win isa HFDMRG.SlicedBasisCachedWindow
        ENV["HFDMRG_BENCH_TIMING"] = env_flag
    end
    (time_per_call = dt / nrep, alloc_per_call = alloc / nrep, timing = timing)
end

function bench_uhf!(win, rhoup, rhodn; nrep)
    Fup = zeros(size(rhoup))
    Fdn = zeros(size(rhoup))
    HFDMRG.vee_add_fock!(Fup, Fdn, rhoup, rhodn, win)

    if win isa HFDMRG.SlicedBasisCachedWindow
        HFDMRG._bench_timing_reset!()
    end

    t0 = time_ns()
    for _ = 1:nrep
        fill!(Fup, 0.0)
        fill!(Fdn, 0.0)
        HFDMRG.vee_add_fock!(Fup, Fdn, rhoup, rhodn, win)
    end
    dt = (time_ns() - t0) * 1e-9

    timing = win isa HFDMRG.SlicedBasisCachedWindow ? HFDMRG._bench_timing_snapshot() : Dict{Symbol, Float64}()

    env_flag = get(ENV, "HFDMRG_BENCH_TIMING", "")
    if win isa HFDMRG.SlicedBasisCachedWindow
        ENV["HFDMRG_BENCH_TIMING"] = ""
    end
    alloc = @allocated begin
        for _ = 1:nrep
            fill!(Fup, 0.0)
            fill!(Fdn, 0.0)
            HFDMRG.vee_add_fock!(Fup, Fdn, rhoup, rhodn, win)
        end
    end
    if win isa HFDMRG.SlicedBasisCachedWindow
        ENV["HFDMRG_BENCH_TIMING"] = env_flag
    end
    (time_per_call = dt / nrep, alloc_per_call = alloc / nrep, timing = timing)
end

function print_timing(label, res)
    println(label, ": ", round(res.time_per_call * 1e3, digits = 3), " ms/call, ",
        round(res.alloc_per_call / 1e6, digits = 3), " MB/call")
end

function print_cached_breakdown(label, timing, nrep)
    isempty(timing) && return
    println(label, " timing breakdown (ms/call):")
    println("  direct_LL_RR: ", round(timing[:direct_LL_RR] / nrep * 1e3, digits = 3))
    println("  direct_CL_CR: ", round(timing[:direct_CL_CR] / nrep * 1e3, digits = 3))
    println("  direct_LR: ", round(timing[:direct_LR] / nrep * 1e3, digits = 3))
    println("  exchange: ", round(timing[:exchange] / nrep * 1e3, digits = 3))
end

function bench_backend(name, win, rho, rhoup, rhodn; nrep)
    println("== ", name, " ==")
    println("  split slices: ", split_label(win))

    rhf = bench_rhf!(win, rho; nrep = nrep)
    print_timing("  RHF", rhf)
    print_cached_breakdown("  RHF", rhf.timing, nrep)

    uhf = bench_uhf!(win, rhoup, rhodn; nrep = nrep)
    print_timing("  UHF", uhf)
    print_cached_breakdown("  UHF", uhf.timing, nrep)
end

function main()
    rng = MersenneTwister(1234)
    ns = 30
    nj = 4
    N = ns * nj
    nrep = 200

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

    win_proj = build_window(backend_proj, Lra, Cra, Rra, Lphi, Rphi, N)
    win_cached = build_window(backend_cached, Lra, Cra, Rra, Lphi, Rphi, N)

    println("N = ", N, ", ns = ", ns, ", nj = ", nj,
        ", ml = ", ml, ", lc = ", length(Cra), ", mr = ", mr,
        ", nrep = ", nrep)
    println("H norm = ", norm(H))

    bench_backend("Projection backend", win_proj, rho, rhoup, rhodn; nrep = nrep)
    bench_backend("Cached backend", win_cached, rho, rhoup, rhodn; nrep = nrep)
end

main()
