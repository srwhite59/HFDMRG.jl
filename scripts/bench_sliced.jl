using LinearAlgebra
using Random

using HFDMRG

const SEED, D, K, BATCHES, TARGET_SECONDS = 20260719, 9, 4, 3, 0.25

function orthonormal_cols(rng, n, m)
    Q = Matrix(qr(randn(rng, n, m)).Q)
    Q[:, 1:m]
end

function cache_plan(rng, ns)
    A = cos(0.2) * Matrix{Float64}(I, K, K); C = sin(0.2) * orthonormal_cols(rng, D, K)
    left, right = [orthonormal_cols(rng, D, K)], [orthonormal_cols(rng, D, K)]
    for _ = 2:(ns - 2)
        push!(left, vcat(left[end] * A, C)); push!(right, vcat(C, right[end] * A))
    end
    (; A, C, left, right)
end

slice_range(s) = ((s - 1) * D + 1):(s * D)

function init_edge(side, backend, plan, ns); ra = slice_range(side === :left ? 1 : ns); phi = side === :left ? plan.left[1] : plan.right[1]; raV = side === :left ? ((last(ra) + 1):(ns * D)) : (1:(first(ra) - 1)); HFDMRG.vee_init_block(side, ra, raV, phi, backend); end

function grow(side, state, backend, plan, ns, nblocks)
    phis = side === :left ? plan.left : plan.right
    for j = 2:nblocks
        cra = slice_range(side === :left ? j : ns - j + 1); phi = phis[j]
        raV = side === :left ? ((last(cra) + 1):(ns * D)) : (1:(first(cra) - 1))
        state = HFDMRG.vee_absorb_block(side, state, cra, plan.A, plan.C, phi, raV, backend)
    end
    state
end

function stage(f); f(); samples = [(GC.gc(); t = @timed f(); (t.time, t.bytes)) for _ = 1:BATCHES]; times = sort(first.(samples)); (; time = times[2], range = extrema(times), bytes = maximum(last.(samples))); end
function counted(f); ENV["HFDMRG_BENCH_TIMING"] = "1"; HFDMRG._bench_timing_reset!(); f(); counts = HFDMRG._bench_timing_snapshot(); ENV["HFDMRG_BENCH_TIMING"] = ""; counts; end
function run_calls(f, n); for _ = 1:n; f(); end; end
function batches(f); run_calls(f, 3); sample = @elapsed run_calls(f, 3); calls = clamp(ceil(Int, 3TARGET_SECONDS / max(sample, eps())), 3, 100_000); times = [(GC.gc(); (@elapsed run_calls(f, calls)) / calls) for _ = 1:BATCHES]; GC.gc(); bytes = @allocated run_calls(f, calls); sort!(times); (; time = times[2], range = extrema(times), calls, bytes); end
function fock_calls(win, rho, rhoup, rhodn); F, Fup, Fdn = zeros(size(rho)), zeros(size(rho)), zeros(size(rho)); (; F, Fup, Fdn, rhf = () -> (HFDMRG.vee_add_fock_r!(F, rho, win); nothing), uhf = () -> (HFDMRG.vee_add_fock!(Fup, Fdn, rhoup, rhodn, win); nothing)); end

relerr(A, B) = norm(A - B) / max(1.0, norm(A), norm(B))
route(c, keys) = Tuple(round(Int, c[key]) for key in keys)
showstage(s, x) = println(s, ": median=", x.time, " s range=", x.range, " allocation=", round(x.bytes / 2.0^20, digits = 3), " MiB")
showfock(s, x) = println(s, ": median=", x.time * 1e3, " ms/call range=", x.range .* 1e3, " calls/batch=", x.calls, " batches=3 allocation=", x.bytes, " bytes/batch")

function bench_size(ns)
    rng = MersenneTwister(SEED + ns); layout = HFDMRG.SliceLayout(fill(D, ns))
    V = randn(rng, D, D, D, D, ns, ns); plan, center = cache_plan(rng, ns), div(ns, 2) + 1
    cached, projection = HFDMRG.SlicedBasisBackendCached(layout, V),
        HFDMRG.SlicedBasisBackend(layout, V)
    initial = () -> begin L, R = init_edge(:left, cached, plan, ns), init_edge(:right, cached, plan, ns); (L, grow(:right, R, cached, plan, ns, ns - 2)); end
    L0, R0 = init_edge(:left, cached, plan, ns), init_edge(:right, cached, plan, ns)
    sweep = () -> (grow(:left, L0, cached, plan, ns, ns - 2),
        grow(:right, R0, cached, plan, ns, ns - 2))
    initial_result, sweep_result = stage(initial), stage(sweep); initial_counts = counted(initial)
    sweep_counts = counted(sweep)

    L = grow(:left, L0, cached, plan, ns, center - 1); R = grow(:right, R0, cached, plan, ns, ns - center)
    Cra = slice_range(center); Lra, Rra = 1:(first(Cra) - 1), (last(Cra) + 1):(ns * D)
    Lp = HFDMRG.vee_init_block(:left, Lra, first(Cra):(ns * D), L.phi, projection)
    Rp = HFDMRG.vee_init_block(:right, Rra, 1:last(Cra), R.phi, projection)
    cwindow, pwindow = () -> HFDMRG.vee_window(L, R, Cra, cached),
        () -> HFDMRG.vee_window(Lp, Rp, Cra, projection)
    cw_result, pw_result = stage(cwindow), stage(pwindow); window_counts = counted(cwindow)
    cw, pw = cwindow(), pwindow()

    n = 2K + D; symmetric() = (x = randn(rng, n, n); (x + x') / 2)
    rho, rhoup, rhodn = symmetric(), symmetric(), symmetric()
    cc, pc = fock_calls(cw, rho, rhoup, rhodn), fock_calls(pw, rho, rhoup, rhodn)
    cc.rhf(); pc.rhf(); cc.uhf(); pc.uhf()
    parity = (relerr(cc.F, pc.F), max(relerr(cc.Fup, pc.Fup), relerr(cc.Fdn, pc.Fdn)))
    cr, pr, cu, pu = batches(cc.rhf), batches(pc.rhf), batches(cc.uhf), batches(pc.uhf)

    ck = (:cache_init_calls, :cache_incremental_calls, :cache_fallback_calls)
    wk = (:window_local_calls, :window_projection_calls)
    ic, sc, wc = route(initial_counts, ck), route(sweep_counts, ck), route(window_counts, wk)
    println("\nS=", ns, " d=9 kL=kR=4 center_slices=1 seed=", SEED + ns)
    showstage("  cache initial", initial_result); println("    route=", ic); showstage("  cache sweep", sweep_result)
    println("    route=", sc)
    showstage("  cached window", cw_result); showstage("  projection window", pw_result)
    println("    window route(local,projection)=", wc, " parity(RHF,UHF)=", parity)
    for (label, result) in (("cached RHF", cr), ("projection RHF", pr),
        ("cached UHF", cu), ("projection UHF", pu))
        showfock("  " * label, result)
    end

    ic == (2, ns - 3, 0) && sc == (0, 2(ns - 3), 0) || error("cache route failed")
    wc == (1, 0) || error("aligned window used fallback")
    max(parity...) <= 1e-12 || error("Fock parity failed")
    cr.bytes == 0 && cu.bytes == 0 || error("cached steady Fock allocated")
    ns == 52 && (initial_result.time <= 16 && sweep_result.time <= 32 &&
        initial_result.bytes + sweep_result.bytes <= 512 * 2^20 || error("cache-chain gate failed"))
    (; rhf = (cached = cr.time, projection = pr.time),
        uhf = (cached = cu.time, projection = pu.time))
end

function main()
    ENV["HFDMRG_BENCH_TIMING"] = ""
    println("host=", gethostname(), " commit=", readchomp(`git rev-parse --short HEAD`), " julia=", VERSION,
        " machine=", Sys.MACHINE, " julia_threads=", Threads.nthreads(), " blas_threads=", BLAS.get_num_threads(), " batches=3 statistic=median(per-call batch means)")
    r8 = bench_size(8); GC.gc(); r52 = bench_size(52)
    for spin in (:rhf, :uhf)
        a, b = getproperty(r8, spin), getproperty(r52, spin)
        speedup, scaling = b.projection / b.cached, b.cached / a.cached
        println(spin, " S52 speedup=", speedup, " cached S52/S8=", scaling)
        speedup >= 20 && scaling <= 1.5 || error("$spin performance gate failed")
    end
end

main()
