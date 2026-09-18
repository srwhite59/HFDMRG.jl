using HFDMRG
using LinearAlgebra
using Random

function fixture(N, nup, ndn)
    x = range(-1.0, 1.0; length = N)
    Hup = Matrix(Diagonal(x))
    Hdn = Matrix(Diagonal(reverse(x) .+ 0.01))
    V = [0.02 / (1 + abs(i - j)) for i = 1:N, j = 1:N]
    rng = MersenneTwister(91702)
    Cup = Matrix(qr(randn(rng, N, nup)).Q)[:, 1:nup]
    Cdn = Matrix(qr(randn(rng, N, ndn)).Q)[:, 1:ndn]
    Hup, Hdn, V, Cup, Cdn
end

function dense_uhf(Hup, Hdn, V, Cup, Cdn)
    Dup, Ddn = Cup * Cup', Cdn * Cdn'
    direct = V * (diag(Dup) + diag(Ddn))
    Fup = Hup + Diagonal(direct) - V .* Dup
    Fdn = Hdn + Diagonal(direct) - V .* Ddn
    FupC, FdnC = Fup * Cup, Fdn * Cdn
    Rup = FupC - Cup * (Cup' * FupC)
    Rdn = FdnC - Cdn * (Cdn' * FdnC)
    0.5dot(Dup, Fup + Hup) + 0.5dot(Ddn, Fdn + Hdn),
        max(norm(Rup), norm(Rdn))
end

function dense_rhf(H, V, C)
    D = C * C'
    d = diag(D)
    J = V * d
    F = H + 2Diagonal(J) - V .* D
    FC = F * C
    2dot(D, H) + 2dot(d, J) - sum(V .* D .* D),
        norm(FC - C * (C' * FC))
end

function sliced_fixture(slices, nup, ndn)
    rng = MersenneTwister(91703)
    dims = fill(2, slices)
    layout = HFDMRG.SliceLayout(dims)
    factors = [cat([Matrix(Symmetric(randn(rng, 2, 2))) for _ = 1:2]...;
        dims = 3) for _ = 1:slices]
    V = zeros(2, 2, 2, 2, slices, slices)
    for s = 1:slices, t = 1:slices
        V[:, :, :, :, s, t] = 0.002reshape(
            reshape(factors[s], 4, 2) * reshape(factors[t], 4, 2)', 2, 2, 2, 2)
    end
    N = 2slices
    A = randn(rng, N, N)
    Hup = Matrix(Symmetric(A))
    Hdn = Hup + Diagonal(range(-0.02, 0.03; length = N))
    Cup = Matrix(qr(randn(rng, N, nup)).Q)[:, 1:nup]
    Cdn = Matrix(qr(randn(rng, N, ndn)).Q)[:, 1:ndn]
    HFDMRG.SlicedBasisBackend(layout, V), Hup, Hdn, Cup, Cdn
end

function dense_sliced_rhf(H, backend, C)
    D, F = C * C', copy(H)
    HFDMRG.sliced_add_fock_r!(F, D, backend.V, backend.layout)
    FC = F * C
    dot(D, F + H), norm(FC - C * (C' * FC))
end

function dense_sliced_uhf(Hup, Hdn, backend, Cup, Cdn)
    Dup, Ddn = Cup * Cup', Cdn * Cdn'
    Fup, Fdn = copy(Hup), copy(Hdn)
    HFDMRG.sliced_add_fock_uhf!(Fup, Fdn, Dup, Ddn, backend.V, backend.layout)
    FupC, FdnC = Fup * Cup, Fdn * Cdn
    Rup = FupC - Cup * (Cup' * FupC)
    Rdn = FdnC - Cdn * (Cdn' * FdnC)
    0.5dot(Dup, Fup + Hup) + 0.5dot(Ddn, Fdn + Hdn),
        max(norm(Rup), norm(Rdn))
end

mode = isempty(ARGS) ? "matched" : ARGS[1]
N = length(ARGS) < 2 ? 600 : parse(Int, ARGS[2])
nup, ndn = 6, 5
println("mode=", mode, " parameter=", N, " nup=", nup, " ndn=", ndn,
    " julia=", VERSION, " blas_threads=", BLAS.get_num_threads())

if mode == "matched"
    Hup, Hdn, V, Cup, Cdn = fixture(N, nup, ndn)
    println("physical_dimension=", N, " permanent_input_bytes=",
        Base.summarysize((Hup, Hdn, V, Cup, Cdn)))
    calls = (("old_rhf", () -> dense_rhf(Hup, V, Cup)),
        ("new_rhf", () -> HFDMRG._history_physical(Hup, V, Cup)),
        ("old_uhf", () -> dense_uhf(Hup, Hdn, V, Cup, Cdn)),
        ("new_uhf", () -> HFDMRG._history_uhf_physical(Hup, Hdn, V, Cup, Cdn)))
    for (_, call) in calls
        call()
    end
    for (name, call) in calls
        GC.gc()
        sample = @timed call()
        println(name, " seconds=", sample.time, " bytes=", sample.bytes,
            " value=", sample.value)
    end
elseif mode == "matched_sliced"
    backend, Hsu, Hsd, Csu, Csd = sliced_fixture(N, nup, ndn)
    println("slices=", N, " physical_dimension=", size(Hsu, 1),
        " permanent_input_bytes=", Base.summarysize((backend, Hsu, Hsd, Csu, Csd)))
    calls = (("old_sliced_rhf", () -> dense_sliced_rhf(Hsu, backend, Csu)),
        ("new_sliced_rhf", () -> HFDMRG._history_physical(Hsu, backend, Csu)),
        ("old_sliced_uhf", () -> dense_sliced_uhf(Hsu, Hsd, backend, Csu, Csd)),
        ("new_sliced_uhf", () -> HFDMRG._history_uhf_physical(
            Hsu, Hsd, backend, Csu, Csd)))
    for (_, call) in calls
        call()
    end
    for (name, call) in calls
        GC.gc()
        sample = @timed call()
        println(name, " seconds=", sample.time, " bytes=", sample.bytes,
            " value=", sample.value)
    end
elseif mode == "old_uhf"
    Hup, Hdn, V, Cup, Cdn = fixture(N, nup, ndn)
    println("physical_dimension=", N, " permanent_input_bytes=",
        Base.summarysize((Hup, Hdn, V, Cup, Cdn)))
    fresh = @timed dense_uhf(Hup, Hdn, V, Cup, Cdn)
    println("seconds=", fresh.time, " bytes=", fresh.bytes, " result=", fresh.value)
    println("process_peak_rss_bytes=", Sys.maxrss())
elseif mode == "new_uhf"
    Hup, Hdn, V, Cup, Cdn = fixture(N, nup, ndn)
    println("physical_dimension=", N, " permanent_input_bytes=",
        Base.summarysize((Hup, Hdn, V, Cup, Cdn)))
    fresh = @timed HFDMRG._history_uhf_physical(Hup, Hdn, V, Cup, Cdn)
    println("seconds=", fresh.time, " bytes=", fresh.bytes,
        " result=", (fresh.value.energy, fresh.value.residual))
    println("process_peak_rss_bytes=", Sys.maxrss())
else
    error("mode must be matched, matched_sliced, old_uhf, or new_uhf")
end
