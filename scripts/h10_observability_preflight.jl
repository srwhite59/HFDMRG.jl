using HFDMRG
using DelimitedFiles
using LinearAlgebra
using Random
using SHA
using Sockets
const CONTRACT_SHA = "d42704cc4312b92b7b6db8c96b871d251d1837b14008374370344cb9a591b445"
const MANIFEST_SHA = "f635299798b2f08cb848604465b133efcd86c026e67055f11500801c32b91cfd"
const V6_SHA = "1ee248862c8f0193f1462452b05fdb9bcdefaec39aafb7ab9e056c4279f9ef2b"
const INPUT_SHA = Dict(
    "data/H.bin" => "66ae521070e41c91543d353ca080f4ccc094570c15d4ab6b2cb993fe2bd16957",
    "data/C.bin" => "396da26cbaa1d26c2c24f519bbdf25334567ed9619a2b54ca4c6de12817c6c51",
    "data/vg.bin" => "4fa57a577153c9d52def4e41789473ef70ac27af378d801295be657f4bb025ea",
    "data/atom_local_af_s1_alpha.tsv" =>
        "21143810922ce920406d4d94c1261f9dc942a69fe115cb330f9ffacbb6177643",
    "data/atom_local_af_s1_beta.tsv" =>
        "d4b0b30d2b69e319a73fdb011c94e9dc5d9cff24523531dadc4267ccf809500a",
    "data/slice_permutations.tsv" =>
        "e56a3df2659568c2476fce55f43d990d8c3bbf0f417b5b296ed83da06fc9ecf5",
    "tools/read_compact_parent.jl" =>
        "d3c646567f23d20d81d3740531c7bb3aa8a4b45412d51beb71d8ffa711f2cc82",
    "tools/materialize_v6_from_compact_parent.jl" =>
        "855172629458beec71c27d651e4fba1560baa3ec07608fcb543f6432ff5f6b2a",
)
const NS, D, N, NOCC = 401, 4, 1604, 5
const SELECTED = Set((1, 200, 399, 796))
const GATE = 2e-11
gate(ok, message) = ok || error(message)
sha256_file(path) = bytes2hex(open(SHA.sha256, path))
array_sha(A) = bytes2hex(SHA.sha256(reinterpret(UInt8, vec(A))))
encode(x) = x isa UnitRange ? "$(first(x)):$(last(x))" :
    x isa AbstractVector ? join(x, ',') : replace(string(x), '\t' => ' ')

function arguments()
    defaults = Dict("mode" => "smoke", "order" => "", "archive" => "",
        "contract" => "", "v6" => "", "outdir" => "", "trials" => "3")
    for arg in ARGS
        key, value = split(lstrip(arg, '-'), '='; limit = 2)
        gate(haskey(defaults, key), "unknown option --$key")
        defaults[key] = value
    end
    mode = Symbol(defaults["mode"])
    gate(mode in (:smoke, :overhead, :check, :h10), "mode must be smoke, overhead, check, or h10")
    gate(!isempty(defaults["outdir"]), "--outdir=PATH is required")
    order = isempty(defaults["order"]) ? nothing : Symbol(defaults["order"])
    if mode in (:check, :h10)
        gate(order in (:natural, :odd_even), "--order=natural|odd_even is required")
        gate(all(!isempty(defaults[k]) for k in ("archive", "contract", "v6")),
            "--archive, --contract, and --v6 are required")
    end
    trials = parse(Int, defaults["trials"])
    gate(trials > 0 && (mode != :overhead || trials >= 3), "overhead mode requires at least 3 trials")
    (; mode, order, archive = abspath(defaults["archive"]),
        contract = abspath(defaults["contract"]), v6 = abspath(defaults["v6"]),
        outdir = abspath(defaults["outdir"]), trials)
end
function write_tsv(path, header, rows)
    tmp = path * ".tmp.$(getpid())"
    open(tmp, "w") do io
        println(io, join(header, '\t'))
        for row in rows
            println(io, join((encode(x) for x in row), '\t'))
        end
    end
    mv(tmp, path; force = true)
end
function read_array(path, dims)
    gate(filesize(path) == 8prod(dims), "unexpected byte count for $path")
    A = Array{Float64}(undef, dims)
    open(path, "r") do io
        read!(io, A)
        gate(eof(io), "trailing bytes in $path")
    end
    A
end
function read_start(path)
    table = readdlm(path, '\t', Float64; skipstart = 1)
    gate(size(table) == (N, NOCC + 1) && Int.(table[:, 1]) == 1:N,
        "invalid frozen start $path")
    Matrix(@view table[:, 2:end])
end
function read_orders(path)
    orders = Dict(:natural => Int[], :odd_even => Int[])
    for line in Iterators.drop(eachline(path), 1)
        fields = split(line, '\t')
        push!(orders[Symbol(fields[3])], parse(Int, fields[2]))
    end
    gate(orders[:natural] == 1:NS, "natural permutation mismatch")
    gate(orders[:odd_even] == vcat(collect(1:2:399), collect(2:2:400), 401),
        "odd-even permutation mismatch")
    orders
end
orbital_order(slices) = reduce(vcat, ((s - 1) * D .+ (1:D) for s in slices))
density_energy(Du, Dd, Fu, Fd, H) =
    0.5 * (dot(transpose(Du), Fu) + dot(transpose(Du), H) +
           dot(transpose(Dd), Fd) + dot(transpose(Dd), H))
# Independent physical contraction. It calls no HFDMRG backend routine.
function physical_fock(H, V6, Du, Dd)
    Fu, Fd = copy(H), copy(H)
    @inbounds for n = 1:NS, m = 1:NS, a = 1:D, b = 1:D, c = 1:D, d = 1:D
        p, r = (n - 1) * D + a, (n - 1) * D + b
        q, s = (m - 1) * D + c, (m - 1) * D + d
        v = V6[a, b, c, d, n, m]
        direct = (Du[q, s] + Dd[q, s]) * v
        Fu[p, r] += direct
        Fd[p, r] += direct
        Fu[p, s] -= Du[r, q] * v
        Fd[p, s] -= Dd[r, q] * v
    end
    Fu, Fd, density_energy(Du, Dd, Fu, Fd, H)
end
function verify_inputs(opts)
    gate(sha256_file(opts.contract) == CONTRACT_SHA, "governing contract hash mismatch")
    gate(sha256_file(joinpath(opts.archive, "SHA256SUMS")) == MANIFEST_SHA,
        "archive manifest hash mismatch")
    for (relative, expected) in INPUT_SHA
        gate(sha256_file(joinpath(opts.archive, relative)) == expected,
            "archive hash mismatch: $relative")
    end
    gate(filesize(opts.v6) == 329_320_448 && sha256_file(opts.v6) == V6_SHA,
        "regenerated V6 hash/size mismatch")
    H = read_array(joinpath(opts.archive, "data", "H.bin"), (N, N))
    Vnatural = read_array(opts.v6, (D, D, D, D, NS, NS))
    Cup = read_start(joinpath(opts.archive, "data", "atom_local_af_s1_alpha.tsv"))
    Cdn = read_start(joinpath(opts.archive, "data", "atom_local_af_s1_beta.tsv"))
    orders = read_orders(joinpath(opts.archive, "data", "slice_permutations.tsv"))
    gate(all(isfinite, H) && all(isfinite, Vnatural) &&
         norm(H - H') <= 512eps(Float64) * max(1, norm(H)), "nonfinite or asymmetric input")
    gate(max(norm(Cup' * Cup - I), norm(Cdn' * Cdn - I)) <= 2e-12,
        "frozen start is not orthonormal")
    slices = orders[opts.order]
    p = orbital_order(slices)
    gate(sort(p) == 1:N && invperm(p)[p] == 1:N, "orbital permutation round trip failed")
    if opts.order === :natural
        return (; H, V6 = Vnatural, Cup, Cdn, slices,
            covariance_energy = 0.0, covariance_fock = 0.0, ordered_v6_sha = V6_SHA)
    end
    Dun, Ddn = Cup * Cup', Cdn * Cdn'
    Fun, Fdn, En = physical_fock(H, Vnatural, Dun, Ddn)
    Hp, Cp, Dp = H[p, p], Cup[p, :], Cdn[p, :]
    Vp = Vnatural[:, :, :, :, slices, slices]
    ip, is = invperm(p), invperm(slices)
    gate(Hp[ip, ip] == H && Cp[ip, :] == Cup && Dp[ip, :] == Cdn,
        "ordered H/start round trip failed")
    gate(all(Vp[a, b, c, d, is[n], is[m]] == Vnatural[a, b, c, d, n, m]
        for n = 1:NS, m = 1:NS, a = 1:D, b = 1:D, c = 1:D, d = 1:D),
        "ordered V6 round trip failed")
    ordered_v6_sha = array_sha(Vp)
    Vnatural = nothing
    GC.gc()
    Dup, Ddp = Cp * Cp', Dp * Dp'
    Fup, Fdp, Ep = physical_fock(Hp, Vp, Dup, Ddp)
    efail = abs(Ep - En)
    ffail = max(norm(Fup - Fun[p, p]), norm(Fdp - Fdn[p, p]))
    gate(efail <= GATE && ffail <= GATE, "ordering covariance gate failed")
    (; H = Hp, V6 = Vp, Cup = Cp, Cdn = Dp, slices,
        covariance_energy = efail, covariance_fock = ffail, ordered_v6_sha)
end
mutable struct Meter{B}
    backend::B
    layout::HFDMRG.SliceLayout
    source::Vector{Int}
    selected::Set{Int}
    phase::Symbol
    times::Dict{Tuple{Symbol, Symbol}, Float64}
    counts::Dict{Tuple{Symbol, Symbol}, Int}
    live::Dict{Tuple{Symbol, Int}, Any}
    window_index::Int
    window_rows::Vector{Any}
    captures::Vector{Any}
    route_initial::Dict{Symbol, Float64}
    phase_started::UInt64
    walls::Dict{Symbol, Float64}
    observer_seconds::Float64
    boundary_payload::Int
    boundary_summary::Int
end
mutable struct MeterWindow{W, M}
    inner::W
    meter::M
    index::Int
    B::Union{Nothing, Matrix{Float64}}
    call::Int
end
function Meter(backend, source; selected = Set{Int}())
    HFDMRG._bench_timing_reset!()
    Meter(backend, backend.layout, source, selected, :initialization,
        Dict{Tuple{Symbol, Symbol}, Float64}(), Dict{Tuple{Symbol, Symbol}, Int}(),
        Dict{Tuple{Symbol, Int}, Any}(), 0, Any[], Any[], Dict{Symbol, Float64}(),
        time_ns(), Dict{Symbol, Float64}(), 0.0, 0, 0)
end
function counted(m, operation, f)
    key = (m.phase, operation)
    m.counts[key] = get(m.counts, key, 0) + 1
    t = time_ns()
    result = f()
    m.times[key] = get(m.times, key, 0.0) + (time_ns() - t) / 1e9
    result
end
counted(f, m::Meter, operation) = counted(m, operation, f)
function state_slot(m, side, ra)
    side === :left && return searchsortedfirst(m.layout.offs, last(ra)) - 1
    searchsortedlast(m.layout.offs, first(ra) - 1)
end
function track!(m, side, state)
    slot = state_slot(m, side, state.ra)
    delete!(m.live, (side === :left ? :right : :left, slot))
    m.live[(side, slot)] = state
    state
end
function HFDMRG.vee_init_block(side, ra, raV, phi, m::Meter)
    track!(m, side, counted(m, :init) do
        HFDMRG.vee_init_block(side, ra, raV, phi, m.backend)
    end)
end
function HFDMRG.vee_absorb_block(side, old, cra, Po, Pc, phi, raV, m::Meter)
    track!(m, side, counted(m, :absorb) do
        HFDMRG.vee_absorb_block(side, old, cra, Po, Pc, phi, raV, m.backend)
    end)
end
function independent_B(L, R, Cra, Nfull)
    ml, mr = size(L.phi, 2), size(R.phi, 2)
    B = zeros(Float64, Nfull, ml + length(Cra) + mr)
    B[L.ra, 1:ml] = L.phi
    for (column, row) in enumerate(Cra)
        B[row, ml + column] = 1
    end
    B[R.ra, ml + length(Cra) + 1:end] = R.phi
    B
end
function source_label(m, ra)
    d = m.layout.dims[1]
    slices = ((first(ra) - 1) ÷ d + 1):((last(ra) - 1) ÷ d + 1)
    join(m.source[slices], ',')
end
function HFDMRG.vee_window(L, R, Cra, m::Meter)
    if m.phase === :initialization
        m.walls[:initialization] = (time_ns() - m.phase_started) / 1e9
        m.phase_started = time_ns()
        m.route_initial = HFDMRG._bench_timing_snapshot()
        m.phase = :sweep
    end
    m.window_index += 1
    index = m.window_index
    ltr = length(m.layout.dims) - 2
    direction = index <= ltr ? :ltr : :rtl
    ordinal = direction === :ltr ? index : index - ltr
    B = index in m.selected ? independent_B(L, R, Cra, m.layout.offs[end]) : nothing
    inner = counted(m, :window) do
        HFDMRG.vee_window(L, R, Cra, m.backend)
    end
    push!(m.window_rows, (index, direction, ordinal, first(L.ra), last(L.ra),
        first(Cra), last(Cra), first(R.ra), last(R.ra), size(L.phi, 2),
        size(R.phi, 2), source_label(m, L.ra), source_label(m, Cra),
        source_label(m, R.ra)))
    MeterWindow(inner, m, index, B, 0)
end
function HFDMRG.vee_add_fock!(Fu, Fd, Du, Dd, w::MeterWindow)
    call = w.call
    w.call += 1
    selected = w.B !== nothing
    Hu = selected ? copy(Fu) : nothing
    Hd = selected ? copy(Fd) : nothing
    Dcu = selected ? copy(Du) : nothing
    Dcd = selected ? copy(Dd) : nothing
    counted(w.meter, :fock) do
        HFDMRG.vee_add_fock!(Fu, Fd, Du, Dd, w.inner)
    end
    selected && push!(w.meter.captures, (window = w.index, call, B = w.B,
        Du = Dcu, Dd = Dcd, Hu, Hd, Fu = copy(Fu), Fd = copy(Fd)))
    nothing
end
function HFDMRG.vee_add_fock_r!(F, Dmat, w::MeterWindow)
    counted(w.meter, :fock) do
        HFDMRG.vee_add_fock_r!(F, Dmat, w.inner)
    end
end
function payload_bytes!(x, seen)
    x isa Array || return 0
    haskey(seen, x) && return 0
    seen[x] = nothing
    bytes = sizeof(x)
    isbitstype(eltype(x)) || (bytes += sum(payload_bytes!(v, seen) for v in x))
    bytes
end
function cache_memory(m)
    states = collect(values(m.live))
    seen = IdDict{Any, Nothing}()
    payload = sum(payload_bytes!(getfield(state, field), seen)
        for state in states for field in fieldnames(typeof(state)))
    payload, Base.summarysize(states)
end
function oracle_rows(packet, captures, windows)
    rows = Any[]
    for (i, c) in enumerate(captures)
        Du = c.B * c.Du * c.B'
        Dd = c.B * c.Dd * c.B'
        Fu, Fd, Eoracle = physical_fock(packet.H, packet.V6, Du, Dd)
        Hlocal = c.B' * packet.H * c.B
        berr = norm(c.B' * c.B - I)
        herr = max(maximum(abs, c.Hu - Hlocal), maximum(abs, c.Hd - Hlocal))
        ferru = norm(c.Fu - c.B' * Fu * c.B)
        ferrd = norm(c.Fd - c.B' * Fd * c.B)
        asym = max(norm(Fu - Fu') / max(1, norm(Fu)),
            norm(Fd - Fd') / max(1, norm(Fd)))
        asym_gate = 512eps(Float64)
        Elive = c.call == 0 ? NaN : density_energy(c.Du, c.Dd, c.Fu, c.Fd, c.Hu)
        Eprod = c.call == 0 ? NaN : windows[c.window].local_energies[c.call]
        eerr = c.call == 0 ? NaN : max(abs(Elive - Eoracle), abs(Eprod - Eoracle))
        push!(rows, (c.window, c.call, size(c.B, 2), berr, herr, ferru, ferrd,
            asym, asym_gate, Elive, Eprod, Eoracle, eerr))
        gate(berr <= GATE && herr <= GATE && max(ferru, ferrd) <= GATE &&
             asym <= asym_gate && (c.call == 0 ||
             (eerr <= GATE && abs(Elive - Eprod) <= GATE)),
             "actual-window oracle gate failed")
        println("oracle ", i, "/", length(captures), " window=", c.window, " call=", c.call)
        flush(stdout)
    end
    rows
end
function write_diagnostics(outdir, diagnostics)
    scalar_rows = Any[]
    for name in propertynames(diagnostics)
        value = getproperty(diagnostics, name)
        if value isa AbstractVector && (isempty(value) || !isbitstype(eltype(value)))
            isempty(value) && continue
            header = propertynames(first(value))
            rows = [Tuple(getproperty(record, key) for key in header) for record in value]
            write_tsv(joinpath(outdir, "diagnostics_$(name).tsv"), header, rows)
        else
            push!(scalar_rows, (name, encode(value)))
        end
    end
    write_tsv(joinpath(outdir, "diagnostics_summary.tsv"), ("quantity", "value"), scalar_rows)
end
function small_fixture()
    rng = MersenneTwister(42)
    ns, d, n = 7, 2, 14
    H = randn(rng, n, n)
    H = 0.5 * (H + H')
    V6 = zeros(d, d, d, d, ns, ns)
    for i = 1:ns, j = 1:ns, a = 1:d, c = 1:d
        V6[a, a, c, c, i, j] = 0.2 / (1 + abs(i - j))
    end
    Q = Matrix(qr(randn(rng, n, 4)).Q)
    (; H, V6, Cup = Q[:, 1:2], Cdn = Q[:, 3:4],
        layout = HFDMRG.SliceLayout(fill(d, ns)))
end
make_diagnostics() = HFDMRG._RunDiagnostics(:detailed)
function solve_private(H, backend, Cup, Cdn, layout, diagnostics; observer = nothing)
    kwargs = diagnostics === nothing ? NamedTuple() : (; _diagnostics = diagnostics)
    HFDMRG.solve_hfdmrg_core(H, backend, Cup, Cdn; restricted = false,
        block_partition = layout, nblockcenter = 1, maxiter = 1, cutoff = 0.0,
        scf_cutoff = 1e-12, observer, verbose = false, kwargs...)
end
function run_small(outdir, mode, trials)
    packet = small_fixture()
    function once(kind)
        raw = HFDMRG.SlicedBasisBackendCached(packet.layout, packet.V6)
        diagnostics = kind === :raw_off || kind === :meter_off ? nothing : make_diagnostics()
        backend = startswith(string(kind), "meter") ?
            Meter(raw, collect(1:length(packet.layout.dims))) : raw
        observer = backend isa Meter ? (_ -> (gate(length(backend.live) ==
            length(packet.layout.dims), "small live-cache count mismatch"); cache_memory(backend); false)) : nothing
        timed = @timed solve_private(packet.H, backend, packet.Cup, packet.Cdn,
            packet.layout, diagnostics; observer)
        timed.value, timed.time, timed.bytes, diagnostics
    end
    warm = once(:raw_off)[1]
    rows = Any[]
    for trial = 1:(mode === :smoke ? 1 : trials), kind in (:raw_off, :meter_off, :meter_on)
        result, wall, bytes, diagnostics = once(kind)
        gate(result == warm, "small-fixture instrumentation changed arithmetic")
        push!(rows, (trial, kind, wall, bytes, result[3],
            diagnostics === nothing ? 0 : length(getproperty(diagnostics, :windows))))
    end
    write_tsv(joinpath(outdir, "small_runs.tsv"),
        ("trial", "route", "wall_seconds", "allocated_bytes", "energy", "diagnostic_windows"), rows)
end
function run_h10(opts, packet)
    orderdir = joinpath(opts.outdir, string(opts.order))
    mkpath(orderdir)
    layout = HFDMRG.SliceLayout(fill(D, NS))
    raw = HFDMRG.SlicedBasisBackendCached(layout, packet.V6)
    meter = Meter(raw, packet.slices; selected = SELECTED)
    diagnostics = make_diagnostics()
    result_ref = Ref{Any}()
    observer = info -> begin
        t = time_ns()
        meter.walls[:sweep] = (t - meter.phase_started) / 1e9
        gate(info.sweep == 1 && !info.converged, "preflight did not have exact forced length")
        gate(length(meter.live) == NS, "live cache state count mismatch")
        meter.boundary_payload, meter.boundary_summary = cache_memory(meter)
        result_ref[] = (info.energy, info.psiup, info.psidn)
        meter.observer_seconds += (time_ns() - t) / 1e9
        false
    end
    total = @timed solve_private(packet.H, meter, packet.Cup, packet.Cdn,
        layout, diagnostics; observer)
    gate(meter.window_index == 796, "first-sweep window count mismatch")
    route_final = HFDMRG._bench_timing_snapshot()
    route_sweep = Dict(key => route_final[key] - get(meter.route_initial, key, 0.0)
        for key in keys(route_final))
    gate(get(route_final, :cache_init_calls, -1) == 2 &&
         get(meter.route_initial, :cache_incremental_calls, -1) == 398 &&
         get(route_sweep, :cache_incremental_calls, -1) == 796 &&
         get(route_sweep, :window_local_calls, -1) == 796 &&
         get(route_final, :cache_fallback_calls, -1) == 0 &&
         get(route_final, :window_projection_calls, -1) == 0, "cache route gate failed")
    expected = Dict(1 => (:ltr, 1, 1, 4, 5, 8, 9, N),
        200 => (:ltr, 200, 1, 800, 801, 804, 805, N),
        399 => (:ltr, 399, 1, 1596, 1597, 1600, 1601, N),
        796 => (:rtl, 397, 1, 8, 9, 12, 13, N))
    gate(all(meter.window_rows[i][2:9] == row for (i, row) in expected),
        "selected window geometry mismatch")
    gate(Set(c.window for c in meter.captures) == SELECTED &&
         length(meter.captures) <= 20, "selected-window capture count mismatch")
    calls = Dict(w => count(c -> c.window == w, meter.captures) for w in SELECTED)
    gate(all(2 <= count <= 5 for count in values(calls)), "selected local-call count mismatch")
    gate(diagnostics.windows[end].final_energy == result_ref[][1],
        "last hidden/reported energy mismatch")
    oracle_timed = @timed oracle_rows(packet, meter.captures, diagnostics.windows)
    oracle = oracle_timed.value
    write_tsv(joinpath(orderdir, "window_routes.tsv"),
        ("global_window", "direction", "ordinal", "L_first", "L_last", "C_first",
         "C_last", "R_first", "R_last", "left_rank", "right_rank",
         "L_source_slices", "C_source_slices", "R_source_slices"), meter.window_rows)
    write_tsv(joinpath(orderdir, "window_oracle.tsv"),
        ("global_window", "fock_call", "window_rank", "basis_orthogonality",
         "one_body_max_error",
         "alpha_fock_frobenius_error", "beta_fock_frobenius_error",
         "fock_antisymmetry_ratio", "antisymmetry_ratio_gate", "live_energy",
         "production_energy", "oracle_energy", "energy_error"), oracle)
    write_tsv(joinpath(orderdir, "oracle_work.tsv"),
        ("wall_seconds", "allocated_bytes", "gc_seconds"),
        [(oracle_timed.time, oracle_timed.bytes, oracle_timed.gctime)])
    phase_rows = vec([(phase, op, get(meter.times, (phase, op), 0.0),
        get(meter.counts, (phase, op), 0)) for phase in (:initialization, :sweep),
        op in (:init, :absorb, :window, :fock)])
    append!(phase_rows, [(phase, :total, meter.walls[phase], 1)
        for phase in (:initialization, :sweep)])
    write_tsv(joinpath(orderdir, "backend_phases.tsv"),
        ("phase", "operation", "seconds", "calls"), phase_rows)
    write_tsv(joinpath(orderdir, "routes.tsv"), ("phase", "route", "calls"),
        vcat([(:initialization, k, v) for (k, v) in sort!(collect(meter.route_initial))],
             [(:sweep, k, v) for (k, v) in sort!(collect(route_sweep))]))
    energy, Cup, Cdn = result_ref[]
    write_tsv(joinpath(orderdir, "preflight_result.tsv"),
        ("label", "reported_energy", "alpha_sha256", "beta_sha256", "wall_seconds",
         "allocated_bytes", "gc_seconds", "observer_seconds", "solver_minus_observer_seconds",
         "retained_cache_array_bytes", "cache_summarysize_bytes", "maxrss_bytes",
         "covariance_energy_error", "covariance_fock_error"),
        [("first_sweep_preflight_only_not_scientific_endpoint", energy, array_sha(Cup),
          array_sha(Cdn), total.time, total.bytes, total.gctime, meter.observer_seconds,
          total.time - meter.observer_seconds, meter.boundary_payload,
          meter.boundary_summary, Sys.maxrss(), packet.covariance_energy,
          packet.covariance_fock)])
    write_diagnostics(orderdir, diagnostics)
end

opts = arguments()
mkpath(opts.outdir)
pidpath = joinpath(opts.outdir, "h10_observability_preflight.pid")
write(pidpath, string(getpid(), '\n'))
println("pid=", getpid(), " outdir=", opts.outdir)
flush(stdout)
BLAS.set_num_threads(1)
gate(Threads.nthreads() == 1, "launch with one Julia thread")
ENV["HFDMRG_BENCH_TIMING"] = "1"
metadata = [("pid", getpid()), ("mode", opts.mode), ("order", opts.order),
    ("command", join([Base.julia_cmd().exec[1], "--project=.", @__FILE__, ARGS...], ' ')),
    ("pwd", pwd()), ("host", gethostname()), ("julia", VERSION),
    ("julia_threads", Threads.nthreads()), ("blas_threads", BLAS.get_num_threads()),
    ("commit", readchomp(`git rev-parse HEAD`)), ("status", readchomp(`git status --short`)),
    ("runner_sha256", sha256_file(@__FILE__)), ("contract_sha256", CONTRACT_SHA),
    ("archive_manifest_sha256", MANIFEST_SHA), ("v6_sha256", V6_SHA)]
write_tsv(joinpath(opts.outdir, "environment.tsv"), ("key", "value"), metadata)
if opts.mode in (:smoke, :overhead)
    run_small(opts.outdir, opts.mode, opts.trials)
else
    packet = verify_inputs(opts)
    write_tsv(joinpath(opts.outdir, "input_check.tsv"), ("quantity", "value"),
        [("ordering", opts.order), ("contract_sha256", CONTRACT_SHA),
         ("archive_manifest_sha256", MANIFEST_SHA), ("v6_sha256", V6_SHA),
         ("ordered_v6_sha256", packet.ordered_v6_sha),
         ("start_alpha_sha256", INPUT_SHA["data/atom_local_af_s1_alpha.tsv"]),
         ("start_beta_sha256", INPUT_SHA["data/atom_local_af_s1_beta.tsv"]),
         ("permutation_sha256", INPUT_SHA["data/slice_permutations.tsv"]),
         ("ordering_covariance_energy_error", packet.covariance_energy),
         ("ordering_covariance_fock_error", packet.covariance_fock)])
    opts.mode === :h10 && run_h10(opts, packet)
end
write(joinpath(opts.outdir, "completed.status"), "exit_status\t0\n")
println("completed mode=", opts.mode, " order=", opts.order, " maxrss=", Sys.maxrss())
