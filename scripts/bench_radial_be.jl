using HFDMRG, LinearAlgebra, Serialization, SHA, Sockets
const EXPECTED_SHA = "fba29651c02cde3cdde90850d90ba30a1abd120ebfa8ac8b64cf2e5a92fea629"
const PACKET_KEYS = (:metadata, :H, :V6, :dims, :seed_up, :seed_dn, :oracle_up, :oracle_dn, :Zop,
    :radial_centers, :radial_weights, :angular_offsets, :oracle_energy, :oracle_residual, :oracle_fock_angle,
    :oracle_s2, :oracle_zspin, :oracle_parity, :random_parity, :current_backend_known_mismatch)
const BRANCH = (projector = 1e-5, minsv = 1 - 1e-10, energy = 1e-8, fingerprint = 1e-6)
gate(ok, message) = ok || error(message)
projector(C) = C * C'
principal(C, reference) = minimum(svdvals(reference' * C))
projector_distance(C, R) = hypot(norm(C - R * (R' * C)), norm(R - C * (C' * R)))
function arguments()
    options = Dict("packet" => "", "outdir" => "", "modes" => "check", "trials" => "3", "max-sweeps" => "50", "energy-cutoff" => "1e-11")
    for arg in ARGS; key, value = split(lstrip(arg, '-'), '='; limit = 2); gate(haskey(options, key), "unknown option --$key"); options[key] = value; end
    gate(all(!isempty(options[k]) for k in ("packet", "outdir")), "use --packet=PATH --outdir=PATH [--modes=...] [--trials=3] [--max-sweeps=50] [--energy-cutoff=1e-11]")
    modes = Set(Symbol.(split(options["modes"], ',')))
    gate(issubset(modes, Set((:check, :fixture, :parity, :timing, :converge, :alignment, :rss, :all))), "unknown mode; use check,fixture,parity,timing,converge,alignment,rss,all")
    gate(:all ∉ modes || modes == Set((:all,)), "all must be the only mode and excludes rss"); :all in modes && (modes = Set((:fixture, :parity, :timing, :converge, :alignment)))
    trials, max_sweeps, energy_cutoff = parse(Int, options["trials"]), parse(Int, options["max-sweeps"]), parse(Float64, options["energy-cutoff"])
    gate(:timing ∉ modes || trials >= 3, "timing requires --trials=3 or greater")
    gate(max_sweeps > 0 && isfinite(energy_cutoff) && energy_cutoff >= 0, "max-sweeps must be positive and energy-cutoff finite and nonnegative (zero disables outer stopping)")
    gate(:rss ∉ modes || (modes == Set((:rss,)) && abspath(PROGRAM_FILE) == abspath(@__FILE__)), "rss must be the only mode in a fresh standalone Julia process (all excludes rss)")
    abspath(options["packet"]), abspath(options["outdir"]), modes, trials, max_sweeps, energy_cutoff
end
function write_tsv(path, names, rows)
    tmp = path * ".tmp.$(getpid())"; open(tmp, "w") do io
        println(io, join(names, '\t'))
        for row in rows; println(io, join((replace(string(x), '\t' => ' ') for x in row), '\t')); end
    end; mv(tmp, path; force = true)
end
atomic_serialize(path, value) = (tmp = path * ".tmp.$(getpid())"; serialize(tmp, value); mv(tmp, path; force = true))
# Independent physical six-index contraction; it does not call either backend.
function physical_fock(packet, Cup, Cdn)
    H, V6 = packet.H, packet.V6; Dup, Ddn = projector(Cup), projector(Cdn)
    Gup, Gdn = zeros(size(H)), zeros(size(H)); nj, ns = packet.metadata.nj, packet.metadata.nr
    @inbounds for n = 1:ns, m = 1:ns, a = 1:nj, b = 1:nj, c = 1:nj, d = 1:nj
        p, r = (n - 1) * nj + a, (n - 1) * nj + b
        q, s = (m - 1) * nj + c, (m - 1) * nj + d
        v = V6[a, b, c, d, n, m]; direct = (Dup[q, s] + Ddn[q, s]) * v
        Gup[p, r] += direct; Gdn[p, r] += direct
        Gup[p, s] -= Dup[r, q] * v; Gdn[p, s] -= Ddn[r, q] * v
    end
    Fup, Fdn = H + Gup, H + Gdn
    energy = 0.5dot(Cup, (Fup + H) * Cup) + 0.5dot(Cdn, (Fdn + H) * Cdn)
    Fup, Fdn, Dup, Ddn, energy, Gup, Gdn
end
function diagnostics(packet, Cup, Cdn, reported = NaN)
    Fup, Fdn, Dup, Ddn, energy = physical_fock(packet, Cup, Cdn)
    FCup, FCdn = Fup * Cup, Fdn * Cdn; Rup, Rdn = FCup - Cup * (Cup' * FCup), FCdn - Cdn * (Cdn' * FCdn)
    nup, ndn = size(Cup, 2), size(Cdn, 2); ou, od = principal(Cup, packet.oracle_up), principal(Cdn, packet.oracle_dn)
    Gu, Gd, cross = Cup' * Cup, Cdn' * Cdn, Cup' * Cdn; overlap = sum(abs2, cross); sz = (nup - ndn) / 2
    (; energy, energy_mismatch = abs(reported - energy), residual = hypot(norm(Rup), norm(Rdn)),
        oracle_min_up = ou, oracle_min_dn = od, oracle_angle = max(acos(clamp(ou, 0, 1)), acos(clamp(od, 0, 1))),
        oracle_projector_up = projector_distance(Cup, packet.oracle_up), oracle_projector_dn = projector_distance(Cdn, packet.oracle_dn),
        energy_delta_oracle = abs(energy - packet.oracle_energy), orthogonality = max(norm(Gu - I), norm(Gd - I)),
        idempotency = max(norm(Gu * Gu - Gu), norm(Gd * Gd - Gd)),
        s2 = sz * (sz + 1) + ndn - overlap, zspin = dot(Cup, packet.Zop * Cup) - dot(Cdn, packet.Zop * Cdn),
        density_spin_norm = sqrt(max(0.0, sum(abs2, Gu) + sum(abs2, Gd) - 2overlap)), spin_overlap = overlap)
end
function determinant_fock_angle(packet, Cup, Cdn)
    Fup, Fdn = physical_fock(packet, Cup, Cdn)[1:2]; Eup, Edn = eigen(Symmetric(Fup)).vectors, eigen(Symmetric(Fdn)).vectors
    fu, fd = principal(Cup, Eup[:, 1:size(Cup, 2)]), principal(Cdn, Edn[:, 1:size(Cdn, 2)])
    max(acos(clamp(fu, 0, 1)), acos(clamp(fd, 0, 1)))
end
mutable struct Meter{B}
    backend::B; layout::HFDMRG.SliceLayout; phase::Int; mark_ns::UInt64; mark_bytes::Int64; wall::Vector{Float64}
    bytes::Vector{Int64}; counts::Dict{Symbol, Vector{Int}}; routes::Vector{Dict{Symbol, Float64}}
end
struct MeterWindow{W, M}; window::W; meter::M; end
function Meter(backend, phases)
    HFDMRG._bench_timing_reset!(); keys = (:init, :absorb, :window, :fock)
    Meter(backend, backend.layout, 1, time_ns(), Base.gc_bytes(), zeros(phases), zeros(Int64, phases), Dict(k => zeros(Int, phases) for k in keys), [Dict{Symbol, Float64}() for _ = 1:phases])
end
counted!(m, key, f) = (m.counts[key][m.phase] += 1; f())
function cut!(m)
    m.wall[m.phase] = (time_ns() - m.mark_ns) / 1e9; m.bytes[m.phase] = Base.gc_bytes() - m.mark_bytes
    m.routes[m.phase] = HFDMRG._bench_timing_snapshot(); HFDMRG._bench_timing_reset!()
end
advance!(m) = (m.phase += 1; m.mark_ns = time_ns(); m.mark_bytes = Base.gc_bytes())
HFDMRG.vee_init_block(side, ra, raV, phi, m::Meter) = counted!(m, :init, () -> HFDMRG.vee_init_block(side, ra, raV, phi, m.backend))
HFDMRG.vee_absorb_block(side, old, cra, Po, Pc, phi, raV, m::Meter) = counted!(m, :absorb, () -> HFDMRG.vee_absorb_block(side, old, cra, Po, Pc, phi, raV, m.backend))
function HFDMRG.vee_window(left, right, center, m::Meter)
    m.phase == 1 && (cut!(m); advance!(m))
    MeterWindow(counted!(m, :window, () -> HFDMRG.vee_window(left, right, center, m.backend)), m)
end
HFDMRG.vee_add_fock_r!(F, D, w::MeterWindow) = counted!(w.meter, :fock, () -> HFDMRG.vee_add_fock_r!(F, D, w.window))
HFDMRG.vee_add_fock!(Fu, Fd, Du, Dd, w::MeterWindow) = counted!(w.meter, :fock, () -> HFDMRG.vee_add_fock!(Fu, Fd, Du, Dd, w.window))
function solve_metered(packet, route; sweeps = 1, aligned = true, converge = false, energy_cutoff = 1e-11, callback = nothing)
    layout = HFDMRG.SliceLayout(packet.dims)
    raw = route == :cached ? HFDMRG.SlicedBasisBackendCached(layout, packet.V6) : HFDMRG.SlicedBasisBackend(layout, packet.V6)
    meter, snapshots = Meter(raw, sweeps + 1), Any[]
    observer = info -> begin
        cut!(meter); push!(snapshots, (sweep = info.sweep, energy = info.energy, converged = info.converged, psiup = copy(info.psiup), psidn = copy(info.psidn)))
        callback !== nothing && callback(info, meter)
        info.sweep < sweeps && advance!(meter); false
    end
    cutoff, scf = converge ? (energy_cutoff, 1e-12) : (0.0, 0.0)
    result = HFDMRG.solve_hfdmrg_core(packet.H, meter, packet.seed_up, packet.seed_dn; restricted = false,
        maxiter = sweeps, blocksize = 9, block_partition = aligned ? layout : nothing,
        cutoff, scf_cutoff = scf, observer, verbose = false)
    (; result, meter, snapshots)
end
function aligned_cached(m, phase; forced = false)
    route = m.routes[phase]
    m.counts[:window][phase] == 98 && (!forced || m.counts[:fock][phase] == 490) &&
        get(route, :cache_incremental_calls, -1) == 98 &&
        get(route, :cache_fallback_calls, -1) == 0 &&
        get(route, :window_local_calls, -1) == 98 &&
        get(route, :window_projection_calls, -1) == 0
end
function window_checks(packet)
    layout = HFDMRG.SliceLayout(packet.dims)
    projected, cached = HFDMRG.SlicedBasisBackend(layout, packet.V6), HFDMRG.SlicedBasisBackendCached(layout, packet.V6)
    both, rows = hcat(packet.seed_up, packet.seed_dn), Any[]
    for center_slice in (2, 26, 51)
        Lra = 1:layout.offs[center_slice]; Cra = layout.offs[center_slice] + 1:layout.offs[center_slice + 1]
        Rra = layout.offs[center_slice + 1] + 1:layout.offs[end]
        Lphi, Rphi = HFDMRG.getphi(both[Lra, :])[1], HFDMRG.getphi(both[Rra, :])[1]
        lp = HFDMRG.vee_init_block(:left, Lra, last(Lra) + 1:last(Rra), Lphi, projected); rp = HFDMRG.vee_init_block(:right, Rra, 1:first(Rra) - 1, Rphi, projected)
        lc = HFDMRG.vee_init_block(:left, Lra, last(Lra) + 1:last(Rra), Lphi, cached); rc = HFDMRG.vee_init_block(:right, Rra, 1:first(Rra) - 1, Rphi, cached)
        wp, wc = HFDMRG.vee_window(lp, rp, Cra, projected), HFDMRG.vee_window(lc, rc, Cra, cached)
        B = wp.B; Cup, Cdn = B' * packet.seed_up, B' * packet.seed_dn
        reconup, recondn = norm(B * Cup - packet.seed_up), norm(B * Cdn - packet.seed_dn)
        Dup, Ddn = projector(Cup), projector(Cdn)
        Gpu, Gpd, Gcu, Gcd = (zeros(size(Dup)) for _ = 1:4)
        HFDMRG.vee_add_fock!(Gpu, Gpd, Dup, Ddn, wp); HFDMRG.vee_add_fock!(Gcu, Gcd, Dup, Ddn, wc)
        Fup, Fdn = physical_fock(packet, B * Cup, B * Cdn)[1:2]; Hlocal = B' * packet.H * B
        Fpu, Fpd, Fcu, Fcd = Hlocal + Gpu, Hlocal + Gpd, Hlocal + Gcu, Hlocal + Gcd
        fullu, fulld = B' * Fup * B, B' * Fdn * B
        push!(rows, (center_slice, reconup, recondn, maximum(abs, Fpu - Fcu),
            maximum(abs, Fpd - Fcd), maximum(abs, Fpu - fullu), maximum(abs, Fpd - fulld)))
    end
    rows
end
packet_path, outdir, modes, trials, max_sweeps, energy_cutoff = arguments(); mkpath(outdir)
ENV["HFDMRG_BENCH_TIMING"] = "1"; sha = bytes2hex(open(SHA.sha256, packet_path))
gate(sha == EXPECTED_SHA, "packet SHA256 mismatch: $sha"); packet = deserialize(packet_path)
gate(keys(packet) == PACKET_KEYS, "packet schema mismatch")
layout = HFDMRG.SliceLayout(packet.dims); metadata = packet.metadata
gate(metadata.name == "Be-radial-Ylm-s015-lmax2-v1" && metadata.nr == 52 && metadata.lmax == 2 && metadata.nj == 9 && metadata.Nup == metadata.Ndn == 2, "fixture identity mismatch")
gate(size(packet.H) == (468, 468) && size(packet.V6) == (9, 9, 9, 9, 52, 52) &&
    packet.dims == fill(9, 52) && all(size(C) == (468, 2) for C in (packet.seed_up, packet.seed_dn, packet.oracle_up, packet.oracle_dn)), "fixture dimensions mismatch")
gate(size(packet.Zop) == (468, 468) && length(packet.radial_centers) == 52 &&
    length(packet.radial_weights) == 52 && packet.angular_offsets == [1, 2, 5, 10] && metadata.V6_permutation == (2, 3, 5, 6, 1, 4), "fixture order mismatch")
arrays = (packet.H, packet.V6, packet.seed_up, packet.seed_dn, packet.oracle_up, packet.oracle_dn, packet.Zop, packet.radial_centers, packet.radial_weights)
gate(all(A -> all(isfinite, A), arrays) && all(>(0), diff(packet.radial_centers)) &&
    all(>(0), packet.radial_weights), "fixture has nonfinite or unordered radial data")
gate(norm(packet.H - packet.H') <= 1e-11 * max(1, norm(packet.H)) &&
    norm(packet.Zop - packet.Zop') <= 1e-11 * max(1, norm(packet.Zop)), "fixture operator symmetry mismatch")
gate(all(isfinite, (packet.oracle_energy, packet.oracle_residual, packet.oracle_fock_angle, packet.oracle_s2, packet.oracle_zspin)) &&
    all(isfinite, values(packet.oracle_parity)) && all(isfinite, values(packet.random_parity)), "fixture has nonfinite oracle metadata")
gate(all(norm(C' * C - I) <= 1e-11 for C in (packet.seed_up, packet.seed_dn, packet.oracle_up, packet.oracle_dn)), "fixture orbitals are not orthonormal")
write(joinpath(outdir, "bench_radial_be.pid"), string(getpid(), '\n'))
command = join([Base.julia_cmd().exec[1], "--project=.", @__FILE__, ARGS...], ' ')
meta = [("pid", getpid()), ("command", command), ("pwd", pwd()), ("commit", readchomp(`git rev-parse HEAD`)),
    ("git_status", readchomp(`git status --short`)), ("host", gethostname()), ("cpu", join(unique(x.model for x in Sys.cpu_info()), ';')),
    ("julia", VERSION), ("julia_threads", Threads.nthreads()), ("blas_threads", BLAS.get_num_threads()), ("blas_config", BLAS.get_config()),
    ("active_project", Base.active_project()), ("packet", packet_path), ("sha256", sha), ("modes", join(sort!(collect(modes)), ',')),
    ("trials", trials), ("warmups", 1), ("max_convergence_sweeps", max_sweeps), ("outer_energy_cutoff", energy_cutoff), ("scf_cutoff", 1e-12), ("branch_projector_gate", BRANCH.projector), ("branch_minsv_gate", BRANCH.minsv),
    ("branch_energy_gate", BRANCH.energy), ("branch_fingerprint_gate", BRANCH.fingerprint), ("allocation_method", "Base.gc_bytes"), ("memory_method", "Sys.maxrss fresh process"), ("timing_statistic", "three or more raw trials; report median")]
for key in ("JULIA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "JULIA_DEPOT_PATH", "JULIA_PROJECT", "HFDMRG_BENCH_TIMING")
    push!(meta, ("env_" * key, get(ENV, key, "")))
end
write_tsv(joinpath(outdir, "environment.tsv"), ("key", "value"), meta)
foreach(row -> println(row[1], "=", row[2]), meta)
if :fixture in modes
    seed = diagnostics(packet, packet.seed_up, packet.seed_dn); oracle = diagnostics(packet, packet.oracle_up, packet.oracle_dn, packet.oracle_energy)
    oracle_fock_angle = determinant_fock_angle(packet, packet.oracle_up, packet.oracle_dn)
    gate(oracle.energy_mismatch <= 1e-10 && oracle.orthogonality <= 1e-11 && oracle.idempotency <= 1e-11, "oracle energy/density gate failed")
    gate(abs(oracle.residual - packet.oracle_residual) <= 1e-10 &&
        abs(oracle.s2 - packet.oracle_s2) <= 1e-12 && abs(oracle.zspin - packet.oracle_zspin) <= 1e-12,
        "oracle branch mismatch")
    gate(abs(oracle_fock_angle - packet.oracle_fock_angle) <= 1e-7, "oracle Fock-angle mismatch")
    write_tsv(joinpath(outdir, "fixture.tsv"), ("state", propertynames(seed)..., "recorded_fock_angle", "computed_fock_angle"), [("seed", values(seed)..., NaN, NaN), ("oracle", values(oracle)..., packet.oracle_fock_angle, oracle_fock_angle)])
end
if :parity in modes
    projection, cached = solve_metered(packet, :projection), solve_metered(packet, :cached); Up, Dp, Ep = projection.result; Uc, Dc, Ec = cached.result
    dp, dc = diagnostics(packet, Up, Dp, Ep), diagnostics(packet, Uc, Dc, Ec); eu, ed = projector_distance(Up, Uc), projector_distance(Dp, Dc)
    gate(abs(Ep - Ec) <= 1e-9 && max(eu, ed) <= 1e-8, "same-range parity gate failed")
    gate(max(dp.energy_mismatch, dc.energy_mismatch) <= 1e-10 &&
        max(dp.orthogonality, dp.idempotency, dc.orthogonality, dc.idempotency) <= 1e-11,
        "returned determinant gate failed")
    mp, mc = projection.meter, cached.meter; initial, sweep = mc.routes[1], mc.routes[2]
    gate(mp.counts[:window][2] == 98 && mp.counts[:fock][2] == 490 &&
        aligned_cached(mc, 2; forced = true), "forced-four count/route gate failed")
    gate(get(initial, :cache_init_calls, -1) == 2 && get(initial, :cache_incremental_calls, -1) == 49 &&
        get(initial, :cache_fallback_calls, -1) == 0, "cached initialization-route gate failed")
    windows = window_checks(packet)
    parity_header = ("energy_projection", "energy_cached", "energy_difference", "projector_up", "projector_dn",
        "projection_energy_mismatch", "cached_energy_mismatch", "projection_windows", "projection_focks", "cached_windows",
        "cached_focks", "cached_init", "cached_initial_incremental", "cached_sweep_incremental", "cached_initial_fallback",
        "cached_sweep_fallback", "cached_window_local", "cached_window_projection")
    write_tsv(joinpath(outdir, "parity.tsv"), parity_header, [(Ep, Ec, abs(Ep - Ec), eu, ed,
        dp.energy_mismatch, dc.energy_mismatch, mp.counts[:window][2], mp.counts[:fock][2],
        mc.counts[:window][2], mc.counts[:fock][2], get(initial, :cache_init_calls, -1),
        get(initial, :cache_incremental_calls, -1), get(sweep, :cache_incremental_calls, -1),
        get(initial, :cache_fallback_calls, -1), get(sweep, :cache_fallback_calls, -1),
        get(sweep, :window_local_calls, -1), get(sweep, :window_projection_calls, -1))])
    write_tsv(joinpath(outdir, "windows.tsv"), ("center_slice", "seed_reconstruction_up", "seed_reconstruction_dn",
        "projection_cached_total_fock_up", "projection_cached_total_fock_dn", "projection_full_total_fock_up", "projection_full_total_fock_dn"), windows)
    gate(maximum(maximum(row[2:3]) for row in windows) <= 1e-11, "retained-basis seed reconstruction gate failed")
    gate(maximum(maximum(row[4:end]) for row in windows) <= 1e-12, "window total-Fock gate failed")
    atomic_serialize(joinpath(outdir, "parity_checkpoint.jls"), (projection = projection.snapshots, cached = cached.snapshots))
end
function meter_rows(trial, m)
    labels = ("initialization", "first_sweep", "second_sweep")
    [(trial, labels[p], m.wall[p], m.bytes[p], (m.counts[k][p] for k in (:init, :absorb, :window, :fock))...,
      (get(m.routes[p], k, 0) for k in (:cache_init_calls, :cache_incremental_calls, :cache_fallback_calls, :window_local_calls, :window_projection_calls))...) for p = 1:3]
end
if :timing in modes
    solve_metered(packet, :cached; sweeps = 2); rows, checkpoints = Any[], Any[]
    for trial = 1:trials
        GC.gc(); run = solve_metered(packet, :cached; sweeps = 2)
        append!(rows, meter_rows(trial, run.meter)); push!(checkpoints, run.snapshots)
        gate(all(run.meter.wall[p] < 40 && run.meter.bytes[p] < 2 * 1024^3 for p = 2:3), "cached timing/allocation gate failed")
        gate(all(aligned_cached(run.meter, p; forced = true) for p = 2:3), "timing count/route gate failed")
    end
    header = ("trial", "stage", "wall_seconds", "gc_bytes", "init_calls", "absorb_calls", "window_calls", "fock_calls",
        "cache_init", "cache_incremental", "cache_fallback", "window_local", "window_projection")
    write_tsv(joinpath(outdir, "timing.tsv"), header, rows)
    atomic_serialize(joinpath(outdir, "timing_checkpoints.jls"), checkpoints)
end
if :converge in modes
    rows, lastdiag = Any[], Ref{Any}()
    callback = (info, m) -> begin
        d = diagnostics(packet, info.psiup, info.psidn, info.energy); lastdiag[] = d; route = m.routes[m.phase]
        gate(aligned_cached(m, m.phase), "convergence route-count gate failed")
        push!(rows, (info.sweep, info.energy, info.converged, values(d)..., m.wall[m.phase], m.bytes[m.phase],
            m.counts[:window][m.phase], m.counts[:fock][m.phase],
            get(route, :cache_incremental_calls, -1), get(route, :cache_fallback_calls, -1),
            get(route, :window_local_calls, -1), get(route, :window_projection_calls, -1)))
        write_tsv(joinpath(outdir, "convergence.tsv"), ("sweep", "reported_energy", "converged", propertynames(d)...,
            "wall_seconds", "gc_bytes", "windows", "focks", "cache_incremental", "cache_fallback", "window_local", "window_projection"), rows)
        atomic_serialize(joinpath(outdir, "convergence_checkpoint.jls"),
            (sweep = info.sweep, energy = info.energy, converged = info.converged,
             diagnostics = d, psiup = info.psiup, psidn = info.psidn))
    end
    run = solve_metered(packet, :cached; sweeps = max_sweeps, converge = true, energy_cutoff, callback)
    if energy_cutoff > 0
        gate(run.snapshots[end].converged, "matched convergence failed within $max_sweeps sweeps")
    else
        gate(length(run.snapshots) == max_sweeps && all(!s.converged for s in run.snapshots), "diagnostic trace did not complete exactly $max_sweeps sweeps")
    end
    gate(lastdiag[].energy_mismatch <= 1e-10 && lastdiag[].orthogonality <= 1e-11 && lastdiag[].idempotency <= 1e-11, "final determinant integrity gate failed")
end
if :alignment in modes
    numeric = solve_metered(packet, :projection; aligned = false); aligned, cached = solve_metered(packet, :projection), solve_metered(packet, :cached)
    gate(numeric.meter.counts[:window][2] == 100 && numeric.meter.counts[:fock][2] == 500 &&
        aligned.meter.counts[:window][2] == 98 && aligned.meter.counts[:fock][2] == 490 &&
        cached.meter.counts[:window][2] == 98 && cached.meter.counts[:fock][2] == 490,
        "alignment geometry/count gate failed")
    gate(aligned_cached(cached.meter, 2; forced = true), "alignment cached-route gate failed")
    rows, results = Any[], Dict{String, Any}()
    for (name, run) in (("numeric_projection", numeric), ("aligned_projection", aligned), ("aligned_cached", cached))
        U, D, E = run.result; d = diagnostics(packet, U, D, E); results[name] = (U, D, E, d)
        push!(rows, (name, E, run.meter.wall[1], run.meter.wall[2], run.meter.bytes[2],
            run.meter.counts[:window][2], run.meter.counts[:fock][2], d.energy_mismatch,
            d.oracle_projector_up, d.oracle_projector_dn, get(run.meter.routes[2], :cache_fallback_calls, 0),
            get(run.meter.routes[2], :window_local_calls, 0),
            get(run.meter.routes[2], :window_projection_calls, 0)))
    end
    gate(all(row[8] <= 1e-10 for row in rows), "alignment returned-energy gate failed")
    Ua, Da, Ea = results["aligned_projection"][1:3]; Uc, Dc, Ec = results["aligned_cached"][1:3]
    gate(abs(Ea - Ec) <= 1e-9 && max(projector_distance(Ua, Uc), projector_distance(Da, Dc)) <= 1e-8, "alignment projection/cached parity gate failed")
    write_tsv(joinpath(outdir, "alignment.tsv"), ("route", "energy", "initialization_seconds", "sweep_seconds", "gc_bytes",
        "windows", "fock_calls", "energy_mismatch", "oracle_projector_up", "oracle_projector_dn", "cache_fallback", "window_local", "window_projection"), rows)
    atomic_serialize(joinpath(outdir, "alignment_checkpoints.jls"), (numeric = numeric.snapshots, aligned = aligned.snapshots, cached = cached.snapshots))
end
if :rss in modes
    run = solve_metered(packet, :cached)
    gate(aligned_cached(run.meter, 2; forced = true), "RSS cached-route gate failed")
    rss = Sys.maxrss()
    write_tsv(joinpath(outdir, "rss.tsv"), ("maxrss_bytes", "maxrss_mib", "sweep_seconds", "gc_bytes"), [(rss, rss / 2.0^20, run.meter.wall[2], run.meter.bytes[2])])
end
println("completed modes=", join(sort!(collect(modes)), ','), " maxrss=", Sys.maxrss())
