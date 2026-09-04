source_path = joinpath(ARGS[3], "validation", "packet8f7a",
    "run_block_native_product.jl")
source = read(source_path, String)
dispatch = findlast("\nproduct_main()", source)
dispatch === nothing && error("frozen PPP fixture changed")
include_string(Main, source[1:first(dispatch)-1], source_path)

include(joinpath(@__DIR__, "..", "..", "src", "HFDMRG.jl"))
using .HFDMRG
using LinearAlgebra: BLAS
using Profile
using SHA: sha256
using TOML

@inline elapsed_seconds(start::UInt64) = (time_ns() - start) * 1.0e-9

function addtime!(times, name, start)
    times[name] = get(times, name, 0.0) + elapsed_seconds(start)
end

function build_case(atoms)
    construction = @timed fixture(atoms)
    _, model, periodic, _ = construction.value
    one_body = HFDMRG.BandedOneBody(model.one_body.bands)
    operator = HFDMRG.UnitCellInteraction(periodic.nsite,
        periodic.phase_at_first, periodic.tail_transfer,
        periodic.residual_transfer, periodic.phase_source,
        periodic.phase_sink, periodic.local_triangle)
    intervals = local_atomic_intervals(model.layout)
    control = HFDMRG.UHFStateControl(HFDMRG.RHFStateControl(16,
        1.0e-12; exact=false))
    workspace = HFDMRG.UHFNonlinearWorkspace(one_body, operator, 16, 16)
    initialized = @timed HFDMRG.initialize_uhf_neel_lifecycle(one_body,
        operator, control, workspace.lifecycle; forward=true,
        atom_intervals=intervals)
    (; construction, one_body, operator, control, workspace,
        lifecycle=initialized.value, initialized)
end

function profiled_half!(case)
    lifecycle, one_body, operator = case.lifecycle, case.one_body,
        case.operator
    control, workspace = case.control, case.workspace
    times = Dict{String,Float64}()
    visits = 0
    starting = lifecycle.forward
    while lifecycle.forward == starting
        start = time_ns()
        HFDMRG._preflight_uhf_nonlinear(lifecycle, one_body, operator,
            control, workspace)
        left, right = HFDMRG._uhf_lifecycle_blocks(lifecycle,
            workspace.lifecycle)
        cell = lifecycle.next_center
        alpha_root, beta_root = lifecycle.root.alpha, lifecycle.root.beta
        ma, mb = alpha_root.rank, beta_root.rank
        addtime!(times, "preflight_and_block_lookup", start)

        start = time_ns()
        HFDMRG.prepare_uhf_lifecycle_center!(workspace.lifecycle.prepared,
            left, right, cell, lifecycle.intervals[cell],
            lifecycle.intervals[cell+1], one_body, operator)
        addtime!(times, "incoming_center_preparation", start)

        start = time_ns()
        incoming = HFDMRG.uhf_center_energy_fock!(
            workspace.lifecycle.prepared,
            @view(alpha_root.covariance[1:ma, 1:ma]),
            @view(beta_root.covariance[1:mb, 1:mb]), alpha_root.occupied,
            beta_root.occupied)
        addtime!(times, "incoming_fock_and_energy", start)

        start = time_ns()
        copyto!(@view(workspace.alpha_incoming[1:ma, 1:ma]),
            @view(alpha_root.covariance[1:ma, 1:ma]))
        copyto!(@view(workspace.beta_incoming[1:mb, 1:mb]),
            @view(beta_root.covariance[1:mb, 1:mb]))
        copyto!(@view(workspace.alpha_fock[1:ma, 1:ma]),
            @view(workspace.lifecycle.prepared.alpha.fock[1:ma, 1:ma]))
        copyto!(@view(workspace.beta_fock[1:mb, 1:mb]),
            @view(workspace.lifecycle.prepared.beta.fock[1:mb, 1:mb]))
        HFDMRG._lowest_projector!(workspace.alpha_proposal,
            workspace.alpha_fock, ma, alpha_root.occupied)
        HFDMRG._lowest_projector!(workspace.beta_proposal,
            workspace.beta_fock, mb, beta_root.occupied)
        workspace.eigensolution_calls += 2
        addtime!(times, "proposal_copies_and_eigensolutions", start)

        start = time_ns()
        HFDMRG._copy_uhf_root!(workspace.trial_root, lifecycle.root)
        copyto!(@view(workspace.trial_root.alpha.covariance[1:ma, 1:ma]),
            @view(workspace.alpha_proposal[1:ma, 1:ma]))
        copyto!(@view(workspace.trial_root.beta.covariance[1:mb, 1:mb]),
            @view(workspace.beta_proposal[1:mb, 1:mb]))
        addtime!(times, "trial_root_copy", start)

        start = time_ns()
        candidate = HFDMRG._uhf_build_candidate!(workspace.trial_root,
            lifecycle, left, right, one_body, operator, control,
            workspace.lifecycle, workspace.trial_alpha_work,
            workspace.trial_beta_work)
        addtime!(times, "factor_entry_build_and_root_move", start)

        start = time_ns()
        next_left, next_right = HFDMRG._uhf_candidate_blocks(candidate,
            lifecycle, workspace.lifecycle)
        next_cell = candidate.next_center
        HFDMRG.prepare_uhf_lifecycle_center!(workspace.lifecycle.prepared,
            next_left, next_right, next_cell,
            lifecycle.intervals[next_cell], lifecycle.intervals[next_cell+1],
            one_body, operator)
        addtime!(times, "candidate_center_preparation", start)

        start = time_ns()
        trial = candidate.root
        proposal_energy = HFDMRG.uhf_center_true_energy!(
            workspace.lifecycle.prepared,
            @view(trial.alpha.covariance[1:trial.alpha.rank,
                1:trial.alpha.rank]),
            @view(trial.beta.covariance[1:trial.beta.rank,
                1:trial.beta.rank]), trial.alpha.occupied,
            trial.beta.occupied).total
        workspace.trial_energy_calls += 1
        addtime!(times, "candidate_true_energy", start)
        proposal_energy <= incoming.total + HFDMRG._rhf_energy_envelope(
            incoming.total, ma + mb) || error(
            "first-half probe unexpectedly entered fallback")

        start = time_ns()
        HFDMRG._publish_uhf_advance_candidate!(lifecycle, candidate)
        workspace.publication_calls += 1
        lifecycle.center_visits += 1
        workspace.transaction_calls += 1
        addtime!(times, "publication", start)
        visits += 1
    end
    lifecycle.half_sweeps += 1
    start = time_ns()
    terminal_energy = HFDMRG._uhf_invariant_energy!(lifecycle, one_body,
        operator, workspace.lifecycle)
    addtime!(times, "terminal_invariant", start)
    (; times, visits, terminal_energy)
end

function main()
    length(ARGS) == 4 || error(
        "usage: OUTPUT.toml PROFILE.txt PPP_ROOT LOG_PATH")
    output, profile_path, _, log_path = ARGS
    pid_path = log_path * ".pid"
    open(pid_path, "w") do io
        println(io, getpid())
    end
    atexit(() -> rm(pid_path; force=true))
    println("pid=", getpid(), " pid_path=", pid_path)
    BLAS.set_num_threads(4)
    small = build_case(10)
    profiled_half!(small)
    GC.gc()
    case = build_case(1000)
    Profile.clear(); Profile.init(n=20_000_000, delay=0.001)
    measured = @timed begin
        result = nothing
        Profile.@profile result = profiled_half!(case)
        result
    end
    open(profile_path, "w") do io
        Profile.print(io; format=:flat, sortedby=:count, mincount=2)
    end
    result = measured.value
    record = Dict(
        "schema" => "hfdmrg-packet10a5r2a-recipient-h1000-half-profile-v1",
        "julia_version" => string(VERSION), "blas_threads" => 4,
        "atoms" => 1000, "sites" => case.operator.sites,
        "visits" => result.visits,
        "terminal_energy" => result.terminal_energy,
        "half_sweep_seconds_with_sampling" => measured.time,
        "half_sweep_allocation_bytes" => measured.bytes,
        "phase_seconds" => result.times,
        "initialization_seconds" => case.initialized.time,
        "collection_bytes" => Base.summarysize(case.lifecycle.collection),
        "workspace_bytes" => Base.summarysize(case.workspace),
        "peak_rss_bytes" => Int(Sys.maxrss()),
        "profile_path" => profile_path,
        "profile_sha256" => bytes2hex(sha256(read(profile_path))),
        "production_source_changed" => false)
    open(output, "w") do io
        TOML.print(io, record; sorted=true)
    end
    open(log_path, "w") do io
        println(io, "pid=$(getpid())")
        println(io, "seconds=$(measured.time) visits=$(result.visits) " *
            "energy=$(result.terminal_energy)")
    end
    rm(pid_path; force=true)
    println(record)
end

main()
