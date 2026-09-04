source_path = joinpath(ARGS[3], "validation", "packet8f7a",
    "run_block_native_product.jl")
source = read(source_path, String)
dispatch = findlast("\nproduct_main()", source)
dispatch === nothing && error("frozen PPP block-native runner changed")
include_string(Main, source[1:first(dispatch)-1], source_path)

using Dates: UTC, now
using LinearAlgebra: BLAS
using SHA: sha256
using TOML

function logline(io, message)
    println(message); flush(stdout)
    println(io, message); flush(io)
end

function build_case(atoms::Int)
    construction = @timed fixture(atoms)
    _, model, periodic, report = construction.value
    one_body = model.one_body
    operator = VIEW._UnitCellInteractionView(periodic)
    reflected_one = CM.reflect_one_body(one_body)
    reflected_operator = VIEW._reflect(operator)
    centers, _ = SWEEP.terminal_complete_two_block_centers(periodic)
    intervals = local_atomic_intervals(model.layout)
    state = STATE._StateLinkControl(min(STATE_MAXDIM, atoms ÷ 2),
        STATE_CUTOFF)
    controls = (state, state)
    fragments = @timed FIXED._uhf_atomic_neel_fragments(one_body, intervals)
    provenance = UInt(hash((:packet10a5r2a, :frozen_ppp, atoms,
        :physical)))
    initialized = @timed FIXED._initialize_local_product_composite(one_body,
        operator, centers, fragments.value, controls, provenance)
    (; construction, report, one_body, operator, reflected_one,
        reflected_operator, centers, controls, fragments, initialized,
        provenance)
end

function compile_paths()
    case = build_case(10)
    run_nonlinear!(:uhf, case.initialized.value.lifecycle, case.centers,
        case.controls, case.one_body, case.operator, case.reflected_one,
        case.reflected_operator, case.provenance, 1)
    nothing
end

function main()
    length(ARGS) == 4 || error(
        "usage: OUTPUT.toml LOG_PATH FROZEN_PPP_ROOT FROZEN_COMMIT")
    output, log_path, frozen_root, frozen_commit = ARGS
    readchomp(`git -C $frozen_root rev-parse HEAD`) == frozen_commit ||
        error("PPP snapshot is not at the frozen donor commit")
    mkpath(dirname(output)); mkpath(dirname(log_path))
    pid_path = log_path * ".pid"
    open(pid_path, "w") do io
        println(io, getpid())
    end
    log_io = open(log_path, "w")
    old_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(4)
    process_start = time()
    try
        logline(log_io, "pid=$(getpid()) julia=$(VERSION) compile_start")
        VERSION == v"1.12.6" || error("qualification requires Julia 1.12.6")
        compile_paths()
        GC.gc()
        wall_context_before = readchomp(`uptime`)
        utc_start = string(now(UTC), "Z")
        case = build_case(1000)
        lifecycle = case.initialized.value.lifecycle
        initialization = case.initialized.value.evidence
        logline(log_io, "fixture_ready sites=$(case.operator.nsite) " *
            "seconds=$(case.construction.time) init=$(case.initialized.time)")
        nonlinear = run_nonlinear!(:uhf, lifecycle, case.centers,
            case.controls, case.one_body, case.operator,
            case.reflected_one, case.reflected_operator, case.provenance, 10)
        result = nonlinear.result
        logline(log_io, "solve_done seconds=$(nonlinear.seconds) " *
            "halves=$(result.half_sweeps) energy=$(result.energy)")
        used = nonlinear.evidence
        full_steps = sum(item["full_steps"] for item in used; init=0)
        geodesic_steps = sum(item["geodesic_steps"] for item in used; init=0)
        no_updates = sum(item["no_variational_updates"] for item in used;
            init=0)
        true_trials = sum(item["true_energy_trials"] for item in used;
            init=0)
        record = Dict{String,Any}(
            "schema" => "hfdmrg-packet10a5r2a-fresh-frozen-ppp-h1000-uhf-v1",
            "frozen_ppp_commit" => frozen_commit,
            "frozen_ppp_snapshot" => frozen_root,
            "julia_version" => string(VERSION),
            "blas_threads" => BLAS.get_num_threads(),
            "atoms" => 1000, "sites" => case.operator.nsite,
            "initialization" => "block_native_atomic_neel",
            "state_maxdim" => STATE_MAXDIM,
            "state_cutoff" => STATE_CUTOFF,
            "energy_tolerance" => ENERGY_TOLERANCE,
            "subspace_tolerance" => SUBSPACE_TOLERANCE,
            "converged" => result.converged,
            "outcome" => String(result.outcome),
            "represented_electronic_energy" => result.energy,
            "half_sweeps" => result.half_sweeps,
            "fixture_seconds" => case.construction.time,
            "fragment_construction_seconds" => case.fragments.time,
            "initialization_seconds" => case.initialized.time,
            "internal_nonlinear_seconds" => nonlinear.seconds,
            "internal_nonlinear_allocation_bytes" => nonlinear.bytes,
            "whole_process_wall_seconds" => time() - process_start,
            "workspace_bytes" => nonlinear.workspace_bytes,
            "collection_bytes" => Base.summarysize(lifecycle.collection),
            "root_bytes" => Base.summarysize(lifecycle.root),
            "peak_rss_bytes" => Int(Sys.maxrss()),
            "physical_memory_bytes" => Int(Sys.total_memory()),
            "free_memory_bytes_at_record" => Int(Sys.free_memory()),
            "fragment_ingestions" => initialization.fragment_ingestions,
            "global_dense_row_reads" => initialization.global_dense_row_reads,
            "dense_initializations" => lifecycle.dense_initializations,
            "later_dense_row_reads" => 0,
            "full_steps" => full_steps,
            "geodesic_steps" => geodesic_steps,
            "no_variational_updates" => no_updates,
            "true_energy_trials" => true_trials,
            "update_records" => length(nonlinear.updates),
            "half_sweep_evidence" => used,
            "wall_context_before" => wall_context_before,
            "wall_context_after" => readchomp(`uptime`),
            "utc_start" => utc_start,
            "utc_finish" => string(now(UTC), "Z"),
            "runner_sha256" => bytes2hex(sha256(read(@__FILE__))),
            "replay_performed" => false,
            "producer_audit_performed" => false)
        open(output, "w") do io
            TOML.print(io, record; sorted=true)
        end
        logline(log_io, "record_written output=$output")
    finally
        BLAS.set_num_threads(old_threads)
        close(log_io)
        rm(pid_path; force=true)
    end
end

main()
