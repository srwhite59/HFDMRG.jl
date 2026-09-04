source_path = joinpath(ARGS[3], "validation", "packet8f7a",
    "run_block_native_product.jl")
source = read(source_path, String)
dispatch = findlast("\nproduct_main()", source)
dispatch === nothing && error("frozen PPP runner changed")
include_string(Main, source[1:first(dispatch)-1], source_path)

using LinearAlgebra: BLAS
using TOML

function build_case(atoms)
    construction = @timed fixture(atoms)
    _, model, periodic, _ = construction.value
    one_body = model.one_body
    operator = VIEW._UnitCellInteractionView(periodic)
    reflected_one = CM.reflect_one_body(one_body)
    reflected_operator = VIEW._reflect(operator)
    centers, _ = SWEEP.terminal_complete_two_block_centers(periodic)
    intervals = local_atomic_intervals(model.layout)
    state = STATE._StateLinkControl(min(STATE_MAXDIM, atoms ÷ 2),
        STATE_CUTOFF)
    controls = (state, state)
    fragments = FIXED._uhf_atomic_neel_fragments(one_body, intervals)
    provenance = UInt(hash((:packet10a5r2a, :ppp_half, atoms)))
    initialized = @timed FIXED._initialize_local_product_composite(one_body,
        operator, centers, fragments, controls, provenance)
    workspace = SWEEP._IndependentSpinStoredNonlinearWorkspace(
        initialized.value.lifecycle, centers, controls, one_body, operator)
    (; construction, one_body, operator, reflected_one, reflected_operator,
        centers, controls, provenance, initialized, workspace)
end

function one_half!(case)
    control = SWEEP._ControlledStoredSweepControl(maximum_full_sweeps=1,
        energy_tolerance=ENERGY_TOLERANCE,
        subspace_tolerance=SUBSPACE_TOLERANCE, strategy=:monotone_one_fock)
    updates = NamedTuple[]
    result = SWEEP._independent_spin_nonlinear_half_sweep!(
        case.initialized.value.lifecycle, case.centers, case.controls,
        case.one_body, case.operator, case.reflected_one,
        case.reflected_operator, case.provenance, case.workspace, NaN, false,
        control; update_evidence=updates)
    result, updates
end

function main()
    length(ARGS) == 4 || error(
        "usage: OUTPUT.toml LOG_PATH FROZEN_PPP_ROOT FROZEN_COMMIT")
    output, log_path, frozen_root, frozen_commit = ARGS
    pid_path = log_path * ".pid"
    open(pid_path, "w") do io
        println(io, getpid())
    end
    atexit(() -> rm(pid_path; force=true))
    println("pid=", getpid(), " pid_path=", pid_path)
    readchomp(`git -C $frozen_root rev-parse HEAD`) == frozen_commit ||
        error("PPP snapshot is not at the frozen commit")
    BLAS.set_num_threads(4)
    one_half!(build_case(10))
    GC.gc()
    case = build_case(1000)
    TG4.reset_timing_report!(); TG4.set_timing!(true)
    measured = @timed one_half!(case)
    TG4.set_timing!(false)
    result, updates = measured.value
    record = Dict(
        "schema" => "hfdmrg-packet10a5r2a-frozen-ppp-h1000-half-profile-v1",
        "frozen_ppp_commit" => frozen_commit,
        "julia_version" => string(VERSION), "blas_threads" => 4,
        "atoms" => 1000, "sites" => case.operator.nsite,
        "visits" => length(updates), "terminal_energy" => result[1],
        "full_steps" => result[4], "geodesic_steps" => result[5],
        "no_variational_updates" => result[6],
        "fock_builds" => result[7], "eigensolutions" => result[8],
        "true_energy_trials" => result[9],
        "half_sweep_seconds" => measured.time,
        "half_sweep_allocation_bytes" => measured.bytes,
        "initialization_seconds" => case.initialized.time,
        "collection_bytes" => Base.summarysize(
            case.initialized.value.lifecycle.collection),
        "workspace_bytes" => Base.summarysize(case.workspace),
        "peak_rss_bytes" => Int(Sys.maxrss()),
        "timing" => timing_summary(TG4.current_timing_report().roots),
        "production_source_changed" => false)
    open(output, "w") do io
        TOML.print(io, record; sorted=true)
    end
    open(log_path, "w") do io
        println(io, "pid=$(getpid()) seconds=$(measured.time) " *
            "visits=$(length(updates)) energy=$(result[1])")
    end
    rm(pid_path; force=true)
    println(record)
end

main()
