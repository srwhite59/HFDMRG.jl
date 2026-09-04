base_path = joinpath(@__DIR__, "..", "packet10a5r2",
    "run_h1000_qualification.jl")
base = read(base_path, String)
main_start = findfirst("\nfunction main()", base)
main_start === nothing && error("R2 qualification support changed")
include_string(Main, base[1:first(main_start)-1], base_path)

using LinearAlgebra: BLAS
using TOML

const FRESH_PPP_SECONDS = 95.051920417

function compile_uhf_path()
    _, _, one_body, operator, atom_intervals = make_inputs(10)
    HFDMRG.solve_hfdmrg(one_body, operator; spin=:uhf,
        initialization=:atomic_neel, state_maxdim=16,
        state_cutoff=1.0e-12, energy_tolerance=1.0e-11,
        maximum_half_sweeps=4, starting_orientation=:physical,
        exact=false, atom_intervals)
    nothing
end

function main_r2b()
    length(ARGS) == 3 || error(
        "usage: OUTPUT.toml LOG_PATH PPP_ROOT")
    output, log_path = ARGS[1], ARGS[2]
    pid_path = log_path * ".pid"
    mkpath(dirname(output)); mkpath(dirname(log_path))
    log_io = open(log_path, "w")
    logline(message) = begin
        println(message); flush(stdout)
        println(log_io, message); flush(log_io)
    end
    open(pid_path, "w") do io
        println(io, getpid())
    end
    old_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(4)
    try
        logline("pid=$(getpid()) julia=$(VERSION) compile_start")
        compile_uhf_path()
        construction = @timed make_inputs(1000)
        _, _, one_body, operator, atom_intervals = construction.value
        logline("fixture_ready sites=$(operator.sites) " *
            "seconds=$(construction.time)")
        qualified = solve_uhf(one_body, operator, atom_intervals)
        result, lifecycle, workspace, outcome, initialized, solved, replay,
            counts, ownership = qualified
        logline("solve_done seconds=$(solved.time) " *
            "halves=$(result.half_sweeps) " *
            "energy=$(result.represented_energy)")

        result.converged || error("H1000 compact UHF did not converge")
        abs(result.represented_energy - REPRESENTED[:uhf]) <= 1.0e-9 ||
            error("represented energy left accepted practical basin")
        solved.time <= 1.5FRESH_PPP_SECONDS ||
            error("H1000 compact UHF hard performance gate failed")
        lifecycle.fragment_ingestions == 1000 ||
            error("block-native initialization count changed")
        lifecycle.original_rows_read == 0 ||
            error("dense rows were read after initialization")

        stationary = uhf_replay(lifecycle, one_body, operator;
            spin_swapped=false)
        stationary_seconds = stationary["seconds"]
        logline("stationary_replay_done " *
            "seconds=$stationary_seconds")
        swapped = uhf_replay(lifecycle, one_body, operator;
            spin_swapped=true)
        swapped_seconds = swapped["seconds"]
        logline("spin_swapped_replay_done " *
            "seconds=$swapped_seconds")

        maximum_raw_increase = maximum((update.published_energy -
            update.baseline_energy for update in outcome.updates); init=0.0)
        record = Dict{String,Any}(
            "schema" => "hfdmrg-packet10a5r2b-h1000-uhf-v1",
            "julia_version" => string(VERSION),
            "blas_threads" => BLAS.get_num_threads(),
            "atoms" => 1000, "sites" => operator.sites,
            "initialization" => "atomic_neel",
            "state_maxdim" => 16, "state_cutoff" => 1.0e-12,
            "energy_tolerance" => 1.0e-11,
            "converged" => result.converged,
            "reason" => String(result.reason),
            "represented_energy" => result.represented_energy,
            "accepted_represented_energy" => REPRESENTED[:uhf],
            "represented_energy_error" =>
                result.represented_energy - REPRESENTED[:uhf],
            "half_sweeps" => result.half_sweeps,
            "full_sweep_changes" => outcome.changes,
            "maximum_raw_published_energy_increase" => maximum_raw_increase,
            "initialization_seconds" => initialized.time,
            "initialization_allocation_bytes" => initialized.bytes,
            "nonlinear_seconds" => solved.time,
            "nonlinear_allocation_bytes" => solved.bytes,
            "fresh_frozen_ppp_nonlinear_seconds" => FRESH_PPP_SECONDS,
            "recipient_to_ppp_time_ratio" =>
                solved.time / FRESH_PPP_SECONDS,
            "performance_target_passed" =>
                solved.time <= 1.25FRESH_PPP_SECONDS,
            "fixture_seconds" => construction.time,
            "fixed_physical_reflected_replay" => replay,
            "stationary_nonlinear_replay" => stationary,
            "spin_swapped_stationary_replay" => swapped,
            "call_counts" => counts,
            "fragment_ingestions" => lifecycle.fragment_ingestions,
            "later_dense_row_reads" => lifecycle.original_rows_read,
            "collection_count" => 1, "root_count" => 1,
            "collection_summary_bytes" =>
                Base.summarysize(lifecycle.collection),
            "root_summary_bytes" => Base.summarysize(lifecycle.root),
            "nonlinear_workspace_summary_bytes" =>
                Base.summarysize(workspace),
            "lifecycle_storage_bytes" => Dict(string(key) => value
                for (key, value) in pairs(ownership)),
            "center_workspace_bytes" =>
                HFDMRG.fast_uhf_center_storage_bytes(
                    workspace.lifecycle.prepared),
            "peak_rss_bytes" => Int(Sys.maxrss()),
            "maximum_alpha_rank" => result.maximum_alpha_rank,
            "maximum_beta_rank" => result.maximum_beta_rank,
            "maximum_alpha_pair_rank" => result.maximum_alpha_pair_rank,
            "maximum_beta_pair_rank" => result.maximum_beta_pair_rank,
            "maximum_alpha_discarded_squared_weight" =>
                result.maximum_alpha_discarded_squared_weight,
            "maximum_beta_discarded_squared_weight" =>
                result.maximum_beta_discarded_squared_weight,
            "producer_audit_performed" => false,
            "producer_residual_evaluated" => false,
            "global_hf_minimum_claim" => false)
        open(output, "w") do io
            TOML.print(io, record; sorted=true)
        end
        logline("qualification_written output=$output")
    finally
        BLAS.set_num_threads(old_threads)
        close(log_io)
        rm(pid_path; force=true)
    end
end

main_r2b()
