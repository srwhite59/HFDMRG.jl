source_path = joinpath(ARGS[3], "validation", "packet8f7a",
    "run_block_native_product.jl")
source = read(source_path, String)
main_start = findfirst("\nfunction product_main()", source)
main_start === nothing && error("frozen PPP product fixture changed")
include_string(Main, source[1:first(main_start)-1], source_path)

include(joinpath(@__DIR__, "..", "..", "src", "HFDMRG.jl"))

using .HFDMRG
using GaussletBases: sliced_h1_bandwidth
using LinearAlgebra: BLAS
using Serialization: serialize
using SHA: sha256
using TOML

const REPRESENTED = Dict(:rhf => -2227.1283869230083,
    :uhf => -2282.2570732618196)
const PRODUCER = Dict(:rhf => -425.61104061524975,
    :uhf => -480.7394819304184)
const PPP_SECONDS = Dict(:rhf => 72.99095425, :uhf => 123.755389166)

function state_hash(state)
    io = IOBuffer()
    serialize(io, state)
    bytes2hex(sha256(take!(io)))
end

function make_inputs(atoms::Int)
    chain, model, periodic, _ = fixture(atoms)
    one_body = HFDMRG.BandedOneBody(model.one_body.bands)
    operator = HFDMRG.UnitCellInteraction(periodic.nsite,
        periodic.phase_at_first, periodic.tail_transfer,
        periodic.residual_transfer, periodic.phase_source,
        periodic.phase_sink, periodic.local_triangle)
    chain, model, one_body, operator, local_atomic_intervals(model.layout)
end

function compile_paths(spin::Symbol)
    chain, model, one_body, operator, atom_intervals = make_inputs(10)
    result = HFDMRG.solve_hfdmrg(one_body, operator; spin,
        initialization=spin === :rhf ? :dimer : :atomic_neel,
        state_maxdim=16, state_cutoff=1.0e-12,
        energy_tolerance=1.0e-11, maximum_half_sweeps=4,
        starting_orientation=:physical, exact=false, atom_intervals)
    producer = HFDMRG.ProducerHamiltonian(sliced_h1(chain),
        sliced_vee(chain); constant=model.nuclear_repulsion,
        h1_bandwidth=sliced_h1_bandwidth(chain))
    workspace = HFDMRG.ProducerEnergyWorkspace(result, producer)
    HFDMRG.producer_energy(result, producer, workspace)
    nothing
end

function rhf_replay(lifecycle, one_body, operator)
    state = deepcopy(lifecycle)
    control = HFDMRG.RHFStateControl(16, 1.0e-12; exact=false)
    workspace = HFDMRG.RHFNonlinearWorkspace(one_body, operator, 16)
    before = HFDMRG._invariant_energy!(state, one_body, operator,
        workspace.lifecycle)
    timed = @timed HFDMRG._solver_loop!(state, one_body, operator, control,
        1.0e-11, 4, workspace)
    outcome = timed.value
    Dict("half_sweeps" => length(outcome[3]),
        "starting_energy" => before, "ending_energy" => last(outcome[3]),
        "energy_change" => last(outcome[3]) - before,
        "full_sweep_changes" => outcome[4], "converged" => outcome[1],
        "seconds" => timed.time, "allocation_bytes" => timed.bytes)
end

function uhf_replay(lifecycle, one_body, operator; spin_swapped::Bool)
    state = deepcopy(lifecycle)
    spin_swapped && HFDMRG._spin_swap_uhf_lifecycle!(state, operator)
    control = HFDMRG.UHFStateControl(HFDMRG.RHFStateControl(16,
        1.0e-12; exact=false))
    workspace = HFDMRG.UHFNonlinearWorkspace(one_body, operator, 16, 16)
    before = HFDMRG._uhf_invariant_energy!(state, one_body, operator,
        workspace.lifecycle)
    timed = @timed HFDMRG._uhf_solver_loop!(state, one_body, operator,
        control, 1.0e-11, 4, workspace)
    outcome = timed.value
    Dict("spin_swapped" => spin_swapped,
        "half_sweeps" => length(outcome.energies),
        "starting_energy" => before,
        "ending_energy" => last(outcome.energies),
        "energy_change" => last(outcome.energies) - before,
        "full_sweep_changes" => outcome.changes,
        "converged" => outcome.converged,
        "seconds" => timed.time, "allocation_bytes" => timed.bytes)
end

function solve_rhf(one_body, operator, atom_intervals)
    control = HFDMRG.RHFStateControl(16, 1.0e-12; exact=false)
    workspace = HFDMRG.RHFNonlinearWorkspace(one_body, operator, 16)
    initialized = @timed HFDMRG.initialize_rhf_dimer_lifecycle(one_body,
        operator, control, workspace.lifecycle; forward=true, atom_intervals)
    lifecycle = initialized.value
    solved = @timed HFDMRG._solver_loop!(lifecycle, one_body, operator,
        control, 1.0e-11, 20, workspace)
    outcome = solved.value
    halves = length(outcome[3])
    visits = outcome[5] + outcome[6] + outcome[7]
    prepared = workspace.lifecycle.prepared
    counts = Dict("visited_centers" => visits,
        "center_preparations" => prepared.preparation_calls,
        "fock_constructions" => prepared.fock_calls,
        "eigensolutions" => workspace.eigensolution_calls,
        "true_energy_evaluations" => workspace.trial_energy_calls,
        "publications" => workspace.publication_calls,
        "transactions" => workspace.transaction_calls)
    counts["center_preparations"] - halves == 2visits ||
        error("redundant RHF center preparation")
    counts["fock_constructions"] - halves == visits ||
        error("redundant RHF Fock construction")
    counts["eigensolutions"] == visits || error("redundant RHF eigensolution")
    counts["true_energy_evaluations"] == visits ||
        error("redundant RHF true-energy evaluation")
    counts["publications"] == visits || error("redundant RHF publication")
    replay_before = HFDMRG._invariant_energy!(lifecycle, one_body, operator,
        workspace.lifecycle)
    fixed = @timed begin
        first = HFDMRG.run_fixed_half_sweep!(lifecycle, one_body, operator,
            control, workspace.lifecycle)
        middle = HFDMRG._invariant_energy!(lifecycle, one_body, operator,
            workspace.lifecycle)
        second = HFDMRG.run_fixed_half_sweep!(lifecycle, one_body, operator,
            control, workspace.lifecycle)
        after = HFDMRG._invariant_energy!(lifecycle, one_body, operator,
            workspace.lifecycle)
        (first=first, middle=middle, second=second, after=after)
    end
    replay = Dict("physical_or_reflected_first" =>
            fixed.value.first[1].forward ? "physical" : "reflected",
        "first_visits" => length(fixed.value.first),
        "first_energy_change" => fixed.value.middle - replay_before,
        "second_visits" => length(fixed.value.second),
        "second_energy_change" => fixed.value.after - fixed.value.middle,
        "total_energy_change" => fixed.value.after - replay_before,
        "seconds" => fixed.time, "allocation_bytes" => fixed.bytes)
    result = HFDMRG.HFDMRGResult(outcome[1], outcome[2], fixed.value.after,
        halves, outcome[5], outcome[6], outcome[7], outcome[8], outcome[9],
        outcome[10], outcome[4], fixed.value.after - replay_before, NaN,
        4sqrt(outcome[10]) + 8192eps(Float64)*max(lifecycle.root.rank, 1),
        HFDMRG.RHFCompactState(lifecycle))
    ownership = HFDMRG.lifecycle_storage_bytes(lifecycle, operator,
        workspace.lifecycle)
    result, lifecycle, workspace, outcome, initialized, solved, replay,
        counts, ownership
end

function solve_uhf(one_body, operator, atom_intervals)
    control = HFDMRG.UHFStateControl(HFDMRG.RHFStateControl(16,
        1.0e-12; exact=false))
    workspace = HFDMRG.UHFNonlinearWorkspace(one_body, operator, 16, 16)
    initialized = @timed HFDMRG.initialize_uhf_neel_lifecycle(one_body,
        operator, control, workspace.lifecycle; forward=true, atom_intervals)
    lifecycle = initialized.value
    solved = @timed HFDMRG._uhf_solver_loop!(lifecycle, one_body, operator,
        control, 1.0e-11, 20, workspace)
    outcome = solved.value
    halves = length(outcome.energies)
    visits = outcome.full + outcome.fallback + outcome.rejected
    prepared = workspace.lifecycle.prepared
    counts = Dict("visited_centers" => visits,
        "center_preparations" => prepared.preparation_calls,
        "fock_constructions" => prepared.fock_calls,
        "eigensolutions" => workspace.eigensolution_calls,
        "true_energy_evaluations" => workspace.trial_energy_calls,
        "publications" => workspace.publication_calls,
        "transactions" => workspace.transaction_calls)
    counts["center_preparations"] - halves == 2visits ||
        error("redundant UHF center preparation")
    counts["fock_constructions"] - halves == visits ||
        error("redundant UHF Fock construction")
    counts["eigensolutions"] == 2visits || error("redundant UHF eigensolution")
    counts["true_energy_evaluations"] == visits ||
        error("redundant UHF true-energy evaluation")
    counts["publications"] == visits || error("redundant UHF publication")
    replay_before = HFDMRG._uhf_invariant_energy!(lifecycle, one_body,
        operator, workspace.lifecycle)
    fixed = @timed begin
        first = HFDMRG.run_uhf_fixed_half_sweep!(lifecycle, one_body,
            operator, control, workspace.lifecycle)
        middle = HFDMRG._uhf_invariant_energy!(lifecycle, one_body, operator,
            workspace.lifecycle)
        second = HFDMRG.run_uhf_fixed_half_sweep!(lifecycle, one_body,
            operator, control, workspace.lifecycle)
        after = HFDMRG._uhf_invariant_energy!(lifecycle, one_body, operator,
            workspace.lifecycle)
        (first=first, middle=middle, second=second, after=after)
    end
    replay = Dict("physical_or_reflected_first" =>
            fixed.value.first[1].forward ? "physical" : "reflected",
        "first_visits" => length(fixed.value.first),
        "first_energy_change" => fixed.value.middle - replay_before,
        "second_visits" => length(fixed.value.second),
        "second_energy_change" => fixed.value.after - fixed.value.middle,
        "total_energy_change" => fixed.value.after - replay_before,
        "seconds" => fixed.time, "allocation_bytes" => fixed.bytes)
    alpha_envelope = 4sqrt(outcome.max_alpha_discarded) +
        8192eps(Float64)*max(lifecycle.root.alpha.rank, 1)
    beta_envelope = 4sqrt(outcome.max_beta_discarded) +
        8192eps(Float64)*max(lifecycle.root.beta.rank, 1)
    result = HFDMRG.UHFDMRGResult(outcome.converged, outcome.reason,
        fixed.value.after, halves, outcome.full, outcome.fallback,
        outcome.rejected, outcome.max_alpha, outcome.max_beta,
        outcome.max_alpha_pair, outcome.max_beta_pair,
        outcome.max_alpha_discarded, outcome.max_beta_discarded,
        outcome.changes, fixed.value.after - replay_before, NaN, NaN,
        alpha_envelope, beta_envelope, HFDMRG.UHFCompactState(lifecycle))
    ownership = HFDMRG.uhf_lifecycle_storage_bytes(lifecycle, operator,
        workspace.lifecycle)
    result, lifecycle, workspace, outcome, initialized, solved, replay,
        counts, ownership
end

function main()
    length(ARGS) == 4 || error(
        "usage: rhf|uhf OUTPUT.toml PPP_ROOT LOG_PATH")
    spin = Symbol(ARGS[1]); output, log_path = ARGS[2], ARGS[4]
    spin in (:rhf, :uhf) || error("spin must be rhf or uhf")
    pid_path = output * ".pid"
    mkpath(dirname(output)); mkpath(dirname(log_path))
    log_io = open(log_path, "w")
    logline(message) = begin
        println(message); flush(stdout); println(log_io, message); flush(log_io)
    end
    open(pid_path, "w") do io
        println(io, getpid())
    end
    old_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(4)
    try
        logline("pid=$(getpid()) spin=$spin julia=$(VERSION) compile_start")
        compile_paths(spin)
        construction = @timed make_inputs(1000)
        chain, model, one_body, operator, atom_intervals = construction.value
        logline("fixture_ready sites=$(operator.sites) seconds=$(construction.time)")
        qualified = spin === :rhf ? solve_rhf(one_body, operator,
            atom_intervals) : solve_uhf(one_body, operator, atom_intervals)
        result, lifecycle, workspace, outcome, initialized, solved, replay,
            counts, ownership = qualified
        logline("solve_done seconds=$(solved.time) halves=$(result.half_sweeps) energy=$(result.represented_energy)")

        result.converged || error("H1000 compact solve did not converge")
        abs(result.represented_energy - REPRESENTED[spin]) <= 1.0e-9 ||
            error("represented energy left accepted practical basin")
        solved.time <= 1.5PPP_SECONDS[spin] ||
            error("H1000 hard performance gate failed")
        lifecycle.fragment_ingestions == (spin === :rhf ? 500 : 1000) ||
            error("initialization ingestion count changed")
        lifecycle.original_rows_read == 0 ||
            error("dense rows were read after block-native initialization")

        stationary = spin === :rhf ? rhf_replay(lifecycle, one_body,
            operator) : uhf_replay(lifecycle, one_body, operator;
            spin_swapped=false)
        swapped = spin === :uhf ? uhf_replay(lifecycle, one_body, operator;
            spin_swapped=true) : Dict("performed" => false)
        stationary_seconds = stationary["seconds"]
        logline("stationary_replay_done seconds=$stationary_seconds")

        producer = HFDMRG.ProducerHamiltonian(sliced_h1(chain),
            sliced_vee(chain); constant=model.nuclear_repulsion,
            h1_bandwidth=sliced_h1_bandwidth(chain))
        producer_workspace = @timed HFDMRG.ProducerEnergyWorkspace(result,
            producer)
        before_hash = state_hash(result.state)
        producer_timed = @timed HFDMRG.producer_energy(result, producer,
            producer_workspace.value)
        after_hash = state_hash(result.state)
        value = producer_timed.value
        before_hash == after_hash || error("producer audit mutated state")
        abs(value.total - PRODUCER[spin]) <= max(value.error_envelope, 1.0e-9) ||
            error("producer energy left accepted practical basin")
        logline("producer_done seconds=$(producer_timed.time) energy=$(value.total)")

        changes = spin === :rhf ? outcome[4] : outcome.changes
        maximum_raw_increase = maximum((update.published_energy -
            update.baseline_energy for update in (spin === :rhf ?
            outcome[11] : outcome.updates)); init=0.0)
        record = Dict{String,Any}(
            "schema" => "hfdmrg-packet10a5r2-h1000-v1",
            "spin" => String(spin), "julia_version" => string(VERSION),
            "blas_threads" => BLAS.get_num_threads(), "atoms" => 1000,
            "sites" => operator.sites,
            "initialization" => spin === :rhf ? "dimer" : "atomic_neel",
            "state_maxdim" => 16, "state_cutoff" => 1.0e-12,
            "energy_tolerance" => 1.0e-11,
            "converged" => result.converged, "reason" => String(result.reason),
            "represented_energy" => result.represented_energy,
            "accepted_represented_energy" => REPRESENTED[spin],
            "represented_energy_error" =>
                result.represented_energy - REPRESENTED[spin],
            "half_sweeps" => result.half_sweeps,
            "full_sweep_changes" => changes,
            "maximum_raw_published_energy_increase" => maximum_raw_increase,
            "initialization_seconds" => initialized.time,
            "initialization_allocation_bytes" => initialized.bytes,
            "nonlinear_seconds" => solved.time,
            "nonlinear_allocation_bytes" => solved.bytes,
            "frozen_ppp_nonlinear_seconds" => PPP_SECONDS[spin],
            "recipient_to_ppp_time_ratio" => solved.time / PPP_SECONDS[spin],
            "performance_target_passed" => solved.time <= 1.25PPP_SECONDS[spin],
            "fixture_seconds" => construction.time,
            "fixed_physical_reflected_replay" => replay,
            "stationary_nonlinear_replay" => stationary,
            "spin_swapped_stationary_replay" => swapped,
            "call_counts" => counts,
            "fragment_ingestions" => lifecycle.fragment_ingestions,
            "later_dense_row_reads" => lifecycle.original_rows_read,
            "collection_count" => 1, "root_count" => 1,
            "collection_summary_bytes" => Base.summarysize(lifecycle.collection),
            "root_summary_bytes" => Base.summarysize(lifecycle.root),
            "nonlinear_workspace_summary_bytes" => Base.summarysize(workspace),
            "lifecycle_storage_bytes" => Dict(string(key) => value
                for (key, value) in pairs(ownership)),
            "center_workspace_bytes" => spin === :rhf ?
                HFDMRG.fast_rhf_center_storage_bytes(
                    workspace.lifecycle.prepared) :
                HFDMRG.fast_uhf_center_storage_bytes(
                    workspace.lifecycle.prepared),
            "peak_rss_bytes" => Int(Sys.maxrss()),
            "producer_energy" => value.total,
            "accepted_producer_energy" => PRODUCER[spin],
            "producer_energy_error" => value.total - PRODUCER[spin],
            "producer_error_envelope" => value.error_envelope,
            "producer_workspace_seconds" => producer_workspace.time,
            "producer_workspace_bytes" => Base.summarysize(
                producer_workspace.value),
            "producer_evaluation_seconds" => producer_timed.time,
            "producer_evaluation_allocation_bytes" => producer_timed.bytes,
            "producer_projector_block_reads" => value.projector_block_reads,
            "producer_state_hash_before" => before_hash,
            "producer_state_hash_after" => after_hash,
            "producer_state_mutated" => false,
            "producer_residual_evaluated" => false,
            "producer_hamiltonian_stationarity_claim" => false,
            "global_hf_minimum_claim" => false)
        if spin === :rhf
            record["maximum_active_rank"] = result.maximum_state_rank
            record["maximum_pair_rank"] = result.maximum_pair_rank
            record["maximum_discarded_squared_weight"] =
                result.maximum_discarded_squared_weight
        else
            record["maximum_alpha_rank"] = result.maximum_alpha_rank
            record["maximum_beta_rank"] = result.maximum_beta_rank
            record["maximum_alpha_pair_rank"] = result.maximum_alpha_pair_rank
            record["maximum_beta_pair_rank"] = result.maximum_beta_pair_rank
            record["maximum_alpha_discarded_squared_weight"] =
                result.maximum_alpha_discarded_squared_weight
            record["maximum_beta_discarded_squared_weight"] =
                result.maximum_beta_discarded_squared_weight
        end
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

main()
