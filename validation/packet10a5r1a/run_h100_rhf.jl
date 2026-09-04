source_path = joinpath(ARGS[1], "validation", "packet8f7a",
    "run_block_native_product.jl")
source = read(source_path, String)
main_start = findfirst("\nfunction product_main()", source)
main_start === nothing && error("frozen PPP product fixture changed")
include_string(Main, source[1:first(main_start)-1], source_path)

include(joinpath(@__DIR__, "..", "..", "src", "HFDMRG.jl"))
using .HFDMRG
using LinearAlgebra: BLAS
using TOML

function main()
    length(ARGS) == 2 || error("usage: PPP_ROOT OUTPUT.toml")
    output = ARGS[2]
    old_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(4)
    _, small_model, small_periodic, _ = fixture(10)
    small_h1 = HFDMRG.BandedOneBody(small_model.one_body.bands)
    small_operator = HFDMRG.UnitCellInteraction(small_periodic.nsite,
        small_periodic.phase_at_first, small_periodic.tail_transfer,
        small_periodic.residual_transfer, small_periodic.phase_source,
        small_periodic.phase_sink, small_periodic.local_triangle)
    small_atoms = local_atomic_intervals(small_model.layout)
    HFDMRG.solve_hfdmrg(small_h1, small_operator; state_maxdim=16,
        state_cutoff=1.0e-12, energy_tolerance=1.0e-11,
        maximum_half_sweeps=4, exact=false, atom_intervals=small_atoms,
        initialization=:dimer)

    construction = @timed fixture(100)
    _, model, periodic, _ = construction.value
    one_body = HFDMRG.BandedOneBody(model.one_body.bands)
    operator = HFDMRG.UnitCellInteraction(periodic.nsite,
        periodic.phase_at_first, periodic.tail_transfer,
        periodic.residual_transfer, periodic.phase_source,
        periodic.phase_sink, periodic.local_triangle)
    atom_intervals = local_atomic_intervals(model.layout)
    energy_tolerance = 1.0e-11
    control = HFDMRG.RHFStateControl(16, 1.0e-12; exact=false)
    workspace = HFDMRG.RHFNonlinearWorkspace(one_body, operator, 16)
    initialized = @timed HFDMRG.initialize_rhf_dimer_lifecycle(one_body,
        operator, control, workspace.lifecycle; forward=true,
        atom_intervals)
    lifecycle = initialized.value
    solved = @timed HFDMRG._solver_loop!(lifecycle, one_body, operator,
        control, energy_tolerance, 20, workspace)
    outcome = solved.value
    visits = outcome[5] + outcome[6] + outcome[7]
    groups = length(lifecycle.intervals)
    expected_visits = (groups - 1) + (length(outcome[3]) - 1) * (groups - 2)
    primary_preparations = workspace.lifecycle.prepared.preparation_calls
    primary_focks = workspace.lifecycle.prepared.fock_calls
    maximum_raw_increase = maximum((update.published_energy -
        update.baseline_energy for update in outcome[11]); init=0.0)
    maximum_envelope_excess = maximum((update.published_energy -
        update.baseline_energy - HFDMRG._rhf_energy_envelope(
            update.baseline_energy, update.center_rank)
        for update in outcome[11]); init=0.0)

    replay = Dict{String,Any}[]
    replay_before = HFDMRG._invariant_energy!(lifecycle, one_body, operator,
        workspace.lifecycle)
    replay_timing = @timed begin
        for _ = 1:2
            orientation = lifecycle.forward ? "physical" : "reflected"
            evidence = HFDMRG.run_fixed_half_sweep!(lifecycle, one_body,
                operator, control, workspace.lifecycle)
            after = HFDMRG._invariant_energy!(lifecycle, one_body, operator,
                workspace.lifecycle)
            push!(replay, Dict("orientation" => orientation,
                "visits" => length(evidence), "energy_before" => replay_before,
                "energy_after" => after,
                "energy_change" => after - replay_before))
            replay_before = after
        end
    end
    represented_energy = replay_before
    ppp_low = 13.913151458
    ppp_high = 15.136047125
    outcome[1] || error("H100 RHF did not converge")
    visits == expected_visits || error("recipient center schedule changed")
    outcome[5] == visits && outcome[6] == 0 && outcome[7] == 0 ||
        error("H100 RHF did not use one accepted full proposal per center")
    primary_preparations - length(outcome[3]) == 2visits ||
        error("redundant accepted-center preparation")
    primary_focks - length(outcome[3]) == visits ||
        error("redundant accepted-center Fock construction")
    workspace.eigensolution_calls == visits ||
        error("redundant accepted-center eigensolution")
    workspace.trial_energy_calls == visits ||
        error("redundant accepted-center true-energy evaluation")
    workspace.publication_calls == visits ||
        error("redundant accepted-center publication")
    maximum_envelope_excess <= 0.0 || error("nonmonotone center update")
    lifecycle.fragment_ingestions == 50 || error("initialization count changed")
    lifecycle.original_rows_read == 0 || error("dense rows reread after initialization")
    solved.time <= 1.25ppp_low || error("H100 performance gate failed")
    all(abs(item["energy_change"]) <= energy_tolerance +
        HFDMRG._rhf_energy_envelope(item["energy_before"],
            lifecycle.root.rank) for item in replay) ||
        error("update-disabled orientation replay failed")
    record = Dict(
        "schema" => "hfdmrg-packet10a5r1a-h100-rhf-v2",
        "julia_version" => string(VERSION),
        "fixture_seconds" => construction.time,
        "initialization_seconds" => initialized.time,
        "initialization_allocation_bytes" => initialized.bytes,
        "solver_seconds" => solved.time,
        "solver_allocation_bytes" => solved.bytes,
        "replay_seconds" => replay_timing.time,
        "replay_allocation_bytes" => replay_timing.bytes,
        "represented_energy" => represented_energy,
        "converged" => outcome[1],
        "reason" => String(outcome[2]),
        "half_sweeps" => length(outcome[3]),
        "half_sweep_energies" => outcome[3],
        "full_sweep_changes" => outcome[4],
        "accepted_full" => outcome[5],
        "accepted_fallback" => outcome[6],
        "rejected_updates" => outcome[7],
        "maximum_state_rank" => outcome[8],
        "maximum_pair_rank" => outcome[9],
        "maximum_discarded_squared_weight" => outcome[10],
        "maximum_raw_published_energy_increase" => maximum_raw_increase,
        "maximum_monotonicity_envelope_excess" => maximum_envelope_excess,
        "energy_tolerance" => energy_tolerance,
        "traversal_groups" => groups,
        "visited_centers" => visits,
        "expected_schedule_visits" => expected_visits,
        "center_preparations" => primary_preparations,
        "fock_constructions" => primary_focks,
        "eigensolutions" => workspace.eigensolution_calls,
        "post_factorization_true_energy_evaluations" =>
            workspace.trial_energy_calls,
        "publications" => workspace.publication_calls,
        "transactions" => workspace.transaction_calls,
        "preparations_per_visited_center_excluding_half_invariants" =>
            (primary_preparations - length(outcome[3])) / visits,
        "focks_per_visited_center_excluding_half_invariants" =>
            (primary_focks - length(outcome[3])) / visits,
        "eigensolutions_per_visited_center" =>
            workspace.eigensolution_calls / visits,
        "true_energies_per_visited_center" =>
            workspace.trial_energy_calls / visits,
        "publications_per_visited_center" =>
            workspace.publication_calls / visits,
        "ppp_visits" => 800,
        "ppp_preparations_per_visited_center" => 2408 / 800,
        "ppp_focks_per_visited_center" => 1.0,
        "ppp_eigensolutions_per_visited_center" => 1.0,
        "ppp_true_energies_per_visited_center" => 1.0,
        "solver_ratio_to_ppp_low" => solved.time / ppp_low,
        "solver_ratio_to_ppp_high" => solved.time / ppp_high,
        "update_disabled_replays" => replay,
        "fragment_ingestions" => lifecycle.fragment_ingestions,
        "later_dense_row_reads" => lifecycle.original_rows_read,
        "center_workspace_bytes" =>
            HFDMRG.fast_rhf_center_storage_bytes(workspace.lifecycle.prepared),
        "peak_rss_bytes" => Int(Sys.maxrss()),
        "blas_threads" => BLAS.get_num_threads())
    open(output, "w") do io
        TOML.print(io, record; sorted=true)
    end
    BLAS.set_num_threads(old_threads)
    println(record)
end

main()
