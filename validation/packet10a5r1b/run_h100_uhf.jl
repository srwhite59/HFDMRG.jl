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
    pid_path = output * ".pid"
    mkpath(dirname(output))
    open(pid_path, "w") do io
        println(io, getpid())
    end
    old_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(4)
    try
        # Compile every production path outside the measured boundary.
        _, small_model, small_periodic, _ = fixture(10)
        small_h1 = HFDMRG.BandedOneBody(small_model.one_body.bands)
        small_operator = HFDMRG.UnitCellInteraction(small_periodic.nsite,
            small_periodic.phase_at_first, small_periodic.tail_transfer,
            small_periodic.residual_transfer, small_periodic.phase_source,
            small_periodic.phase_sink, small_periodic.local_triangle)
        small_atoms = local_atomic_intervals(small_model.layout)
        HFDMRG.solve_hfdmrg(small_h1, small_operator; spin=:uhf,
            initialization=:atomic_neel, state_maxdim=16,
            state_cutoff=1.0e-12, energy_tolerance=1.0e-11,
            maximum_half_sweeps=4, exact=false,
            atom_intervals=small_atoms)

        construction = @timed fixture(100)
        _, model, periodic, _ = construction.value
        one_body = HFDMRG.BandedOneBody(model.one_body.bands)
        operator = HFDMRG.UnitCellInteraction(periodic.nsite,
            periodic.phase_at_first, periodic.tail_transfer,
            periodic.residual_transfer, periodic.phase_source,
            periodic.phase_sink, periodic.local_triangle)
        atom_intervals = local_atomic_intervals(model.layout)
        control = HFDMRG.UHFStateControl(HFDMRG.RHFStateControl(16,
            1.0e-12; exact=false))
        workspace = HFDMRG.UHFNonlinearWorkspace(one_body, operator, 16, 16)
        initialized = @timed HFDMRG.initialize_uhf_neel_lifecycle(one_body,
            operator, control, workspace.lifecycle; forward=true,
            atom_intervals)
        lifecycle = initialized.value
        solved = @timed HFDMRG._uhf_solver_loop!(lifecycle, one_body,
            operator, control, 1.0e-11, 20, workspace)
        outcome = solved.value
        visits = outcome.full + outcome.fallback + outcome.rejected
        halves = length(outcome.energies)
        groups = length(lifecycle.intervals)
        expected_visits = (groups - 1) + (halves - 1) * (groups - 2)
        prepared = workspace.lifecycle.prepared
        maximum_envelope_excess = maximum((update.published_energy -
            update.baseline_energy - HFDMRG._rhf_energy_envelope(
            update.baseline_energy,
            update.alpha_center_rank + update.beta_center_rank)
            for update in outcome.updates); init=0.0)

        outcome.converged || error("H100 UHF did not converge")
        visits == expected_visits || error("recipient center schedule changed")
        outcome.full == visits && outcome.fallback == 0 &&
            outcome.rejected == 0 || error(
            "H100 UHF did not use one accepted full proposal per center")
        prepared.preparation_calls - halves == 2visits || error(
            "redundant accepted-center preparation")
        prepared.fock_calls - halves == visits || error(
            "redundant accepted-center Fock construction")
        workspace.eigensolution_calls == 2visits || error(
            "redundant accepted-center eigensolution")
        workspace.trial_energy_calls == visits || error(
            "redundant accepted-center true-energy evaluation")
        workspace.publication_calls == visits || error(
            "redundant accepted-center publication")
        maximum_envelope_excess <= 0.0 || error("nonmonotone center update")
        lifecycle.fragment_ingestions == 100 || error(
            "initialization fragment count changed")
        lifecycle.original_rows_read == 0 || error(
            "dense rows reread after initialization")

        record = Dict(
            "schema" => "hfdmrg-packet10a5r1b-h100-uhf-v1",
            "julia_version" => string(VERSION),
            "fixture_seconds" => construction.time,
            "initialization_seconds" => initialized.time,
            "initialization_allocation_bytes" => initialized.bytes,
            "solver_seconds" => solved.time,
            "solver_allocation_bytes" => solved.bytes,
            "represented_energy" => last(outcome.energies),
            "converged" => outcome.converged,
            "reason" => String(outcome.reason),
            "half_sweeps" => halves,
            "half_sweep_energies" => outcome.energies,
            "full_sweep_changes" => outcome.changes,
            "accepted_full" => outcome.full,
            "accepted_fallback" => outcome.fallback,
            "rejected_updates" => outcome.rejected,
            "maximum_alpha_rank" => outcome.max_alpha,
            "maximum_beta_rank" => outcome.max_beta,
            "maximum_alpha_pair_rank" => outcome.max_alpha_pair,
            "maximum_beta_pair_rank" => outcome.max_beta_pair,
            "maximum_alpha_discarded_squared_weight" =>
                outcome.max_alpha_discarded,
            "maximum_beta_discarded_squared_weight" =>
                outcome.max_beta_discarded,
            "maximum_monotonicity_envelope_excess" =>
                maximum_envelope_excess,
            "energy_tolerance" => 1.0e-11,
            "traversal_groups" => groups,
            "visited_centers" => visits,
            "expected_schedule_visits" => expected_visits,
            "center_preparations" => prepared.preparation_calls,
            "fock_constructions" => prepared.fock_calls,
            "eigensolutions" => workspace.eigensolution_calls,
            "post_factorization_true_energy_evaluations" =>
                workspace.trial_energy_calls,
            "publications" => workspace.publication_calls,
            "transactions" => workspace.transaction_calls,
            "preparations_per_visited_center_excluding_half_invariants" =>
                (prepared.preparation_calls - halves) / visits,
            "focks_per_visited_center_excluding_half_invariants" =>
                (prepared.fock_calls - halves) / visits,
            "eigensolutions_per_spin_per_visited_center" =>
                workspace.eigensolution_calls / (2visits),
            "true_energies_per_visited_center" =>
                workspace.trial_energy_calls / visits,
            "publications_per_visited_center" =>
                workspace.publication_calls / visits,
            "fragment_ingestions" => lifecycle.fragment_ingestions,
            "later_dense_row_reads" => lifecycle.original_rows_read,
            "center_workspace_bytes" =>
                HFDMRG.fast_uhf_center_storage_bytes(prepared),
            "lifecycle_storage" => Dict(string(key) => value
                for (key, value) in pairs(
                    HFDMRG.uhf_lifecycle_storage_bytes(lifecycle, operator,
                        workspace.lifecycle))),
            "workspace_summary_bytes" => Base.summarysize(workspace),
            "peak_rss_bytes" => Int(Sys.maxrss()),
            "blas_threads" => BLAS.get_num_threads())
        println(record)
        flush(stdout)
        open(output, "w") do io
            TOML.print(io, record; sorted=true)
        end
    finally
        BLAS.set_num_threads(old_threads)
        rm(pid_path; force=true)
    end
end

main()
