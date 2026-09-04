source_path = joinpath(ARGS[1], "validation", "packet8f7a",
    "run_block_native_product.jl")
source = read(source_path, String)
main_start = findfirst("\nfunction product_main()", source)
main_start === nothing && error("frozen PPP product fixture changed")
include_string(Main, source[1:first(main_start)-1], source_path)

include(joinpath(@__DIR__, "..", "..", "src", "HFDMRG.jl"))
using .HFDMRG
using LinearAlgebra: BLAS, norm
using TOML

function projectors(result, one_body, operator)
    lifecycle = deepcopy(result.state.lifecycle)
    capacity = max(result.maximum_alpha_rank, result.maximum_beta_rank, 1)
    workspace = HFDMRG.UHFLifecycleWorkspace(one_body, operator, capacity,
        capacity)
    alpha = zeros(operator.sites, lifecycle.root.alpha.total_occupied)
    beta = zeros(operator.sites, lifecycle.root.beta.total_occupied)
    HFDMRG.reconstruct_uhf_terminal!(alpha, beta, lifecycle, workspace)
    alpha * transpose(alpha), beta * transpose(beta)
end

function run_case(one_body, operator, atom_intervals, orientation, exact;
        spin_swapped=false)
    timed = @timed HFDMRG.solve_hfdmrg(one_body, operator; spin=:uhf,
        initialization=:atomic_neel, state_maxdim=16,
        state_cutoff=1.0e-12, energy_tolerance=1.0e-11,
        maximum_half_sweeps=20, starting_orientation=orientation, exact,
        atom_intervals, spin_swapped)
    result = timed.value
    result.converged || error("H20 UHF case did not converge")
    result.state.lifecycle.fragment_ingestions == 20 || error(
        "H20 initialization count changed")
    result.state.lifecycle.original_rows_read == 0 || error(
        "H20 dense rows reread after initialization")
    Pa, Pb = projectors(result, one_body, operator)
    record = Dict("orientation" => String(orientation), "exact" => exact,
        "spin_swapped" => spin_swapped, "seconds" => timed.time,
        "allocation_bytes" => timed.bytes,
        "represented_energy" => result.represented_energy,
        "half_sweeps" => result.half_sweeps,
        "accepted_full" => result.accepted_full,
        "accepted_fallback" => result.accepted_fallback,
        "rejected_updates" => result.rejected_updates,
        "replay_energy_change" => result.replay_energy_change,
        "replay_alpha_projector_change" =>
            result.replay_alpha_projector_change,
        "replay_beta_projector_change" =>
            result.replay_beta_projector_change,
        "alpha_projector_envelope" => result.alpha_projector_envelope,
        "beta_projector_envelope" => result.beta_projector_envelope,
        "maximum_alpha_rank" => result.maximum_alpha_rank,
        "maximum_beta_rank" => result.maximum_beta_rank,
        "maximum_alpha_pair_rank" => result.maximum_alpha_pair_rank,
        "maximum_beta_pair_rank" => result.maximum_beta_pair_rank,
        "later_dense_row_reads" => result.state.lifecycle.original_rows_read)
    record, Pa, Pb
end

function main()
    length(ARGS) == 2 || error("usage: PPP_ROOT OUTPUT.toml")
    output = ARGS[2]
    old_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(4)
    try
        construction = @timed fixture(20)
        _, model, periodic, _ = construction.value
        one_body = HFDMRG.BandedOneBody(model.one_body.bands)
        operator = HFDMRG.UnitCellInteraction(periodic.nsite,
            periodic.phase_at_first, periodic.tail_transfer,
            periodic.residual_transfer, periodic.phase_source,
            periodic.phase_sink, periodic.local_triangle)
        atom_intervals = local_atomic_intervals(model.layout)
        cases = Dict{String,Any}()
        projectors_by_case = Dict{String,Tuple{Matrix{Float64},Matrix{Float64}}}()
        for exact in (false, true), orientation in (:physical, :reflected)
            key = string(exact ? "exact" : "finite", "_", orientation)
            record, Pa, Pb = run_case(one_body, operator, atom_intervals,
                orientation, exact)
            cases[key] = record
            projectors_by_case[key] = (Pa, Pb)
        end
        swapped, Sa, Sb = run_case(one_body, operator, atom_intervals,
            :physical, false; spin_swapped=true)
        cases["finite_physical_spin_swapped"] = swapped
        projectors_by_case["finite_physical_spin_swapped"] = (Sa, Sb)

        for prefix in ("finite", "exact")
            physical = cases[prefix * "_physical"]
            reflected = cases[prefix * "_reflected"]
            abs(physical["represented_energy"] -
                reflected["represented_energy"]) <= 1.0e-9 || error(
                "H20 orientation energies leave the finite envelope")
        end
        ordinary = cases["finite_physical"]
        abs(ordinary["represented_energy"] - swapped["represented_energy"]) <=
            1.0e-9 || error("H20 spin-swap energy changed")
        Pa, Pb = projectors_by_case["finite_physical"]
        norm(Pa - Sb) <= 2.0e-5 && norm(Pb - Sa) <= 2.0e-5 || error(
            "H20 spin-swap projectors changed basin")

        record = Dict("schema" => "hfdmrg-packet10a5r1b-h20-uhf-v1",
            "julia_version" => string(VERSION),
            "fixture_seconds" => construction.time,
            "blas_threads" => BLAS.get_num_threads(), "cases" => cases,
            "finite_orientation_energy_spread" => abs(
                cases["finite_physical"]["represented_energy"] -
                cases["finite_reflected"]["represented_energy"]),
            "exact_orientation_energy_spread" => abs(
                cases["exact_physical"]["represented_energy"] -
                cases["exact_reflected"]["represented_energy"]),
            "spin_swap_alpha_projector_distance" => norm(Pa - Sb),
            "spin_swap_beta_projector_distance" => norm(Pb - Sa))
        open(output, "w") do io
            TOML.print(io, record; sorted=true)
        end
        println(record)
    finally
        BLAS.set_num_threads(old_threads)
    end
end

main()
