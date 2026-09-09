using HFDMRG
using LinearAlgebra: BLAS
using TOML

include(joinpath(@__DIR__, "hydrogen_chain_fixture.jl"))
using .HydrogenChainFixtures: load_hydrogen_chain_fixture

function result_record(result::HFDMRGResult)
    Dict("converged" => result.converged, "reason" => String(result.reason),
        "represented_energy_ha" => result.represented_energy,
        "half_sweeps" => result.half_sweeps,
        "accepted_full" => result.accepted_full,
        "accepted_fallback" => result.accepted_fallback,
        "rejected_updates" => result.rejected_updates,
        "maximum_state_rank" => result.maximum_state_rank,
        "maximum_pair_rank" => result.maximum_pair_rank,
        "maximum_discarded_squared_weight" =>
            result.maximum_discarded_squared_weight,
        "replay_energy_change_ha" => result.replay_energy_change,
        "replay_projector_change" => result.replay_projector_change,
        "projector_envelope" => result.projector_envelope,
        "requested_energy_tolerance_ha" =>
            result.requested_energy_tolerance,
        "effective_energy_tolerance_ha" =>
            result.effective_energy_tolerance,
        "stopping_tolerance_reason" =>
            String(result.stopping_tolerance_reason))
end

function result_record(result::UHFDMRGResult)
    Dict("converged" => result.converged, "reason" => String(result.reason),
        "represented_energy_ha" => result.represented_energy,
        "half_sweeps" => result.half_sweeps,
        "accepted_full" => result.accepted_full,
        "accepted_fallback" => result.accepted_fallback,
        "rejected_updates" => result.rejected_updates,
        "maximum_alpha_rank" => result.maximum_alpha_rank,
        "maximum_beta_rank" => result.maximum_beta_rank,
        "maximum_alpha_pair_rank" => result.maximum_alpha_pair_rank,
        "maximum_beta_pair_rank" => result.maximum_beta_pair_rank,
        "maximum_alpha_discarded_squared_weight" =>
            result.maximum_alpha_discarded_squared_weight,
        "maximum_beta_discarded_squared_weight" =>
            result.maximum_beta_discarded_squared_weight,
        "replay_energy_change_ha" => result.replay_energy_change,
        "replay_alpha_projector_change" =>
            result.replay_alpha_projector_change,
        "replay_beta_projector_change" => result.replay_beta_projector_change,
        "alpha_projector_envelope" => result.alpha_projector_envelope,
        "beta_projector_envelope" => result.beta_projector_envelope,
        "requested_energy_tolerance_ha" =>
            result.requested_energy_tolerance,
        "effective_energy_tolerance_ha" =>
            result.effective_energy_tolerance,
        "stopping_tolerance_reason" =>
            String(result.stopping_tolerance_reason))
end

function producer_record(value::ProducerEnergy)
    Dict("one_body_alpha_ha" => value.one_body_alpha,
        "one_body_beta_ha" => value.one_body_beta,
        "hartree_ha" => value.hartree,
        "exchange_alpha_ha" => value.exchange_alpha,
        "exchange_beta_ha" => value.exchange_beta,
        "constant_ha" => value.constant, "total_ha" => value.total,
        "error_envelope_ha" => value.error_envelope,
        "projector_block_reads" => value.projector_block_reads)
end

function main()
    length(ARGS) in (2, 3) || error(
        "usage: physical_hydrogen_chain.jl H10_OR_H20 RHF_OR_UHF [OUTPUT.toml]")
    case_name = Symbol(lowercase(ARGS[1]))
    spin = Symbol(lowercase(ARGS[2]))
    spin in (:rhf, :uhf) || error("spin must be rhf or uhf")
    fixture = load_hydrogen_chain_fixture(case_name)
    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(4)
    try
        solved = @timed solve_hfdmrg(fixture.one_body, fixture.interaction;
            spin, initialization=spin === :rhf ? :dimer : :atomic_neel,
            state_maxdim=16, state_cutoff=1.0e-12,
            energy_tolerance=1.0e-11, maximum_half_sweeps=20,
            starting_orientation=:physical, exact=false,
            atom_intervals=fixture.atom_intervals)
        result = solved.value
        observed = @timed producer_energy(result, fixture.producer)
        value = observed.value
        record = Dict("schema" => "hfdmrg-physical-hydrogen-result-v1",
            "case" => String(case_name), "spin" => String(spin),
            "julia_version" => string(VERSION),
            "blas_threads" => BLAS.get_num_threads(),
            "solver_seconds" => solved.time,
            "solver_allocation_bytes" => solved.bytes,
            "producer_seconds" => observed.time,
            "producer_allocation_bytes" => observed.bytes,
            "result" => result_record(result),
            "producer" => producer_record(value),
            "represented_total_with_nuclear_ha" =>
                result.represented_energy + value.constant,
            "represented_minus_producer_ha" =>
                result.represented_energy + value.constant - value.total,
            "reference" => fixture.metadata["references"])
        if length(ARGS) == 3
            open(ARGS[3], "w") do io
                TOML.print(io, record; sorted=true)
            end
        else
            TOML.print(stdout, record; sorted=true)
        end
        result.converged || error("physical compact solve did not converge")
    finally
        BLAS.set_num_threads(old_blas)
    end
end

main()
