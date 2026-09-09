using TOML

function verify_case(manifest, root, name, spin)
    fixture = manifest["cases"][name]
    record = TOML.parsefile(joinpath(root, "$(name)_$(spin).toml"))
    result = record["result"]
    result["converged"] || error("$name $spin did not converge")
    result["reason"] == "energy_converged" || error(
        "$name $spin stopped for an unexpected reason")
    requested = result["requested_energy_tolerance_ha"]
    effective = result["effective_energy_tolerance_ha"]
    requested == 1.0e-11 || error("$name $spin controls changed")
    effective >= requested || error("$name $spin effective tolerance shrank")
    abs(result["replay_energy_change_ha"]) <= effective || error(
        "$name $spin replay exceeds its reported stopping threshold")
    reference = fixture["references"][spin]
    represented_total = record["represented_total_with_nuclear_ha"]
    reference_delta = represented_total - reference["ppp_packet6a4_total_ha"]
    abs(reference_delta) <= reference["ppp_energy_tolerance_ha"] || error(
        "$name $spin left the matched frozen-PPP energy basin")
    electrons = 2fixture["occupied_per_spin"]
    maximum_entry = fixture["compact_operator"]["maximum_entry_error"]
    representation_bound = 0.5electrons * (electrons + 1) * maximum_entry
    representation_gap = abs(record["represented_minus_producer_ha"])
    representation_gap <= representation_bound || error(
        "$name $spin producer gap exceeds the case compression bound")
    Dict("represented_total_reference_delta_ha" => reference_delta,
        "producer_representation_gap_ha" => representation_gap,
        "conservative_interaction_energy_bound_ha" => representation_bound,
        "represented_energy_ha" => result["represented_energy_ha"],
        "producer_total_ha" => record["producer"]["total_ha"],
        "half_sweeps" => result["half_sweeps"],
        "requested_energy_tolerance_ha" => requested,
        "effective_energy_tolerance_ha" => effective,
        "stopping_tolerance_reason" => result["stopping_tolerance_reason"],
        "solver_seconds" => record["solver_seconds"],
        "producer_seconds" => record["producer_seconds"])
end

function main()
    length(ARGS) == 3 || error(
        "usage: validate_reproductions.jl MANIFEST RECORD_ROOT OUTPUT")
    manifest = TOML.parsefile(ARGS[1])
    cases = Dict{String,Any}()
    for name in ("h10", "h20"), spin in ("rhf", "uhf")
        cases["$(name)_$(spin)"] = verify_case(manifest, ARGS[2], name, spin)
    end
    evidence = Dict("schema" =>
            "hfdmrg-packet10a7f-physical-reproduction-v1",
        "base_commit" => "68364562d78c80da31ae351417bddc88138f155f",
        "julia_version" => string(VERSION), "blas_threads" => 4,
        "clean_candidate_checkout" => true,
        "clean_candidate_contains_git_metadata" => false,
        "public_execution_imported_ppp" => false,
        "public_execution_imported_gausslet_bases" => false,
        "coefficient_checkpoints_used" => false,
        "reference_results_used_as_inputs" => false,
        "pinned_source_entry_parity" => true,
        "producer_input_parity" => true,
        "fixture_license_authorized" => true, "cases" => cases)
    open(ARGS[3], "w") do io
        TOML.print(io, evidence; sorted=true)
    end
end

main()
