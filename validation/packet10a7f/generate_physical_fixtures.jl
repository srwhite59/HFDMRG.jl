using DensityMPOLab
using GaussletBases: sliced_h1, sliced_h1_bandwidth, sliced_vee
using SHA: sha256
using TOML

const EXPECTED_PPP = "0bcd4f11831a98ae4bf4fb3d09eac5968df9d96c"
const EXPECTED_GAUSSLET = "c680bbbcce5e7ff9ef6405c60fec0d3c2fabb2f0"
const BASIS_FINGERPRINT =
    "73cc87798a9c78e7cc0f41ec56fc078bf8ccb8500fa6acb9655da8941bc788c7"

hex(bytes) = bytes2hex(sha256(bytes))

function require_sha(path, expected)
    actual = hex(read(path))
    actual == expected || error("pinned source checksum changed: $path")
end

function source_prefix(ppp_root)
    path = joinpath(ppp_root, "validation", "packet8f7a",
        "run_block_native_product.jl")
    source = read(path, String)
    start = findfirst("\nfunction product_main()", source)
    start === nothing && error("frozen PPP fixture runner changed")
    include_string(Main, source[1:first(start)-1], path)
end

length(ARGS) == 3 && source_prefix(abspath(ARGS[1]))

function array_record(name, array, offset)
    bytes = reinterpret(UInt8, vec(array))
    Dict("name" => name, "offset_bytes" => offset,
        "element_count" => length(array),
        "dimensions" => collect(size(array)), "sha256" => hex(bytes)), bytes
end

function producer_matrix(matrix)
    sites = size(matrix, 1)
    size(matrix, 2) == sites || error("producer interaction must be square")
    [Float64(matrix[first, second]) for first in 1:sites, second in 1:sites]
end

function reference_record(ppp_root, atoms, spin)
    name = "h$(atoms)_$(spin).toml"
    record = TOML.parsefile(joinpath(ppp_root, "validation", "packet6a4",
        name))
    Dict("ppp_packet6a4_total_ha" => record["global_audit"]["energy"],
        "ppp_record_generation_commit" => record["git"]["commit"],
        "ppp_record_path" => "validation/packet6a4/$name",
        "ppp_energy_tolerance_ha" => record["controls"]["energy_tolerance"])
end

function make_case(ppp_root, output_root, atoms)
    chain, model, periodic, report = Base.invokelatest(fixture, atoms)
    chain.diagnostics.basis_fingerprint == BASIS_FINGERPRINT || error(
        "physical basis fingerprint changed")
    sites = model.layout.nsites
    h1_view = sliced_h1(chain)
    bandwidth = sliced_h1_bandwidth(chain)
    bands = Matrix{Float64}(model.one_body.bands)
    size(bands) == (bandwidth + 1, sites) || error(
        "producer/adapter H1 dimensions changed")
    for offset in 0:bandwidth, first in 1:sites-offset
        bands[offset + 1, first] == h1_view[first, first + offset] || error(
            "producer/adapter H1 value changed")
    end
    vee = sliced_vee(chain)
    arrays = (bands, Vector{Float64}(periodic.tail_transfer),
        Matrix{Float64}(periodic.residual_transfer),
        Matrix{Float64}(periodic.phase_source),
        Matrix{Float64}(periodic.phase_sink),
        Vector{Float64}(periodic.local_triangle), producer_matrix(vee))
    names = ("one_body_bands", "tail_transfer", "residual_transfer",
        "phase_source", "phase_sink", "local_triangle",
        "producer_vee")
    filename = "h$(atoms).f64le"
    path = joinpath(output_root, filename)
    records = Dict{String,Any}[]
    open(path, "w") do io
        offset = 0
        for (name, array) in zip(names, arrays)
            record, bytes = array_record(name, array, offset)
            write(io, bytes)
            push!(records, record)
            offset += length(bytes)
        end
    end
    intervals = local_atomic_intervals(model.layout)
    Dict("atoms" => atoms, "sites" => sites,
        "occupied_per_spin" => atoms ÷ 2,
        "phase_at_first" => periodic.phase_at_first,
        "atom_starts" => first.(intervals), "atom_stops" => last.(intervals),
        "nuclear_repulsion_ha" => model.nuclear_repulsion,
        "element_type" => "Float64", "byte_order" => "little",
        "data_file" => filename, "data_sha256" => hex(read(path)),
        "data_bytes" => filesize(path), "arrays" => records,
        "compact_operator" => Dict("tail_rank" => periodic.tail_rank,
            "residual_rank" => periodic.residual_rank,
            "maximum_entry_error" => report.maximum_entry_error,
            "finite_chain_frobenius_error" =>
                report.finite_chain_frobenius_error,
            "requested_accuracy" => report.requested_accuracy,
            "training_sites" => report.training_nsites),
        "references" => Dict("rhf" => reference_record(ppp_root, atoms,
            "rhf"), "uhf" => reference_record(ppp_root, atoms, "uhf")))
end

function main()
    length(ARGS) == 3 || error(
        "usage: generate_physical_fixtures.jl PPP_ROOT GAUSSLET_ROOT OUTPUT_ROOT")
    ppp_root, gausslet_root, output_root = abspath.(ARGS)
    ENDIAN_BOM == 0x04030201 || error("fixture generation requires little endian")
    require_sha(joinpath(ppp_root, "validation", "packet6a4", "run_anchor.jl"),
        "613274e5331518e78fd7d7fbefb1c9baf7020c1648ff1613b4a0e7c2b71b29f5")
    require_sha(joinpath(ppp_root, "validation", "packet8f7a",
        "run_block_native_product.jl"),
        "5e0ed8f59b34cce852125911514aa4cbfc4dd2a1a8c397de565d1c6eacbc0cbe")
    require_sha(joinpath(ppp_root, "src", "SlicedHydrogenModels.jl"),
        "cb7a8b227f1c715234283f926eea63657ce8f4fc39cd93fe0bb88474071715b3")
    require_sha(joinpath(gausslet_root, "src", "sliced_hydrogen_chain.jl"),
        "35597a94781e2555d379e9b193c0be9db53d019ec14476684a7edd3dd67cbec3")
    require_sha(joinpath(gausslet_root, "LICENSE"),
        "be9b446af13b3c0df4e86d13a4481311f17b5309cfd50248ccf0d19180e4b19d")
    project = TOML.parsefile(joinpath(ppp_root, "Project.toml"))
    project["sources"]["GaussletBases"]["rev"] == EXPECTED_GAUSSLET ||
        error("frozen PPP GaussletBases revision changed")
    mkpath(output_root)
    cases = Dict("h10" => make_case(ppp_root, output_root, 10),
        "h20" => make_case(ppp_root, output_root, 20))
    manifest = Dict("schema" => "hfdmrg-physical-hydrogen-chain-v1",
        "binary_format" => "contiguous-f64le-v1",
        "license" => "MIT; Copyright (c) 2026 Steven R. White",
        "provenance" => Dict("ppp_oracle_commit" => EXPECTED_PPP,
            "gausslet_commit" => EXPECTED_GAUSSLET,
            "basis_fingerprint" => BASIS_FINGERPRINT,
            "geometry" => Dict("separation_bohr" => 3.6,
                "sites_per_atom" => 18, "spacing_bohr" => 0.2,
                "phase_bohr" => 0.0, "left_padding_bohr" => 2.0,
                "right_padding_bohr" => 2.0, "alignment" => "atom_centered",
                "boundary" => "open", "transverse_mode" =>
                    "periodic_template", "basis" => "H/STO-6G"),
            "one_body_tolerance" => 1.0e-8,
            "interaction_tolerance" => 1.0e-5,
            "energy_unit" => "hartree", "length_unit" => "bohr",
            "generated_data_terms" => "HFDMRG MIT distribution"),
        "cases" => cases)
    open(joinpath(output_root, "manifest.toml"), "w") do io
        TOML.print(io, manifest; sorted=true)
    end
end

main()
