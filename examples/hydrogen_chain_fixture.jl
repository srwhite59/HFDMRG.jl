module HydrogenChainFixtures

using HFDMRG
using SHA: sha256
using TOML

export DEFAULT_MANIFEST, load_hydrogen_chain_fixture

const DEFAULT_MANIFEST = joinpath(@__DIR__, "data", "hydrogen_chains_v1",
    "manifest.toml")
const _SCHEMA = "hfdmrg-physical-hydrogen-chain-v1"
const _FORMAT = "contiguous-f64le-v1"
const _ARRAY_NAMES = ("one_body_bands", "tail_transfer",
    "residual_transfer", "phase_source", "phase_sink", "local_triangle",
    "producer_vee")

_hex(bytes) = bytes2hex(sha256(bytes))

function _integer(value, label)
    value isa Integer || throw(ArgumentError("$label must be an integer"))
    Int(value)
end

function _dimensions(value, label)
    value isa Vector || throw(ArgumentError("$label dimensions must be a list"))
    dimensions = [_integer(item, label) for item in value]
    !isempty(dimensions) && all(>(0), dimensions) || throw(ArgumentError(
        "$label dimensions must be positive"))
    dimensions
end

function _element_count(dimensions, label)
    count = 1
    for dimension in dimensions
        count = try
            Base.Checked.checked_mul(count, dimension)
        catch
            throw(ArgumentError("$label dimensions overflow"))
        end
    end
    count
end

function _load_arrays(root, record)
    record["element_type"] == "Float64" || throw(ArgumentError(
        "fixture element type must be Float64"))
    record["byte_order"] == "little" || throw(ArgumentError(
        "fixture byte order must be little endian"))
    ENDIAN_BOM == 0x04030201 || throw(ArgumentError(
        "this fixture loader currently requires a little-endian host"))
    filename = String(record["data_file"])
    filename == basename(filename) || throw(ArgumentError(
        "fixture data_file must be a local basename"))
    path = joinpath(root, filename)
    bytes = read(path)
    length(bytes) == _integer(record["data_bytes"], "data_bytes") ||
        throw(ArgumentError("fixture payload byte count changed"))
    _hex(bytes) == record["data_sha256"] || throw(ArgumentError(
        "fixture payload checksum mismatch"))
    specifications = record["arrays"]
    specifications isa Vector || throw(ArgumentError(
        "fixture arrays must be a list"))
    length(specifications) == length(_ARRAY_NAMES) || throw(ArgumentError(
        "fixture array count changed"))
    arrays = Dict{String,Array{Float64}}()
    offset = 0
    for (expected, specification) in zip(_ARRAY_NAMES, specifications)
        name = String(specification["name"])
        name == expected || throw(ArgumentError(
            "fixture array order changed at $expected"))
        haskey(arrays, name) && throw(ArgumentError(
            "duplicate fixture array $name"))
        stored_offset = _integer(specification["offset_bytes"], name)
        stored_offset == offset || throw(ArgumentError(
            "fixture arrays must be contiguous"))
        dimensions = _dimensions(specification["dimensions"], name)
        count = _element_count(dimensions, name)
        count == _integer(specification["element_count"], name) ||
            throw(ArgumentError("fixture element count changed for $name"))
        byte_count = try
            Base.Checked.checked_mul(count, sizeof(Float64))
        catch
            throw(ArgumentError("fixture byte count overflows for $name"))
        end
        stop = try
            Base.Checked.checked_add(offset, byte_count)
        catch
            throw(ArgumentError("fixture byte range overflows for $name"))
        end
        stop <= length(bytes) || throw(ArgumentError(
            "fixture payload is truncated at $name"))
        segment = bytes[offset + 1:stop]
        _hex(segment) == specification["sha256"] || throw(ArgumentError(
            "fixture array checksum mismatch for $name"))
        values = copy(reinterpret(Float64, segment))
        all(isfinite, values) || throw(ArgumentError(
            "fixture array $name contains a nonfinite value"))
        arrays[name] = reshape(values, Tuple(dimensions))
        offset = stop
    end
    offset == length(bytes) || throw(ArgumentError(
        "fixture payload has trailing bytes"))
    arrays
end

function _atom_intervals(record, sites)
    starts = [_integer(item, "atom start") for item in record["atom_starts"]]
    stops = [_integer(item, "atom stop") for item in record["atom_stops"]]
    atoms = _integer(record["atoms"], "atoms")
    length(starts) == length(stops) == atoms || throw(ArgumentError(
        "fixture atom partition length changed"))
    intervals = [starts[index]:stops[index] for index in 1:atoms]
    first(first(intervals)) == 1 && last(last(intervals)) == sites ||
        throw(ArgumentError("fixture atom partition must cover all sites"))
    for index in eachindex(intervals)
        !isempty(intervals[index]) || throw(ArgumentError(
            "fixture atom intervals must be nonempty"))
        index == 1 || first(intervals[index]) == last(intervals[index-1]) + 1 ||
            throw(ArgumentError("fixture atom intervals must be contiguous"))
    end
    intervals
end

function _producer_matrices(bands, interaction, sites)
    h1 = zeros(Float64, sites, sites)
    bandwidth = size(bands, 1) - 1
    for offset in 0:bandwidth, first in 1:sites-offset
        second = first + offset
        h1[first, second] = bands[offset + 1, first]
        h1[second, first] = bands[offset + 1, first]
    end
    size(interaction) == (sites, sites) || throw(DimensionMismatch(
        "producer interaction dimensions changed"))
    h1, interaction, bandwidth
end

"""Load one verified, distributable H10 or H20 compact/producer fixture."""
function load_hydrogen_chain_fixture(case_name::Union{Symbol,String};
        manifest_path::AbstractString=DEFAULT_MANIFEST)
    manifest = TOML.parsefile(manifest_path)
    get(manifest, "schema", nothing) == _SCHEMA || throw(ArgumentError(
        "unknown hydrogen fixture schema"))
    get(manifest, "binary_format", nothing) == _FORMAT ||
        throw(ArgumentError("unknown hydrogen fixture binary format"))
    key = lowercase(String(case_name))
    key in ("h10", "h20") || throw(ArgumentError(
        "fixture case must be h10 or h20"))
    record = manifest["cases"][key]
    sites = _integer(record["sites"], "sites")
    arrays = _load_arrays(dirname(abspath(manifest_path)), record)
    bands = arrays["one_body_bands"]
    size(bands, 2) == sites || throw(DimensionMismatch(
        "fixture one-body site count changed"))
    one_body = HFDMRG.BandedOneBody(bands)
    interaction = HFDMRG.UnitCellInteraction(sites,
        _integer(record["phase_at_first"], "phase_at_first"),
        vec(arrays["tail_transfer"]), arrays["residual_transfer"],
        arrays["phase_source"], arrays["phase_sink"],
        vec(arrays["local_triangle"]))
    intervals = _atom_intervals(record, sites)
    h1, vee, bandwidth = _producer_matrices(bands,
        arrays["producer_vee"], sites)
    constant = Float64(record["nuclear_repulsion_ha"])
    isfinite(constant) || throw(ArgumentError(
        "fixture nuclear repulsion must be finite"))
    producer = HFDMRG.ProducerHamiltonian(h1, vee; constant,
        h1_bandwidth=bandwidth)
    (; name=Symbol(key), one_body, interaction, producer,
        atom_intervals=intervals, metadata=record,
        provenance=manifest["provenance"])
end

end
