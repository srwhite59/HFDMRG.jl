using HFDMRG
using LinearAlgebra
using Test
using TOML

const EXPECTED_EXPORTS = Set((:BandedOneBody, :HFDMRG, :HFDMRGResult,
    :RHFCompactState, :UHFDMRGResult, :UHFCompactState,
    :UnitCellInteraction, :solve_hfdmrg))
const COMPACT_COMMON_RHF_FILES = (
    "unit_cell_interaction.jl", "pairs.jl", "state.jl", "one_body.jl",
    "interaction_blocks.jl", "interaction_centers.jl",
    "fixed_state_lifecycle.jl", "nonlinear_rhf.jl",
)
const COMPACT_UHF_FILES = (
    "uhf_blocks.jl", "uhf_centers.jl", "uhf_lifecycle.jl",
    "nonlinear_uhf.jl",
)

line_count(path) = count(==(UInt8('\n')), read(path))

function compact_fixture(atoms::Int)
    sites = 18atoms
    bands = zeros(2, sites)
    @inbounds for site = 1:sites
        phase = 1 + (site - 1) % 18
        bands[1, site] = 0.001 * (phase - 9.5)^2
        site < sites && (bands[2, site] = site % 18 == 0 ? -0.035 : -0.12)
    end
    one_body = HFDMRG.BandedOneBody(bands)
    decay, amplitude = 0.24, 0.018
    transfer = [exp(-18decay)]
    source = reshape([amplitude * exp(decay * phase) for phase = 1:18], 18, 1)
    sink = reshape([exp(-decay * (18 + phase)) for phase = 1:18], 1, 18)
    triangle = [first == second ? 0.08 :
        amplitude * exp(-decay * (second - first))
        for second = 1:18 for first = 1:second]
    operator = HFDMRG.UnitCellInteraction(sites, 1, transfer, zeros(0, 0),
        source, sink, triangle)
    one_body, operator
end

function main()
    Set(names(HFDMRG)) == EXPECTED_EXPORTS || error("unexpected public exports")
    length(methods(HFDMRG.solve_hfdmrg)) == 19 ||
        error("unexpected solve_hfdmrg method count")
    isempty(Test.detect_ambiguities(HFDMRG; recursive=true)) ||
        error("public method table contains ambiguities")
    project = TOML.parsefile(joinpath(@__DIR__, "..", "..", "Project.toml"))
    project["deps"] == Dict(
        "LinearAlgebra" => "37e2e46d-f89d-539d-b4ee-838fcccc9c8e") ||
        error("unexpected runtime dependency")

    src = normpath(joinpath(@__DIR__, "..", "..", "src"))
    for file in (COMPACT_COMMON_RHF_FILES..., COMPACT_UHF_FILES...)
        source = read(joinpath(src, file), String)
        for forbidden in ("solve_hfdmrg_core", "DensityDensityBackend",
                "SlicedBasisBackend")
            occursin(forbidden, source) &&
                error("compact source reaches historical backend: $file")
        end
    end

    sites = 12
    h = Matrix(Diagonal(vcat(-2.0, -1.0, collect(1.0:sites-2))))
    v = zeros(sites, sites)
    alpha = reshape(vcat(1.0, zeros(sites-1)), sites, 1)
    beta = reshape(vcat(0.0, 1.0, zeros(sites-2)), sites, 1)
    historical = HFDMRG.solve_hfdmrg(h, v, alpha, beta; maxiter=1,
        blocksize=2, cutoff=0.0, scf_cutoff=Inf, verbose=false)
    historical isa Tuple{Matrix{Float64},Matrix{Float64},Float64} ||
        error("historical dense UHF dispatch lost tuple contract")

    one_body, operator = compact_fixture(2)
    compact_method = which(HFDMRG.solve_hfdmrg,
        (typeof(one_body), typeof(operator)))
    historical_method = which(HFDMRG.solve_hfdmrg,
        (typeof(h), typeof(v), typeof(alpha), typeof(beta)))
    compact_method == historical_method &&
        error("compact and historical UHF calls share a numerical method")
    blas_before = BLAS.get_num_threads()
    compact = HFDMRG.solve_hfdmrg(one_body, operator; spin=:uhf,
        initialization=:atomic_neel, state_maxdim=1, state_cutoff=1e-12,
        energy_tolerance=1e-7, maximum_half_sweeps=20)
    compact isa HFDMRG.UHFDMRGResult ||
        error("compact UHF dispatch lost compact result contract")
    BLAS.get_num_threads() == blas_before || error("caller BLAS changed")
    lifecycle = compact.state.lifecycle
    control = HFDMRG.UHFStateControl(HFDMRG.RHFStateControl(1, 1e-12))
    workspace = HFDMRG.UHFNonlinearWorkspace(one_body, operator, 1, 1)
    HFDMRG._preflight_uhf_nonlinear(lifecycle, one_body, operator, control,
        workspace)
    bytes = @allocated HFDMRG._preflight_uhf_nonlinear(
        lifecycle, one_body, operator, control, workspace)
    bytes == 0 || error("warmed compact UHF preflight allocated")

    common_rhf_lines = sum(line_count(joinpath(src, file))
        for file in COMPACT_COMMON_RHF_FILES)
    compact_uhf_lines = sum(line_count(joinpath(src, file))
        for file in COMPACT_UHF_FILES)
    total_lines = sum(line_count(joinpath(src, file)) for file in readdir(src)
        if endswith(file, ".jl")) +
        sum(line_count(joinpath(src, "backends", file))
            for file in readdir(joinpath(src, "backends"))
            if endswith(file, ".jl"))
    common_rhf_lines == 2799 || error("unexpected compact common/RHF line count")
    compact_uhf_lines == 1785 || error("unexpected compact UHF line count")

    println("julia_version=", VERSION)
    println("method_count=19")
    println("ambiguities=0")
    println("compact_common_rhf_lines=", common_rhf_lines)
    println("compact_uhf_lines=", compact_uhf_lines)
    println("loaded_production_lines=", total_lines)
    println("blas_threads=", blas_before)
    println("compact_uhf_preflight_allocation_bytes=", bytes)
    println("historical_uhf_tuple_contract=true")
    println("compact_uhf_result_contract=true")
    println("compact_historical_callbacks=0")
end

main()
