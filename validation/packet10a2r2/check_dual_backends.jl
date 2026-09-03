using HFDMRG
using LinearAlgebra
using Test
using TOML

const EXPECTED_EXPORTS = Set((:BandedOneBody, :HFDMRG, :HFDMRGResult,
    :RHFCompactState, :UnitCellInteraction, :solve_hfdmrg))
const COMPACT_FILES = (
    "unit_cell_interaction.jl", "pairs.jl", "state.jl", "one_body.jl",
    "interaction_blocks.jl", "interaction_centers.jl",
    "fixed_state_lifecycle.jl", "nonlinear_rhf.jl",
)

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

line_count(path) = count(==(UInt8('\n')), read(path))

function main()
    Set(names(HFDMRG)) == EXPECTED_EXPORTS || error("unexpected public exports")
    isempty(Test.detect_ambiguities(HFDMRG; recursive=true)) ||
        error("public method table contains ambiguities")
    project = TOML.parsefile(joinpath(@__DIR__, "..", "..", "Project.toml"))
    project["deps"] == Dict(
        "LinearAlgebra" => "37e2e46d-f89d-539d-b4ee-838fcccc9c8e") ||
        error("unexpected runtime dependency")
    src = normpath(joinpath(@__DIR__, "..", "..", "src"))
    for file in COMPACT_FILES
        source = read(joinpath(src, file), String)
        for forbidden in ("solve_hfdmrg_core", "DensityDensityBackend",
                "SlicedBasisBackend")
            occursin(forbidden, source) &&
                error("compact source reaches historical backend: $file")
        end
    end

    sites = 12
    h = Matrix(Diagonal(vcat(-2.0, collect(1.0:sites-1))))
    v = zeros(sites, sites)
    occupied = reshape(vcat(1.0, zeros(sites-1)), sites, 1)
    historical = HFDMRG.solve_hfdmrg(h, v, occupied; maxiter=1,
        blocksize=2, cutoff=0.0, scf_cutoff=Inf, verbose=false)
    historical isa Tuple{Matrix{Float64},Matrix{Float64},Float64} ||
        error("historical dense dispatch lost tuple contract")

    one_body, operator = compact_fixture(2)
    compact_method = which(HFDMRG.solve_hfdmrg,
        (typeof(one_body), typeof(operator)))
    historical_method = which(HFDMRG.solve_hfdmrg,
        (typeof(h), typeof(v), typeof(occupied)))
    compact_method == historical_method &&
        error("compact and historical calls share a numerical method")
    blas_before = BLAS.get_num_threads()
    compact = HFDMRG.solve_hfdmrg(one_body, operator; state_maxdim=1,
        state_cutoff=1e-12, energy_tolerance=1e-7,
        maximum_half_sweeps=20)
    compact isa HFDMRG.HFDMRGResult ||
        error("compact dispatch lost compact result contract")
    BLAS.get_num_threads() == blas_before || error("caller BLAS changed")
    lifecycle = compact.state.lifecycle
    control = HFDMRG.RHFStateControl(1, 1e-12)
    workspace = HFDMRG.RHFNonlinearWorkspace(one_body, operator, 1)
    HFDMRG._preflight_nonlinear(lifecycle, one_body, operator, control,
        workspace)
    bytes = @allocated HFDMRG._preflight_nonlinear(
        lifecycle, one_body, operator, control, workspace)
    bytes == 0 || error("warmed compact preflight allocated")

    compact_lines = sum(line_count(joinpath(src, file)) for file in COMPACT_FILES)
    total_lines = sum(line_count(joinpath(src, file)) for file in readdir(src)
        if endswith(file, ".jl")) +
        sum(line_count(joinpath(src, "backends", file))
            for file in readdir(joinpath(src, "backends"))
            if endswith(file, ".jl"))
    compact_lines == 2793 || error("unexpected compact line count")

    println("julia_version=", VERSION)
    println("method_count=", length(methods(HFDMRG.solve_hfdmrg)))
    println("ambiguities=0")
    println("compact_lines=", compact_lines)
    println("loaded_production_lines=", total_lines)
    println("blas_threads=", blas_before)
    println("compact_preflight_allocation_bytes=", bytes)
    println("historical_tuple_contract=true")
    println("compact_result_contract=true")
    println("compact_historical_callbacks=0")
end

main()
