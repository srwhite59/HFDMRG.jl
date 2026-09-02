using HFDMRG
using LinearAlgebra
using TOML

const EXPECTED_EXPORTS = Set((:BandedOneBody, :HFDMRG, :HFDMRGResult,
    :RHFCompactState, :UnitCellInteraction, :solve_hfdmrg))

function fixture(atoms::Int)
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
    names(HFDMRG) |> Set == EXPECTED_EXPORTS || error("unexpected public exports")
    project = TOML.parsefile(joinpath(@__DIR__, "..", "..", "Project.toml"))
    project["deps"] == Dict(
        "LinearAlgebra" => "37e2e46d-f89d-539d-b4ee-838fcccc9c8e") ||
        error("unexpected runtime dependency")

    one_body, operator = fixture(2)
    blas_before = BLAS.get_num_threads()
    result = HFDMRG.solve_hfdmrg(one_body, operator; state_maxdim=1,
        state_cutoff=1.0e-12, energy_tolerance=1.0e-7,
        maximum_half_sweeps=20)
    BLAS.get_num_threads() == blas_before || error("caller BLAS setting changed")
    result.converged || error("compact H2 solve did not converge")
    lifecycle = result.state.lifecycle
    control = HFDMRG.RHFStateControl(1, 1.0e-12)
    workspace = HFDMRG.RHFNonlinearWorkspace(one_body, operator, 1)
    HFDMRG._preflight_nonlinear(lifecycle, one_body, operator, control,
        workspace)
    preflight_bytes = @allocated HFDMRG._preflight_nonlinear(
        lifecycle, one_body, operator, control, workspace)
    preflight_bytes == 0 || error("warmed nonlinear preflight allocated")

    println("julia_version=", VERSION)
    println("method_count=", length(methods(HFDMRG.solve_hfdmrg)))
    println("export_count=", length(EXPECTED_EXPORTS))
    println("runtime_dependencies=1")
    println("blas_threads=", blas_before)
    println("preflight_allocation_bytes=", preflight_bytes)
    println("h2_energy=", result.represented_energy)
    println("h2_half_sweeps=", result.half_sweeps)
end

main()
