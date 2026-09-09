using HFDMRG
using LinearAlgebra

BLAS.set_num_threads(4)

function scaling_fixture(groups::Int)
    sites = 3 + 18(groups - 1)
    bands = zeros(4, sites)
    @inbounds for site = 1:sites
        phase = 1 + mod(15 + site - 1, 18)
        bands[1, site] = 0.0015 * (phase - 9.0)^2 - 0.04
        for offset = 1:min(3, sites-site)
            bands[offset+1, site] = -0.11 / offset +
                0.002sin(0.13site + 0.07offset)
        end
    end
    one_body = HFDMRG.BandedOneBody(bands)
    tail = [0.17]
    residual = [0.31 -0.06; 0.08 0.23]
    source = [0.012sin(0.19phase + 0.31channel) +
        0.004cos(0.07phase*channel)
        for phase = 1:18, channel = 1:3]
    sink = [0.21cos(0.11channel + 0.17phase) -
        0.03sin(0.09channel*phase)
        for channel = 1:3, phase = 1:18]
    triangle = [first == second ? 0.09 + 0.001first :
        0.018exp(-0.16(second-first)) + 0.0002(first+second)
        for second = 1:18 for first = 1:second]
    operator = HFDMRG.UnitCellInteraction(sites, 16, tail, residual,
        source, sink, triangle)
    one_body, operator, HFDMRG._traversal_intervals(operator)
end

function measure_case(groups::Int, half_sweeps::Int)
    one_body, operator, atoms = scaling_fixture(groups)
    control = HFDMRG.RHFStateControl(2, 1.0e-12)
    workspace = HFDMRG.RHFLifecycleWorkspace(one_body, operator, 2)
    lifecycle = HFDMRG.initialize_rhf_dimer_lifecycle(one_body, operator,
        control, workspace; atom_intervals=atoms)
    HFDMRG.run_fixed_half_sweep!(lifecycle, one_body, operator, control,
        workspace)
    calls_before = workspace.builder.channel_transfer_calls
    steps_before = workspace.builder.residual_transfer_steps
    visits_before = lifecycle.center_visits
    legacy_steps = 0
    GC.gc()
    measurement = @timed begin
        for _ = 1:half_sweeps
            evidence = HFDMRG.run_fixed_half_sweep!(lifecycle, one_body,
                operator, control, workspace)
            for item in evidence
                outward = item.forward ?
                    (item.cell == groups-1 ? 0 : item.cell-1) :
                    (item.cell == 1 ? 0 : groups-item.cell-1)
                legacy_steps += 1 + outward
            end
        end
    end
    visits = lifecycle.center_visits - visits_before
    calls = workspace.builder.channel_transfer_calls - calls_before
    steps = workspace.builder.residual_transfer_steps - steps_before
    (; groups, sites=operator.sites, half_sweeps, visits, calls, steps,
        legacy_steps,
        maximum_power=workspace.builder.maximum_transfer_power,
        maximum_block_rank=maximum(block.rank for block in lifecycle.collection),
        seconds=measurement.time, bytes=measurement.bytes,
        seconds_per_visit=measurement.time/visits,
        bytes_per_visit=measurement.bytes/visits)
end

println("julia_version = \"", VERSION, "\"")
println("blas_threads = ", BLAS.get_num_threads())
for groups in (10, 20)
    result = measure_case(groups, 12)
    println()
    println("[[case]]")
    for name in propertynames(result)
        value = getproperty(result, name)
        println(name, " = ", value isa AbstractString ? repr(value) : value)
    end
end
