function occupation_default_fixture(sites::Int; phase_at_first::Int=1)
    bands = zeros(2, sites)
    @inbounds for site = 1:sites
        bands[1, site] = 0.002 * (site - (sites + 1) / 2)^2
        site < sites && (bands[2, site] = -0.08)
    end
    one_body = HFDMRG.BandedOneBody(bands)
    decay, amplitude = 0.21, 0.015
    transfer = [exp(-18decay)]
    source = reshape([amplitude * exp(decay * phase) for phase = 1:18],
        18, 1)
    sink = reshape([exp(-decay * (18 + phase)) for phase = 1:18], 1, 18)
    triangle = [first == second ? 0.1 :
        amplitude * exp(-decay * (second - first))
        for second = 1:18 for first = 1:second]
    operator = HFDMRG.UnitCellInteraction(sites, phase_at_first, transfer,
        zeros(0, 0), source, sink, triangle)
    one_body, operator
end

function occupation_default_error(f)
    try
        f()
        nothing
    catch error
        error
    end
end

@testset "safe compact occupation defaults" begin
    integral_h1, integral_operator = occupation_default_fixture(36)
    rhf_control = HFDMRG.RHFStateControl(1, 1.0e-12; exact=true)
    uhf_control = HFDMRG.UHFStateControl(rhf_control)
    for forward in (true, false)
        rhf_work = HFDMRG.RHFLifecycleWorkspace(integral_h1,
            integral_operator, 1)
        rhf = HFDMRG.initialize_rhf_dimer_lifecycle(integral_h1,
            integral_operator, rhf_control, rhf_work; forward)
        @test rhf.atom_intervals == [1:18, 19:36]
        @test rhf.root.total_occupied == 1

        uhf_work = HFDMRG.UHFLifecycleWorkspace(integral_h1,
            integral_operator, 1, 1)
        uhf = HFDMRG.initialize_uhf_neel_lifecycle(integral_h1,
            integral_operator, uhf_control, uhf_work; forward)
        @test uhf.atom_intervals == [1:18, 19:36]
        @test uhf.root.alpha.total_occupied == 1
        @test uhf.root.beta.total_occupied == 1
    end

    cases = (
        (name=:partial_suffix,
            data=occupation_default_fixture(39),
            atoms=[1:18, 19:39]),
        (name=:partial_prefix_nontrivial_phase,
            data=occupation_default_fixture(39; phase_at_first=16),
            atoms=[1:21, 22:39]),
    )
    for case in cases
        one_body, operator = case.data
        for orientation in (:physical, :reflected), spin in (:rhf, :uhf)
            error = occupation_default_error() do
                HFDMRG.solve_hfdmrg(one_body, operator; spin,
                    state_maxdim=1, state_cutoff=1.0e-12,
                    energy_tolerance=1.0e-7, maximum_half_sweeps=1,
                    starting_orientation=orientation, exact=true)
            end
            @test error isa ArgumentError
            @test occursin("atom_intervals is required", sprint(showerror,
                error))
        end

        for forward in (true, false)
            rhf_work = HFDMRG.RHFLifecycleWorkspace(one_body, operator, 1)
            rhf_physical = copy(rhf_work.physical.h1)
            rhf_builder = copy(rhf_work.builder.h1)
            rhf_provenance = rhf_work.empty.provenance
            @test_throws ArgumentError HFDMRG.initialize_rhf_dimer_lifecycle(
                one_body, operator, rhf_control, rhf_work; forward)
            @test isequal(rhf_work.physical.h1, rhf_physical)
            @test isequal(rhf_work.builder.h1, rhf_builder)
            @test rhf_work.empty.provenance == rhf_provenance

            rhf = HFDMRG.initialize_rhf_dimer_lifecycle(one_body, operator,
                rhf_control, rhf_work; forward,
                atom_intervals=case.atoms)
            @test rhf.atom_intervals == case.atoms
            @test rhf.root.total_occupied == 1

            uhf_work = HFDMRG.UHFLifecycleWorkspace(one_body, operator, 1, 1)
            uhf_builder = copy(uhf_work.builder.local_interaction)
            uhf_provenance = uhf_work.empty.provenance
            @test_throws ArgumentError HFDMRG.initialize_uhf_neel_lifecycle(
                one_body, operator, uhf_control, uhf_work; forward)
            @test isequal(uhf_work.builder.local_interaction, uhf_builder)
            @test uhf_work.empty.provenance == uhf_provenance

            uhf = HFDMRG.initialize_uhf_neel_lifecycle(one_body, operator,
                uhf_control, uhf_work; forward,
                atom_intervals=case.atoms)
            @test uhf.atom_intervals == case.atoms
            @test uhf.root.alpha.total_occupied == 1
            @test uhf.root.beta.total_occupied == 1
        end
    end
end
