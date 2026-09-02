function fragment_fixture(remainder::Int)
    1 <= remainder <= 17 || throw(ArgumentError("invalid fixture remainder"))
    sites = 2HFDMRG.UNIT_CELL_WIDTH + remainder
    bands = zeros(2, sites)
    @inbounds for site = 1:sites
        phase = 1 + (site - 1) % HFDMRG.UNIT_CELL_WIDTH
        bands[1, site] = 0.001 * (phase - 9.5)^2
        site < sites && (bands[2, site] = site % 18 == 0 ? -0.035 : -0.12)
    end
    one_body = HFDMRG.BandedOneBody(bands)
    decay, amplitude, onsite = 0.24, 0.018, 0.08
    transfer = [exp(-18decay)]
    source = reshape([amplitude * exp(decay * phase) for phase = 1:18],
        18, 1)
    sink = reshape([exp(-decay * (18 + phase)) for phase = 1:18], 1, 18)
    triangle = [first == second ? onsite :
        amplitude * exp(-decay * (second - first))
        for second = 1:18 for first = 1:second]
    operator = HFDMRG.UnitCellInteraction(sites, 1, transfer, zeros(0, 0),
        source, sink, triangle)
    atoms = [1:18, 19:sites]
    one_body, operator, atoms
end

function fragment_coefficients(one_body, atoms)
    reshape(HFDMRG._dimer_orbital(one_body,
        first(atoms[1]):last(atoms[2])), :, 1)
end

@testset "exact-sized terminal fragment schedules and ownership" begin
    integral_one_body, integral_operator = nucleus_fixture(2)
    integral_control = HFDMRG.RHFStateControl(1, 1.0e-12; exact=true)
    integral_workspace = HFDMRG.RHFLifecycleWorkspace(integral_one_body,
        integral_operator, 1)
    integral = HFDMRG.initialize_rhf_dimer_lifecycle(integral_one_body,
        integral_operator, integral_control, integral_workspace)
    @test integral.intervals == [1:18, 19:36]
    @test all(length(interval) == HFDMRG.UNIT_CELL_WIDTH
        for interval in integral.intervals)

    for remainder = 1:17, forward in (true, false)
        one_body, operator, atoms = fragment_fixture(remainder)
        control = HFDMRG.RHFStateControl(1, 1.0e-12; exact=true)
        workspace = HFDMRG.RHFLifecycleWorkspace(one_body, operator, 1)
        lifecycle = HFDMRG.initialize_rhf_dimer_lifecycle(one_body,
            operator, control, workspace; forward, atom_intervals=atoms)
        @test length(lifecycle.intervals) == 3
        @test length(lifecycle.intervals[end]) == remainder
        @test lifecycle.intervals == [1:18, 19:36, 37:36+remainder]
        @test all(block -> size(block.h1) == (block.rank, block.rank),
            lifecycle.collection)
        @test all(block -> size(block.pair_field) ==
            (HFDMRG.pair_dimension(block.rank),
                HFDMRG.pair_dimension(block.rank)), lifecycle.collection)
        @test lifecycle.fragment_ingestions == 1
        @test lifecycle.original_rows_read == 0
        @test lifecycle.root.rank == (forward ? 37 : 19 + remainder)
    end
end

@testset "terminal fragment numerical participation and turnaround" begin
    for remainder in (1, 2, 3, 17), forward in (true, false),
            exact in (true, false)
        one_body, operator, atoms = fragment_fixture(remainder)
        control = HFDMRG.RHFStateControl(1, 1.0e-12; exact)
        workspace = HFDMRG.RHFLifecycleWorkspace(one_body, operator, 1;
            maximum_rhs=2)
        lifecycle = HFDMRG.initialize_rhf_dimer_lifecycle(one_body,
            operator, control, workspace; forward, atom_intervals=atoms)
        coefficients = fragment_coefficients(one_body, atoms)
        projector, oracle_energy, _ = dense_rhf_state(one_body, operator,
            coefficients)
        energies = Float64[]
        for _ = 1:4
            append!(energies, getproperty.(HFDMRG.run_fixed_half_sweep!(
                lifecycle, one_body, operator, control, workspace), :energy))
        end
        exact && @test maximum(abs.(energies .- oracle_energy)) <= 8e-11
        output = zeros(operator.sites, 1)
        HFDMRG.reconstruct_terminal!(output, lifecycle, workspace)
        @test norm(transpose(output) * output - I, Inf) <= 3e-12
        exact && @test norm(output * transpose(output) - projector, Inf) <= 3e-11
        @test lifecycle.root.turnarounds == 4
        @test lifecycle.generation == 4

        endpoint = HFDMRG.initialize_rhf_dimer_lifecycle(one_body,
            operator, control, workspace; forward=false,
            atom_intervals=atoms)
        cell = endpoint.next_center
        left, right = HFDMRG._lifecycle_blocks(endpoint, workspace)
        HFDMRG._prepare_lifecycle_center!(workspace.prepared, left, right,
            cell, lifecycle.intervals[cell], lifecycle.intervals[cell + 1],
            one_body, operator)
        root = @view endpoint.root.covariance[1:endpoint.root.rank,
            1:endpoint.root.rank]
        before = HFDMRG._center_energy_fock!(workspace.prepared, root,
            endpoint.root.occupied).total
        perturbed = copy(one_body.bands)
        perturbed[1, end] += 0.02
        changed = HFDMRG.BandedOneBody(perturbed)
        HFDMRG._prepare_lifecycle_center!(workspace.prepared, left, right,
            cell, lifecycle.intervals[cell], lifecycle.intervals[cell + 1],
            changed, operator)
        after = HFDMRG._center_energy_fock!(workspace.prepared, root,
            endpoint.root.occupied).total
        @test abs(after - before) > 1e-12
    end
end

@testset "terminal fragment validation and warmed preparation" begin
    one_body, operator, atoms = fragment_fixture(3)
    control = HFDMRG.RHFStateControl(1, 1.0e-12; exact=true)
    workspace = HFDMRG.RHFLifecycleWorkspace(one_body, operator, 1)
    lifecycle = HFDMRG.initialize_rhf_dimer_lifecycle(one_body, operator,
        control, workspace; atom_intervals=atoms)
    left, right = HFDMRG._lifecycle_blocks(lifecycle, workspace)
    cell = lifecycle.next_center
    args = (workspace.prepared, left, right, cell,
        lifecycle.intervals[cell], lifecycle.intervals[cell + 1], one_body,
        operator)
    HFDMRG._prepare_lifecycle_center!(args...)
    @test @allocated(HFDMRG._prepare_lifecycle_center!(args...)) == 0
    root_before = copy(lifecycle.root.covariance)
    collection_before = copy(lifecycle.collection)
    @test_throws ArgumentError HFDMRG.initialize_rhf_dimer_lifecycle(
        one_body, operator, control, workspace;
        atom_intervals=[1:18, 20:operator.sites])
    @test isequal(root_before, lifecycle.root.covariance)
    @test collection_before == lifecycle.collection
end
