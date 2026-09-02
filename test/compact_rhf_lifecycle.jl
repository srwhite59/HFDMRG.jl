function dimer_coefficients(one_body, atoms)
    output = zeros(18atoms, atoms ÷ 2)
    for dimer = 1:atoms÷2
        rows = (2dimer-2)*18+1:2dimer*18
        output[rows, dimer] .= HFDMRG._dimer_orbital(one_body, 2dimer-1)
    end
    output
end

function dense_rhf_state(one_body, operator, coefficients)
    density = coefficients * transpose(coefficients)
    h1 = dense_h1(one_body)
    interaction = dense_interaction(operator)
    fock = h1 + Diagonal(2interaction * diag(density)) -
        interaction .* density
    energy = 2dot(h1, density) +
        2dot(diag(density), interaction * diag(density)) -
        dot(density, interaction .* density)
    density, energy, fock
end

function lifecycle_center_allocation(lifecycle, workspace, one_body, operator)
    left, right = HFDMRG._lifecycle_blocks(lifecycle, workspace)
    cell = lifecycle.next_center
    HFDMRG._prepare_lifecycle_center!(workspace.prepared, left, right, cell,
        one_body, operator)
    @allocated HFDMRG._prepare_lifecycle_center!(workspace.prepared,
        left, right, cell, one_body, operator)
end

function root_roundtrip_allocation(root, link, local_rows, workspace)
    HFDMRG._close_root!(root, link, local_rows, true, workspace; exact=true)
    HFDMRG._open_root!(root, link, true, workspace)
    @allocated begin
        HFDMRG._close_root!(root, link, local_rows, true, workspace; exact=true)
        HFDMRG._open_root!(root, link, true, workspace)
        nothing
    end
end

@testset "terminal-complete RHF fixed-state lifecycle" begin
    for atoms in (2, 10, 20), forward in (true, false)
        one_body, operator = nucleus_fixture(atoms)
        control = HFDMRG.RHFStateControl(1, 1.0e-12; exact=true)
        workspace = HFDMRG.RHFLifecycleWorkspace(one_body, operator, 1)
        lifecycle = HFDMRG.initialize_rhf_dimer_lifecycle(one_body,
            operator, control, workspace; forward=forward)
        coefficients = dimer_coefficients(one_body, atoms)
        projector, oracle_energy, _ = dense_rhf_state(one_body, operator,
            coefficients)
        all_evidence = NamedTuple[]
        for _ = 1:4
            append!(all_evidence, HFDMRG.run_fixed_half_sweep!(lifecycle,
                one_body, operator, control, workspace))
        end
        @test maximum(abs(item.energy - oracle_energy)
            for item in all_evidence) <= 5e-11
        reconstructed = zeros(size(coefficients))
        HFDMRG.reconstruct_terminal!(reconstructed, lifecycle, workspace)
        @test norm(transpose(reconstructed) * reconstructed - I, Inf) <= 2e-12
        @test norm(reconstructed * transpose(reconstructed) - projector,
            Inf) <= 2e-11
        @test lifecycle.fragment_ingestions == atoms ÷ 2
        @test lifecycle.original_rows_read == 0
        @test lifecycle.half_sweeps == 4
        @test lifecycle.generation == 4
        @test lifecycle.root.turnarounds == 4
        @test lifecycle.replacements == length(all_evidence)
        @test length(lifecycle.collection) == atoms - 1
        @test all(block -> size(block.h1) == (block.rank, block.rank),
            lifecycle.collection)
        @test lifecycle_center_allocation(lifecycle, workspace, one_body,
            operator) == 0
    end
end

@testset "lifecycle Fock, action, finite selection, and endpoints" begin
    one_body, operator = nucleus_fixture(10)
    exact = HFDMRG.RHFStateControl(1, 1.0e-12; exact=true)
    workspace = HFDMRG.RHFLifecycleWorkspace(one_body, operator, 1;
        maximum_rhs=3)
    lifecycle = HFDMRG.initialize_rhf_dimer_lifecycle(one_body, operator,
        exact, workspace)
    coefficients = dimer_coefficients(one_body, 10)
    _, oracle_energy, physical_fock = dense_rhf_state(one_body, operator,
        coefficients)
    left, right = HFDMRG._lifecycle_blocks(lifecycle, workspace)
    HFDMRG._prepare_lifecycle_center!(workspace.prepared, left, right, 1,
        one_body, operator)
    energy = HFDMRG._center_energy_fock!(workspace.prepared,
        @view(lifecycle.root.covariance[1:36, 1:36]), 1)
    @test energy.total ≈ oracle_energy atol=5e-12
    @test workspace.prepared.fock[1:36, 1:36] ≈
        physical_fock[1:36, 1:36] atol=5e-12 rtol=5e-12
    panel = reshape(collect(range(-0.4, 0.7; length=108)), 36, 3)
    output = similar(panel)
    HFDMRG.apply_center_fock!(output, workspace.prepared, panel)
    @test output ≈ physical_fock[1:36, 1:36] * panel atol=5e-12 rtol=5e-12
    vector = @view panel[:, 1]
    vector_output = @view output[:, 1]
    HFDMRG.apply_center_fock!(vector_output, workspace.prepared, vector)
    @test vector_output ≈ physical_fock[1:36, 1:36] * vector atol=5e-12 rtol=5e-12

    finite = HFDMRG.RHFStateControl(0, 1.0e-12)
    finite_workspace = HFDMRG.RHFLifecycleWorkspace(one_body, operator, 0)
    finite_lifecycle = HFDMRG.initialize_rhf_dimer_lifecycle(one_body,
        operator, finite, finite_workspace)
    first = HFDMRG.run_fixed_half_sweep!(finite_lifecycle, one_body,
        operator, finite, finite_workspace)
    @test sum(item.discarded_squared_weight for item in first) > 0.0
    for _ = 2:4
        HFDMRG.run_fixed_half_sweep!(finite_lifecycle, one_body, operator,
            finite, finite_workspace)
    end
    finite_output = zeros(180, 5)
    HFDMRG.reconstruct_terminal!(finite_output, finite_lifecycle,
        finite_workspace)
    @test norm(transpose(finite_output) * finite_output - I, Inf) <= 2e-12
    @test norm(finite_output * transpose(finite_output) -
        coefficients * transpose(coefficients), Inf) > 1e-3

    perturbed_bands = copy(one_body.bands)
    perturbed_bands[1, 1] += 0.02
    perturbed_bands[1, end] -= 0.015
    perturbed_h1 = HFDMRG.BandedOneBody(perturbed_bands)
    perturbed_coefficients = dimer_coefficients(perturbed_h1, 10)
    _, perturbed_energy, _ = dense_rhf_state(perturbed_h1, operator,
        perturbed_coefficients)
    @test abs(perturbed_energy - oracle_energy) > 1e-5
    for orientation in (true, false)
        local_workspace = HFDMRG.RHFLifecycleWorkspace(perturbed_h1,
            operator, 1)
        local_lifecycle = HFDMRG.initialize_rhf_dimer_lifecycle(perturbed_h1,
            operator, exact, local_workspace; forward=orientation)
        evidence = HFDMRG.run_fixed_half_sweep!(local_lifecycle,
            perturbed_h1, operator, exact, local_workspace)
        @test maximum(abs(item.energy - perturbed_energy)
            for item in evidence) <= 5e-11
    end
end

@testset "lifecycle atomicity, hostile capacity, and allocation" begin
    one_body, operator = nucleus_fixture(10)
    control = HFDMRG.RHFStateControl(1, 1.0e-12; exact=true)
    workspace = HFDMRG.RHFLifecycleWorkspace(one_body, operator, 1)
    lifecycle = HFDMRG.initialize_rhf_dimer_lifecycle(one_body, operator,
        control, workspace)
    root_before = copy(lifecycle.root.covariance)
    collection_before = copy(lifecycle.collection)
    wrong = HFDMRG.RHFStateControl(1, 1.0e-10; exact=true)
    @test_throws ArgumentError HFDMRG.run_fixed_half_sweep!(lifecycle,
        one_body, operator, wrong, workspace)
    @test isequal(root_before, lifecycle.root.covariance)
    @test collection_before == lifecycle.collection

    malformed = lifecycle.collection[2]
    old_provenance = malformed.provenance
    malformed.provenance = xor(old_provenance, UInt(1))
    root_rank_before = lifecycle.root.rank
    next_center_before = lifecycle.next_center
    replacements_before = lifecycle.replacements
    @test_throws ArgumentError HFDMRG.run_fixed_half_sweep!(lifecycle,
        one_body, operator, control, workspace)
    @test lifecycle.root.rank == root_rank_before
    @test lifecycle.next_center == next_center_before
    @test lifecycle.replacements == replacements_before
    malformed.provenance = old_provenance

    original_link = malformed.link
    malformed.link = HFDMRG.RHFStateLink(original_link.old_rank - 1,
        original_link.physical_width + 1, original_link.close_map,
        original_link.completed_map, original_link.selection)
    reconstruction = fill(7.25, 180, 5)
    @test_throws DimensionMismatch HFDMRG.reconstruct_terminal!(
        reconstruction, lifecycle, workspace)
    @test all(==(7.25), reconstruction)
    malformed.link = original_link

    left, right = HFDMRG._lifecycle_blocks(lifecycle, workspace)
    local_rows = left.rank + 18
    link = HFDMRG._factor_state_link(@view(lifecycle.root.covariance[
        1:local_rows, 1:local_rows]), left.rank, 18, control)
    inactive_first = reinterpret(UInt64, copy(workspace.root.first[
        lifecycle.root.rank+1:end, :]))
    allocation = root_roundtrip_allocation(lifecycle.root, link, local_rows,
        workspace.root)
    @test allocation == 0
    @test reinterpret(UInt64, copy(workspace.root.first[
        lifecycle.root.rank+1:end, :])) == inactive_first
    @test BLAS.get_num_threads() >= 1

    tiny = 1.0e-24
    exact_link = HFDMRG._factor_state_link(reshape([tiny], 1, 1), 1, 0,
        control; cross=reshape([sqrt(tiny * (1.0-tiny))], 1, 1))
    @test exact_link.active_rank == 1
    cluster = HFDMRG._factor_state_link(Diagonal([0.5, 0.5]), 2, 0,
        HFDMRG.RHFStateControl(1, 0.0))
    @test cluster.active_rank == 0
    @test cluster.selection.limiting_reason == :maximum_rank
    aggregate = HFDMRG._factor_state_link(Diagonal([0.8, 0.9]), 2, 0,
        HFDMRG.RHFStateControl(2, 0.2))
    @test aggregate.active_rank == 1
    @test aggregate.selection.completed == 1
    @test aggregate.selection.discarded_squared_weight ≈ 0.19
    endpoint_weight = 4.0e-13
    endpoint = HFDMRG._factor_state_link(
        Diagonal([1.0-endpoint_weight]), 1, 0,
        HFDMRG.RHFStateControl(1, 1.0e-12))
    @test endpoint.active_rank == 0
    @test endpoint.selection.completed == 1
    @test endpoint.selection.discarded_squared_weight ≈
        1.0 - (1.0-endpoint_weight)^2 rtol=1.0e-4
    endpoint_retained = HFDMRG._factor_state_link(
        Diagonal([1.0-endpoint_weight]), 1, 0,
        HFDMRG.RHFStateControl(1, 1.0e-14))
    @test endpoint_retained.active_rank == 1
    @test endpoint_retained.selection.discarded_squared_weight == 0.0
    limited = HFDMRG._factor_state_link(Diagonal([0.8, 0.9]), 2, 0,
        HFDMRG.RHFStateControl(1, 0.0))
    @test limited.active_rank == 1
    @test limited.selection.limiting_reason == :maximum_rank
end
