function dense_uhf_oracle(one_body, operator, basis_alpha, basis_beta,
        density_alpha, density_beta)
    h1 = dense_h1(one_body); interaction = dense_interaction(operator)
    physical_alpha = basis_alpha * density_alpha * transpose(basis_alpha)
    physical_beta = basis_beta * density_beta * transpose(basis_beta)
    total_diagonal = diag(physical_alpha + physical_beta)
    fock_alpha = transpose(basis_alpha) *
        (h1 + Diagonal(interaction * total_diagonal) -
            interaction .* physical_alpha) * basis_alpha
    fock_beta = transpose(basis_beta) *
        (h1 + Diagonal(interaction * total_diagonal) -
            interaction .* physical_beta) * basis_beta
    energy = dot(h1, physical_alpha) + dot(h1, physical_beta) +
        0.5dot(total_diagonal, interaction * total_diagonal) -
        0.5dot(physical_alpha, interaction .* physical_alpha) -
        0.5dot(physical_beta, interaction .* physical_beta)
    energy, fock_alpha, fock_beta, physical_alpha, physical_beta
end

function unequal_uhf_center_fixture()
    one_body, operator = nucleus_fixture(4); provenance = UInt(0x9a3a)
    first = HFDMRG.seed_uhf_interval_block(one_body, operator, 1, 1:18;
        provenance)
    second = HFDMRG.seed_uhf_interval_block(one_body, operator, 2, 19:36;
        provenance)
    alpha_map = deterministic_isometry(36, 3, 4)
    beta_map = deterministic_isometry(36, 2, 8)
    alpha_link = HFDMRG.RHFStateLink(18, 18, alpha_map,
        HFDMRG.StateSelection(0, 3, 33, 0.0, :exact))
    beta_link = HFDMRG.RHFStateLink(18, 18, beta_map,
        HFDMRG.StateSelection(0, 2, 34, 0.0, :exact))
    builder = HFDMRG.UHFBlockBuildWorkspace(40, 40,
        HFDMRG._maximum_channel_rank(operator))
    left = HFDMRG.build_uhf_outer_block(first, second, alpha_link, beta_link,
        one_body, operator, builder)
    right = HFDMRG._empty_uhf_block(operator, provenance)
    prepared = HFDMRG.PreparedUHFCenter(39, 38,
        HFDMRG._maximum_channel_rank(operator); maximum_rhs=3)
    HFDMRG.prepare_uhf_lifecycle_center!(prepared, left, right, 3, 37:54,
        55:72, one_body, operator)
    basis_alpha = zeros(72, 39); basis_beta = zeros(72, 38)
    basis_alpha[1:36, 1:3] .= alpha_map
    basis_alpha[37:72, 4:39] .= Matrix{Float64}(I, 36, 36)
    basis_beta[1:36, 1:2] .= beta_map
    basis_beta[37:72, 3:38] .= Matrix{Float64}(I, 36, 36)
    one_body, operator, left, right, prepared, basis_alpha, basis_beta
end

@testset "independent-spin UHF center arithmetic" begin
    unbalanced_alpha = [1.0 0.0; 0.0 1e-14; 0.0 0.0]
    unbalanced_beta = [1.0 0.0; 0.0 1e14; 0.0 0.0]
    unbalanced = HFDMRG._compress_cross_hartree(copy(unbalanced_alpha),
        copy(unbalanced_beta), ones(2), 2, 2)
    @test unbalanced.alpha * Diagonal(unbalanced.weight) *
        transpose(unbalanced.beta) ≈ unbalanced_alpha *
        transpose(unbalanced_beta) atol=3e-14 rtol=3e-14

    one_body, operator, left, right, prepared, basis_alpha, basis_beta =
        unequal_uhf_center_fixture()
    alpha_orbitals = deterministic_isometry(39, 2, 5)
    beta_orbitals = deterministic_isometry(38, 2, 7)
    alpha_density = alpha_orbitals * transpose(alpha_orbitals)
    beta_density = beta_orbitals * transpose(beta_orbitals)
    expected_energy, expected_alpha, expected_beta, _, _ = dense_uhf_oracle(
        one_body, operator, basis_alpha, basis_beta, alpha_density,
        beta_density)
    energy = HFDMRG.uhf_center_energy_fock!(prepared, alpha_density,
        beta_density, 2, 2)
    @test energy.total ≈ expected_energy atol=3e-11
    interaction = dense_interaction(operator)
    physical_alpha = basis_alpha * alpha_density * transpose(basis_alpha)
    physical_beta = basis_beta * beta_density * transpose(basis_beta)
    total_diagonal = diag(physical_alpha + physical_beta)
    expected_hartree = 0.5dot(total_diagonal,
        interaction * total_diagonal)
    expected_exchange = -0.5dot(physical_alpha,
        interaction .* physical_alpha) - 0.5dot(physical_beta,
        interaction .* physical_beta)
    @test energy.hartree ≈ expected_hartree atol=3e-11
    @test energy.exchange ≈ expected_exchange atol=3e-11
    @test prepared.alpha.fock[1:39, 1:39] ≈ expected_alpha atol=3e-11
    @test prepared.beta.fock[1:38, 1:38] ≈ expected_beta atol=3e-11

    alpha_vector = collect(range(-0.4, 0.6; length=39))
    beta_vector = collect(range(0.7, -0.2; length=38))
    alpha_output = similar(alpha_vector); beta_output = similar(beta_vector)
    HFDMRG.apply_uhf_center_fock!(alpha_output, beta_output, prepared,
        alpha_vector, beta_vector)
    @test alpha_output ≈ expected_alpha * alpha_vector atol=3e-11
    @test beta_output ≈ expected_beta * beta_vector atol=3e-11
    alpha_panel = reshape(collect(range(-0.3, 0.8; length=117)), 39, 3)
    beta_panel = reshape(collect(range(0.5, -0.7; length=114)), 38, 3)
    alpha_panel_output = similar(alpha_panel); beta_panel_output = similar(beta_panel)
    HFDMRG.apply_uhf_center_fock!(alpha_panel_output, beta_panel_output,
        prepared, alpha_panel, beta_panel)
    @test alpha_panel_output ≈ expected_alpha * alpha_panel atol=3e-11
    @test beta_panel_output ≈ expected_beta * beta_panel atol=3e-11

    @test left.alpha.rank == 3 && left.beta.rank == 2
    @test left.alpha.left_core === left.beta.left_core
    @test HFDMRG.cross_hartree_rank(left.cross_hartree) <= 6
    @test size(left.cross_hartree.alpha, 2) ==
        size(left.cross_hartree.beta, 2) ==
        HFDMRG.cross_hartree_rank(left.cross_hartree)
    swapped = HFDMRG.spin_swap_uhf_block(left, operator)
    @test swapped.alpha.rank == left.beta.rank
    @test swapped.beta.rank == left.alpha.rank
    @test swapped.cross_hartree.alpha ≈ left.cross_hartree.beta
    @test swapped.cross_hartree.beta ≈ left.cross_hartree.alpha
    reflected = HFDMRG.reflect_uhf_block(left, 4, 72, operator)
    restored = HFDMRG.reflect_uhf_block(reflected, 4, 72, operator)
    @test restored.alpha.h1 == left.alpha.h1
    @test restored.beta.pair_field == left.beta.pair_field
    @test restored.cross_hartree.alpha == left.cross_hartree.alpha

    alpha_snapshot = copy(prepared.alpha.h1)
    beta_snapshot = copy(prepared.beta.h1)
    small = HFDMRG.PreparedUHFCenter(38, 38, 1)
    @test_throws DimensionMismatch HFDMRG.prepare_uhf_lifecycle_center!(small,
        left, right, 3, 37:54, 55:72, one_body, operator)
    @test all(isnan, small.alpha.h1)
    @test all(isnan, small.beta.h1)
    @test_throws ArgumentError HFDMRG.apply_uhf_center_fock!(alpha_vector,
        beta_output, prepared, alpha_vector, beta_vector)
    @test prepared.alpha.h1 == alpha_snapshot
    @test prepared.beta.h1 == beta_snapshot

    # The direct dense fixture above is the controlling RHF-limit oracle; the
    # equality of the two spin Focks is also checked on the physical H2 center.
    h2, v2 = nucleus_fixture(2); empty = HFDMRG._empty_uhf_block(v2, UInt(3))
    rhf_center = HFDMRG.PreparedUHFCenter(36, 36, 1)
    HFDMRG.prepare_uhf_lifecycle_center!(rhf_center, empty, empty, 1, 1:18,
        19:36, h2, v2)
    orbital = deterministic_isometry(36, 1, 11)
    density = orbital * transpose(orbital)
    rhf_result = HFDMRG.uhf_center_energy_fock!(rhf_center, density, density,
        1, 1)
    dense_rhf, dense_fock, _, _, _ = dense_uhf_oracle(h2, v2,
        Matrix{Float64}(I, 36, 36), Matrix{Float64}(I, 36, 36), density,
        density)
    @test rhf_result.total ≈ dense_rhf atol=3e-11
    @test rhf_center.alpha.fock[1:36, 1:36] ≈ dense_fock atol=3e-11
    @test rhf_center.alpha.fock[1:36, 1:36] ≈
        rhf_center.beta.fock[1:36, 1:36] atol=3e-11
end

function finite_fragment_fixture()
    one_body_full, operator_full = nucleus_fixture(3)
    sites = 37
    one_body = HFDMRG.BandedOneBody(one_body_full.bands[:, 1:sites])
    operator = HFDMRG.UnitCellInteraction(sites, 1,
        operator_full.tail_transfer, operator_full.residual_transfer,
        operator_full.phase_source, operator_full.phase_sink,
        operator_full.local_triangle)
    atoms = UnitRange{Int}[1:18, 19:36, 37:37]
    one_body, operator, atoms
end

function check_uhf_lifecycle(atoms::Int, exact::Bool, forward::Bool,
        spin_swapped::Bool)
    one_body, operator = nucleus_fixture(atoms)
    atom_rows = [18(index-1)+1:18index for index = 1:atoms]
    maximum = max(atoms ÷ 2 + 1, 4)
    state = HFDMRG.RHFStateControl(maximum, exact ? 0.0 : 1e-12; exact)
    control = HFDMRG.UHFStateControl(state)
    workspace = HFDMRG.UHFLifecycleWorkspace(one_body, operator, maximum,
        maximum)
    lifecycle = HFDMRG.initialize_uhf_neel_lifecycle(one_body, operator,
        control, workspace; atom_intervals=atom_rows, forward, spin_swapped)
    evidence = [HFDMRG.run_uhf_fixed_half_sweep!(lifecycle, one_body,
        operator, control, workspace) for _ = 1:4]
    nalpha = spin_swapped ? atoms ÷ 2 : (atoms + 1) ÷ 2
    nbeta = atoms - nalpha
    alpha = zeros(18atoms, nalpha); beta = zeros(18atoms, nbeta)
    HFDMRG.reconstruct_uhf_terminal!(alpha, beta, lifecycle, workspace)
    Da = alpha * transpose(alpha); Db = beta * transpose(beta)
    expected, _, _, _, _ = dense_uhf_oracle(one_body, operator,
        Matrix{Float64}(I, 18atoms, 18atoms),
        Matrix{Float64}(I, 18atoms, 18atoms), Da, Db)
    left, right = HFDMRG._uhf_lifecycle_blocks(lifecycle, workspace)
    cell = lifecycle.next_center
    HFDMRG.prepare_uhf_lifecycle_center!(workspace.prepared, left, right,
        cell, lifecycle.intervals[cell], lifecycle.intervals[cell+1],
        one_body, operator)
    represented = HFDMRG.uhf_center_energy_fock!(workspace.prepared,
        @view(lifecycle.root.alpha.covariance[1:lifecycle.root.alpha.rank,
            1:lifecycle.root.alpha.rank]),
        @view(lifecycle.root.beta.covariance[1:lifecycle.root.beta.rank,
            1:lifecycle.root.beta.rank]), lifecycle.root.alpha.occupied,
        lifecycle.root.beta.occupied)
    @test represented.total ≈ expected atol=5e-10
    @test transpose(alpha) * alpha ≈ I atol=3e-11
    @test transpose(beta) * beta ≈ I atol=3e-11
    @test lifecycle.fragment_ingestions == atoms
    @test lifecycle.original_rows_read == 0
    @test lifecycle.half_sweeps == 4
    @test lifecycle.root.alpha.turnarounds == 4
    @test lifecycle.root.beta.turnarounds == 4
    @test all(length(half) > 0 for half in evidence)
    @test all(block.alpha.left_core === block.beta.left_core
        for block in lifecycle.collection)
    @test all(isapprox(transpose(block.alpha.link.close_map) *
            block.alpha.link.close_map, I; atol=3e-11) &&
        isapprox(transpose(block.beta.link.close_map) *
            block.beta.link.close_map, I; atol=3e-11)
        for block in lifecycle.collection)
    lifecycle, workspace, represented.total, alpha, beta
end

function cross_factor_from_pair_field(field, rank)
    rank == 0 && return HFDMRG.UHFCrossHartree(zeros(0, 0), zeros(0, 0),
        Float64[], 0, 0)
    factors = eigen(Symmetric(field))
    scale = maximum(abs, factors.values; init=0.0)
    tolerance = 128eps(Float64) * max(length(factors.values), 1) *
        max(scale, 1.0)
    selected = findall(value -> abs(value) > tolerance, factors.values)
    basis = Matrix(@view factors.vectors[:, selected])
    HFDMRG.UHFCrossHartree(basis, copy(basis),
        Vector(@view(factors.values[selected])), rank, rank)
end

function uhf_rhf_limit_lifecycle(rhf, control, operator)
    provenance = hash((operator.provenance, rhf.intervals,
        rhf.atom_intervals, control.maximum_active, control.cutoff,
        control.exact, control.maximum_active, control.cutoff, control.exact,
        false))
    collection = HFDMRG.UHFOuterBlock[]
    for source in rhf.collection
        alpha = deepcopy(source); beta = deepcopy(source)
        alpha.provenance = provenance; beta.provenance = provenance
        beta.left_core = alpha.left_core; beta.right_core = alpha.right_core
        push!(collection, HFDMRG.UHFOuterBlock(alpha, beta,
            cross_factor_from_pair_field(source.pair_field, source.rank),
            provenance))
    end
    alpha_root = deepcopy(rhf.root); beta_root = deepcopy(rhf.root)
    HFDMRG.UHFFixedLifecycle(collection,
        HFDMRG.UHFMovingRoot(alpha_root, beta_root), copy(rhf.intervals),
        copy(rhf.atom_intervals), rhf.cells, rhf.forward, rhf.next_center,
        provenance, rhf.generation, rhf.half_sweeps,
        2rhf.fragment_ingestions, 0, 0, 0)
end

@testset "independent-spin fixed-state lifecycle" begin
    for atoms in (2, 10, 20), exact in (false, true), forward in (false, true)
        lifecycle, workspace, energy, alpha, beta = check_uhf_lifecycle(atoms,
            exact, forward, false)
        @test isfinite(energy)
        @test length(lifecycle.collection) == atoms - 1
        @test HFDMRG.uhf_lifecycle_storage_bytes(lifecycle,
            nucleus_fixture(atoms)[2], workspace).collection > 0
    end
    physical = check_uhf_lifecycle(10, false, true, false)
    swapped = check_uhf_lifecycle(10, false, true, true)
    @test physical[3] ≈ swapped[3] atol=5e-11
    @test physical[4] * transpose(physical[4]) ≈
        swapped[5] * transpose(swapped[5]) atol=3e-10
    @test physical[5] * transpose(physical[5]) ≈
        swapped[4] * transpose(swapped[4]) atol=3e-10
end


@testset "H10/H20 RHF-limit fixed-state parity" begin
    for atoms in (10, 20), exact in (false, true)
        one_body, operator = nucleus_fixture(atoms)
        atom_rows = [18(index-1)+1:18index for index = 1:atoms]
        maximum = atoms ÷ 2 + 1
        control = HFDMRG.RHFStateControl(maximum, exact ? 0.0 : 1e-12;
            exact)
        rhf_workspace = HFDMRG.RHFLifecycleWorkspace(one_body, operator,
            maximum)
        rhf = HFDMRG.initialize_rhf_dimer_lifecycle(one_body, operator,
            control, rhf_workspace; atom_intervals=atom_rows)
        uhf = uhf_rhf_limit_lifecycle(rhf, control, operator)
        uhf_workspace = HFDMRG.UHFLifecycleWorkspace(one_body, operator,
            maximum, maximum)
        uhf_control = HFDMRG.UHFStateControl(control)
        for _ = 1:4
            HFDMRG.run_fixed_half_sweep!(rhf, one_body, operator, control,
                rhf_workspace)
            HFDMRG.run_uhf_fixed_half_sweep!(uhf, one_body, operator,
                uhf_control, uhf_workspace)
        end
        restricted = zeros(18atoms, atoms ÷ 2)
        alpha = similar(restricted); beta = similar(restricted)
        HFDMRG.reconstruct_terminal!(restricted, rhf, rhf_workspace)
        HFDMRG.reconstruct_uhf_terminal!(alpha, beta, uhf, uhf_workspace)
        restricted_projector = restricted * transpose(restricted)
        @test alpha * transpose(alpha) ≈ restricted_projector atol=3e-10
        @test beta * transpose(beta) ≈ restricted_projector atol=3e-10
        left, right = HFDMRG._uhf_lifecycle_blocks(uhf, uhf_workspace)
        cell = uhf.next_center
        HFDMRG.prepare_uhf_lifecycle_center!(uhf_workspace.prepared, left,
            right, cell, uhf.intervals[cell], uhf.intervals[cell+1],
            one_body, operator)
        uhf_energy = HFDMRG.uhf_center_energy_fock!(uhf_workspace.prepared,
            @view(uhf.root.alpha.covariance[1:uhf.root.alpha.rank,
                1:uhf.root.alpha.rank]),
            @view(uhf.root.beta.covariance[1:uhf.root.beta.rank,
                1:uhf.root.beta.rank]), uhf.root.alpha.occupied,
            uhf.root.beta.occupied)
        rhf_left, rhf_right = HFDMRG._lifecycle_blocks(rhf, rhf_workspace)
        HFDMRG._prepare_lifecycle_center!(rhf_workspace.prepared, rhf_left,
            rhf_right, rhf.next_center, rhf.intervals[rhf.next_center],
            rhf.intervals[rhf.next_center+1], one_body, operator)
        rhf_energy = HFDMRG._center_energy_fock!(rhf_workspace.prepared,
            @view(rhf.root.covariance[1:rhf.root.rank, 1:rhf.root.rank]),
            rhf.root.occupied)
        @test uhf_energy.total ≈ rhf_energy.total atol=5e-10
    end
end

@testset "asymmetric fragment and zero spin sector" begin
    one_body, operator, atoms = finite_fragment_fixture()
    control = HFDMRG.UHFStateControl(HFDMRG.RHFStateControl(4, 1e-12))
    for forward in (false, true), swap in (false, true)
        workspace = HFDMRG.UHFLifecycleWorkspace(one_body, operator, 4, 4)
        lifecycle = HFDMRG.initialize_uhf_neel_lifecycle(one_body, operator,
            control, workspace; atom_intervals=atoms, forward,
            spin_swapped=swap)
        for _ = 1:4
            HFDMRG.run_uhf_fixed_half_sweep!(lifecycle, one_body, operator,
                control, workspace)
        end
        nalpha = swap ? 1 : 2; nbeta = swap ? 2 : 1
        alpha = zeros(37, nalpha); beta = zeros(37, nbeta)
        HFDMRG.reconstruct_uhf_terminal!(alpha, beta, lifecycle, workspace)
        @test transpose(alpha) * alpha ≈ I atol=3e-11
        @test transpose(beta) * beta ≈ I atol=3e-11
        @test lifecycle.intervals[end] == 37:37
    end
    empty = HFDMRG._empty_uhf_block(operator, UInt(9))
    @test empty.alpha.rank == empty.beta.rank == 0
    @test HFDMRG.cross_hartree_rank(empty.cross_hartree) == 0

    h2, v2 = nucleus_fixture(2); provenance = UInt(12)
    first = HFDMRG.seed_uhf_interval_block(h2, v2, 1, 1:18; provenance)
    second = HFDMRG.seed_uhf_interval_block(h2, v2, 2, 19:36; provenance)
    alpha_completed = deterministic_isometry(36, 1, 17)
    alpha_link = HFDMRG.RHFStateLink(18, 18, zeros(36, 0),
        alpha_completed, HFDMRG.StateSelection(1, 0, 35, 0.0, :exact))
    beta_map = deterministic_isometry(36, 2, 19)
    beta_link = HFDMRG.RHFStateLink(18, 18, beta_map,
        HFDMRG.StateSelection(0, 2, 34, 0.0, :exact))
    builder = HFDMRG.UHFBlockBuildWorkspace(36, 36, 1)
    zero_alpha = HFDMRG.build_uhf_outer_block(first, second, alpha_link,
        beta_link, h2, v2, builder)
    @test zero_alpha.alpha.rank == 0
    @test zero_alpha.beta.rank == 2
    @test HFDMRG.cross_hartree_rank(zero_alpha.cross_hartree) == 0
    @test zero_alpha.alpha.completed_count == 1

    hostile = HFDMRG.PreparedUHFCenter(40, 41, 1)
    alpha_inactive = reinterpret(UInt64, copy(hostile.alpha.h1[37:end, :]))
    beta_inactive = reinterpret(UInt64, copy(hostile.beta.h1[37:end, :]))
    cross_alpha_inactive = reinterpret(UInt64,
        copy(hostile.cross_alpha[HFDMRG.pair_dimension(36)+1:end, :]))
    cross_beta_inactive = reinterpret(UInt64,
        copy(hostile.cross_beta[HFDMRG.pair_dimension(36)+1:end, :]))
    empty_h2 = HFDMRG._empty_uhf_block(v2, UInt(13))
    HFDMRG.prepare_uhf_lifecycle_center!(hostile, empty_h2, empty_h2, 1,
        1:18, 19:36, h2, v2)
    @test reinterpret(UInt64, copy(hostile.alpha.h1[37:end, :])) ==
        alpha_inactive
    @test reinterpret(UInt64, copy(hostile.beta.h1[37:end, :])) ==
        beta_inactive
    @test reinterpret(UInt64,
        copy(hostile.cross_alpha[HFDMRG.pair_dimension(36)+1:end, :])) ==
        cross_alpha_inactive
    @test reinterpret(UInt64,
        copy(hostile.cross_beta[HFDMRG.pair_dimension(36)+1:end, :])) ==
        cross_beta_inactive

    atomic_workspace = HFDMRG.UHFLifecycleWorkspace(h2, v2, 2, 2)
    atomic_control = HFDMRG.UHFStateControl(HFDMRG.RHFStateControl(2, 1e-12))
    atomic = HFDMRG.initialize_uhf_neel_lifecycle(h2, v2, atomic_control,
        atomic_workspace; atom_intervals=[1:18, 19:36])
    alpha_before = copy(atomic.root.alpha.covariance)
    beta_before = copy(atomic.root.beta.covariance)
    collection_before = copy(atomic.collection)
    bad_control = HFDMRG.UHFStateControl(HFDMRG.RHFStateControl(1, 1e-12))
    @test_throws ArgumentError HFDMRG.run_uhf_fixed_half_sweep!(atomic, h2,
        v2, bad_control, atomic_workspace)
    @test isequal(atomic.root.alpha.covariance, alpha_before)
    @test isequal(atomic.root.beta.covariance, beta_before)
    @test all(atomic.collection[index] === collection_before[index]
        for index in eachindex(collection_before))

    midflight_workspace = HFDMRG.UHFLifecycleWorkspace(h2, v2, 2, 2)
    midflight = HFDMRG.initialize_uhf_neel_lifecycle(h2, v2, atomic_control,
        midflight_workspace; atom_intervals=[1:18, 19:36])
    alpha_before = copy(midflight.root.alpha.covariance)
    beta_before = copy(midflight.root.beta.covariance)
    alpha_state = (midflight.root.alpha.rank, midflight.root.alpha.occupied,
        midflight.root.alpha.closes, midflight.root.alpha.opens,
        midflight.root.alpha.turnarounds)
    beta_state = (midflight.root.beta.rank, midflight.root.beta.occupied,
        midflight.root.beta.closes, midflight.root.beta.opens,
        midflight.root.beta.turnarounds)
    collection_before = copy(midflight.collection)
    midflight_workspace.beta_root = HFDMRG.RHFRootWorkspace(zeros(0, 0),
        midflight_workspace.beta_root.second)
    @test_throws BoundsError HFDMRG.run_uhf_fixed_half_sweep!(midflight, h2,
        v2, atomic_control, midflight_workspace)
    @test isequal(midflight.root.alpha.covariance, alpha_before)
    @test isequal(midflight.root.beta.covariance, beta_before)
    @test (midflight.root.alpha.rank, midflight.root.alpha.occupied,
        midflight.root.alpha.closes, midflight.root.alpha.opens,
        midflight.root.alpha.turnarounds) == alpha_state
    @test (midflight.root.beta.rank, midflight.root.beta.occupied,
        midflight.root.beta.closes, midflight.root.beta.opens,
        midflight.root.beta.turnarounds) == beta_state
    @test all(midflight.collection[index] === collection_before[index]
        for index in eachindex(collection_before))
end

function uhf_nonfactor_allocations(prepared, left, right, one_body, operator,
        alpha, beta, alpha_output, beta_output, alpha_input, beta_input)
    preparation = @allocated HFDMRG.prepare_uhf_lifecycle_center!(prepared,
        left, right, 1, 1:18, 19:36, one_body, operator)
    fock = @allocated HFDMRG.uhf_center_energy_fock!(prepared, alpha, beta,
        1, 1)
    action = @allocated HFDMRG.apply_uhf_center_fock!(alpha_output,
        beta_output, prepared, alpha_input, beta_input)
    preparation, fock, action
end

function uhf_root_allocations(alpha_root, beta_root, link, alpha_workspace,
        beta_workspace)
    alpha_close = @allocated HFDMRG._close_root!(alpha_root, link, 2, true,
        alpha_workspace; exact=true)
    alpha_open = @allocated HFDMRG._open_root!(alpha_root, link, true,
        alpha_workspace)
    beta_close = @allocated HFDMRG._close_root!(beta_root, link, 2, true,
        beta_workspace; exact=true)
    beta_open = @allocated HFDMRG._open_root!(beta_root, link, true,
        beta_workspace)
    alpha_close, alpha_open, beta_close, beta_open
end

@testset "UHF warmed nonfactor allocation and BLAS" begin
    threads = BLAS.get_num_threads()
    one_body, operator = nucleus_fixture(2); provenance = UInt(4)
    left = HFDMRG._empty_uhf_block(operator, provenance)
    right = HFDMRG._empty_uhf_block(operator, provenance)
    prepared = HFDMRG.PreparedUHFCenter(36, 36, 1; maximum_rhs=2)
    alpha = zeros(36, 36); beta = zeros(36, 36)
    alpha[1, 1] = 1.0; beta[2, 2] = 1.0
    HFDMRG.prepare_uhf_lifecycle_center!(prepared, left, right, 1, 1:18,
        19:36, one_body, operator)
    HFDMRG.uhf_center_energy_fock!(prepared, alpha, beta, 1, 1)
    alpha_input = ones(36, 2); beta_input = fill(0.5, 36, 2)
    alpha_output = similar(alpha_input); beta_output = similar(beta_input)
    HFDMRG.apply_uhf_center_fock!(alpha_output, beta_output, prepared,
        alpha_input, beta_input)
    @test uhf_nonfactor_allocations(prepared, left, right, one_body, operator,
        alpha, beta, alpha_output, beta_output, alpha_input, beta_input) ==
        (0, 0, 0)
    root_covariance = fill(NaN, 4, 4)
    root_covariance[1:2, 1:2] .= [1.0 0.0; 0.0 0.0]
    alpha_root = HFDMRG.RHFMovingRoot(copy(root_covariance), 2, 1, 1,
        0, 0, 0)
    beta_root = HFDMRG.RHFMovingRoot(copy(root_covariance), 2, 1, 1,
        0, 0, 0)
    identity_link = HFDMRG.RHFStateLink(1, 1,
        Matrix{Float64}(I, 2, 2),
        HFDMRG.StateSelection(0, 2, 0, 0.0, :exact))
    alpha_root_workspace = HFDMRG.RHFRootWorkspace(zeros(4, 4), zeros(4, 4))
    beta_root_workspace = HFDMRG.RHFRootWorkspace(zeros(4, 4), zeros(4, 4))
    HFDMRG._close_root!(alpha_root, identity_link, 2, true,
        alpha_root_workspace; exact=true)
    HFDMRG._open_root!(alpha_root, identity_link, true, alpha_root_workspace)
    HFDMRG._close_root!(beta_root, identity_link, 2, true,
        beta_root_workspace; exact=true)
    HFDMRG._open_root!(beta_root, identity_link, true, beta_root_workspace)
    @test uhf_root_allocations(alpha_root, beta_root, identity_link,
        alpha_root_workspace, beta_root_workspace) == (0, 0, 0, 0)
    @test BLAS.get_num_threads() == threads
end
