# Packet 10A7C: one-cell compact-RHF residual-channel propagation.

function residual_fragment_fixture(groups::Int)
    groups >= 4 && iseven(groups) || throw(ArgumentError(
        "fixture requires an even group count of at least four"))
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
    intervals = HFDMRG._traversal_intervals(operator)
    atoms = groups == 4 ? [1:21, 22:39, 40:48, 49:57] : copy(intervals)
    one_body, operator, intervals, atoms
end

function residual_center_oracle(one_body, operator, left_basis, occupied)
    sites = operator.sites
    rank = size(left_basis, 2) + sites - size(left_basis, 1)
    basis = zeros(sites, rank)
    basis[1:size(left_basis, 1), 1:size(left_basis, 2)] .= left_basis
    offset = size(left_basis, 2)
    @inbounds for row = size(left_basis, 1)+1:sites
        basis[row, offset + row - size(left_basis, 1)] = 1.0
    end
    orbitals = deterministic_isometry(rank, occupied, 73)
    density = orbitals * transpose(orbitals)
    physical_density = basis * density * transpose(basis)
    h1 = dense_h1(one_body)
    interaction = dense_interaction(operator)
    physical_fock = h1 + Diagonal(2interaction * diag(physical_density)) -
        interaction .* physical_density
    energy = 2dot(h1, physical_density) +
        2dot(diag(physical_density), interaction * diag(physical_density)) -
        dot(physical_density, interaction .* physical_density)
    density, energy, transpose(basis) * physical_fock * basis
end

@testset "compact RHF one-face residual transfer parity" begin
    one_body, operator, intervals, atoms = residual_fragment_fixture(4)
    provenance = UInt(0x10a7c)
    first = HFDMRG.seed_interval_block(one_body, operator, 1, intervals[1];
        provenance)
    second = HFDMRG.seed_interval_block(one_body, operator, 2, intervals[2];
        provenance)
    map = deterministic_isometry(21, 3, 29)
    link = HFDMRG.RHFStateLink(3, 18, map,
        HFDMRG.StateSelection(0, 3, 18, 0.0, :maximum_rank))
    full_workspace = HFDMRG.BlockBuildWorkspace(21, 3)
    side_workspace = HFDMRG.BlockBuildWorkspace(21, 3)
    full = HFDMRG.build_outer_block(first, second, link, one_body, operator,
        full_workspace)
    side = HFDMRG.build_outer_block(first, second, link, one_body, operator,
        side_workspace; required_face=:right)
    @test side.h1 ≈ full.h1 atol=2e-13 rtol=2e-13
    @test side.pair_field ≈ full.pair_field atol=2e-13 rtol=2e-13
    @test side.right_feature ≈ full.right_feature atol=2e-13 rtol=2e-13
    @test side.right_core ≈ full.right_core atol=2e-13 rtol=2e-13
    @test iszero(side.left_feature)
    @test iszero(side.left_core)
    @test side_workspace.channel_transfer_calls == 1
    @test side_workspace.residual_transfer_steps == 1
    @test side_workspace.maximum_transfer_power == 1
    @test full_workspace.channel_transfer_calls == 2
    @test full_workspace.residual_transfer_steps == 2

    empty = HFDMRG._empty_outer_block(operator, provenance)
    rank = side.rank + length(intervals[3]) + length(intervals[4])
    density, oracle_energy, oracle_fock = residual_center_oracle(one_body,
        operator, map, 2)
    prepared_full = HFDMRG.FastRHFPreparedCenter(rank, operator;
        maximum_rhs=2, maximum_side_rank=3)
    prepared_side = HFDMRG.FastRHFPreparedCenter(rank, operator;
        maximum_rhs=2, maximum_side_rank=3)
    for (prepared, block) in ((prepared_full, full), (prepared_side, side))
        HFDMRG._prepare_lifecycle_center!(prepared, block, empty, 3,
            intervals[3], intervals[4], one_body, operator)
    end
    energy_full = HFDMRG._center_energy_fock!(prepared_full, density, 2)
    energy_side = HFDMRG._center_energy_fock!(prepared_side, density, 2)
    @test energy_side.total ≈ energy_full.total atol=3e-12
    @test energy_side.total ≈ oracle_energy atol=3e-11
    @test prepared_side.fock[1:rank, 1:rank] ≈
        prepared_full.fock[1:rank, 1:rank] atol=3e-12 rtol=3e-12
    @test prepared_side.fock[1:rank, 1:rank] ≈ oracle_fock atol=3e-11 rtol=3e-11
    panel = reshape(collect(range(-0.2, 0.4; length=2rank)), rank, 2)
    output = similar(panel)
    HFDMRG.apply_center_fock!(output, prepared_side, panel)
    @test output ≈ oracle_fock * panel atol=3e-11 rtol=3e-11

    reflected = HFDMRG.reflect_block(side, length(intervals), operator.sites)
    roundtrip = HFDMRG.reflect_block(reflected, length(intervals),
        operator.sites)
    @test roundtrip.right_feature == side.right_feature
    @test roundtrip.left_feature == side.left_feature
    @test roundtrip.right_core == side.right_core
    @test roundtrip.left_core == side.left_core

    calls = side_workspace.channel_transfer_calls
    steps = side_workspace.residual_transfer_steps
    scratch = copy(side_workspace.h1)
    @test_throws ArgumentError HFDMRG._fill_merged!(side_workspace, first,
        second, one_body, operator, :outward)
    @test side_workspace.channel_transfer_calls == calls
    @test side_workspace.residual_transfer_steps == steps
    @test side_workspace.h1 == scratch
    too_small = HFDMRG.BlockBuildWorkspace(20, 3)
    @test_throws DimensionMismatch HFDMRG._fill_merged!(too_small, first,
        second, one_body, operator, :right)
    @test too_small.channel_transfer_calls == 0
end

@testset "compact RHF residual lifecycle and producer observation" begin
    one_body, operator, intervals, atoms = residual_fragment_fixture(4)
    coefficients = reshape(HFDMRG._dimer_orbital(one_body, 1:39), :, 1)
    second = reshape(HFDMRG._dimer_orbital(one_body, 40:57), :, 1)
    coefficients = [coefficients zeros(39); zeros(18) second]
    projector, oracle_energy, _ = dense_rhf_state(one_body, operator,
        coefficients)
    for forward in (true, false)
        control = HFDMRG.RHFStateControl(2, 1.0e-12; exact=true)
        workspace = HFDMRG.RHFLifecycleWorkspace(one_body, operator, 2;
            maximum_rhs=2)
        lifecycle = HFDMRG.initialize_rhf_dimer_lifecycle(one_body,
            operator, control, workspace; forward, atom_intervals=atoms)
        initialization_calls = workspace.builder.channel_transfer_calls
        initialization_steps = workspace.builder.residual_transfer_steps
        @test initialization_calls == length(intervals) - 2
        @test initialization_steps == initialization_calls
        @test workspace.builder.maximum_transfer_power == 1
        evidence = NamedTuple[]
        calls_before = workspace.builder.channel_transfer_calls
        steps_before = workspace.builder.residual_transfer_steps
        for _ = 1:4
            append!(evidence, HFDMRG.run_fixed_half_sweep!(lifecycle,
                one_body, operator, control, workspace))
        end
        builds = length(evidence)
        @test workspace.builder.channel_transfer_calls - calls_before == builds
        @test workspace.builder.residual_transfer_steps - steps_before == builds
        @test workspace.builder.maximum_transfer_power == 1
        @test maximum(abs(item.energy-oracle_energy) for item in evidence) <=
            8e-11
        reconstructed = zeros(operator.sites, 2)
        HFDMRG.reconstruct_terminal!(reconstructed, lifecycle, workspace)
        @test norm(reconstructed * transpose(reconstructed) - projector, Inf) <=
            5e-11
        result = HFDMRG.HFDMRGResult(true, :fixed_state, oracle_energy,
            lifecycle.half_sweeps, 0, 0, 0,
            maximum(block.rank for block in lifecycle.collection), 0, 0.0,
            Float64[], 0.0, 0.0, 0.0, HFDMRG.RHFCompactState(lifecycle))
        producer = HFDMRG.ProducerHamiltonian(dense_h1(one_body),
            dense_interaction(operator); h1_bandwidth=3)
        root_before = copy(lifecycle.root.covariance)
        observed = HFDMRG.producer_energy(result, producer)
        @test abs(observed.total-oracle_energy) <= observed.error_envelope
        @test isequal(root_before, lifecycle.root.covariance)
    end
end
