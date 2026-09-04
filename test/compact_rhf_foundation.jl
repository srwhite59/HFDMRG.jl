# Packet 10A1 compact common/RHF numerical foundation gates.
function nucleus_fixture(atoms::Int)
    sites = 18atoms
    bands = zeros(2, sites)
    @inbounds for site = 1:sites
        phase = 1 + (site - 1) % 18
        bands[1, site] = 0.001 * (phase - 9.5)^2
        site < sites && (bands[2, site] = site % 18 == 0 ? -0.035 : -0.12)
    end
    one_body = HFDMRG.BandedOneBody(bands)
    decay, amplitude, onsite = 0.24, 0.018, 0.08
    transfer = [exp(-18decay)]
    source = reshape([amplitude * exp(decay * phase) for phase = 1:18], 18, 1)
    sink = reshape([exp(-decay * (18 + phase)) for phase = 1:18], 1, 18)
    triangle = [first == second ? onsite :
        amplitude * exp(-decay * (second - first))
        for second = 1:18 for first = 1:second]
    operator = HFDMRG.UnitCellInteraction(sites, 1, transfer, zeros(0, 0),
        source, sink, triangle)
    one_body, operator
end

function dense_h1(one_body)
    sites = HFDMRG.h1_sites(one_body)
    [HFDMRG.h1_entry(one_body, i, j) for i = 1:sites, j = 1:sites]
end

function dense_interaction(operator)
    rank = HFDMRG._maximum_channel_rank(operator)
    first = zeros(rank); second = zeros(rank)
    [HFDMRG._interaction_entry(operator, i, j, first, second)
        for i = 1:operator.sites, j = 1:operator.sites]
end

function deterministic_isometry(rows, columns, seed)
    matrix = [sin(0.17 * (row + seed) * column) +
        cos(0.11 * row * (column + seed))
        for row = 1:rows, column = 1:columns]
    Matrix(qr(matrix).Q[:, 1:columns])
end

function grow_cells(one_body, operator, cells, target_rank, provenance)
    block = HFDMRG.seed_cell_block(one_body, operator, first(cells);
        provenance=provenance)
    basis = Matrix{Float64}(I, 18, 18)
    workspace = HFDMRG.BlockBuildWorkspace(max(36, target_rank + 18),
        HFDMRG._maximum_channel_rank(operator))
    for cell in Iterators.drop(cells, 1)
        physical = HFDMRG.seed_cell_block(one_body, operator, cell;
            provenance=provenance)
        rows = block.rank + physical.rank
        rank = min(target_rank, rows)
        map = deterministic_isometry(rows, rank, cell)
        link = HFDMRG.RHFStateLink(block.rank, physical.rank, map,
            HFDMRG.StateSelection(0, rank, rows-rank,
                rows == rank ? 0.0 : 1.0e-12, :maximum_rank))
        block = HFDMRG.build_outer_block(block, physical, link, one_body,
            operator, workspace)
        old_sites = size(basis, 1)
        raw = zeros(old_sites + 18, rows)
        raw[1:old_sites, 1:size(basis, 2)] .= basis
        @inbounds for index = 1:18
            raw[old_sites + index, size(basis, 2) + index] = 1.0
        end
        basis = raw * map
    end
    block, basis
end

function center_oracle(one_body, operator, blocks, bases, occupied)
    sites = operator.sites
    rank = sum(block.rank for block in blocks)
    basis = zeros(sites, rank)
    offset = 0
    for (block, local_basis) in zip(blocks, bases)
        rows = (first(block.cells)-1)*18+1:last(block.cells)*18
        basis[rows, offset+1:offset+block.rank] .= local_basis
        offset += block.rank
    end
    orbitals = deterministic_isometry(rank, occupied, 41)
    density = orbitals * transpose(orbitals)
    physical_density = basis * density * transpose(basis)
    h1 = dense_h1(one_body)
    interaction = dense_interaction(operator)
    physical_fock = h1 + Diagonal(2interaction * diag(physical_density)) -
        interaction .* physical_density
    fock = transpose(basis) * physical_fock * basis
    one = 2dot(h1, physical_density)
    hartree = 2dot(diag(physical_density), interaction * diag(physical_density))
    exchange = -dot(physical_density, interaction .* physical_density)
    density, (total=one+hartree+exchange, one_body=one,
        hartree=hartree, exchange=exchange), fock, basis
end

function warmed_nonfactor_allocation(prepared, left, right, cell,
        one_body, operator, density, occupied, output, input)
    HFDMRG._prepare_lifecycle_center!(prepared, left, right, cell,
        one_body, operator)
    HFDMRG._center_energy_fock!(prepared, density, occupied)
    HFDMRG.apply_center_fock!(output, prepared, input)
    @allocated begin
        HFDMRG._prepare_lifecycle_center!(prepared, left, right, cell,
            one_body, operator)
        HFDMRG._center_energy_fock!(prepared, density, occupied)
        HFDMRG.apply_center_fock!(output, prepared, input)
        nothing
    end
end

@testset "compact RHF block-state nucleus" begin
    for atoms in (10, 20)
        one_body, operator = nucleus_fixture(atoms)
        provenance = UInt(atoms)
        center_cell = atoms ÷ 2
        left, left_basis = grow_cells(one_body, operator, 1:center_cell-1,
            3, provenance)
        first_center = HFDMRG.seed_cell_block(one_body, operator, center_cell;
            provenance)
        second_center = HFDMRG.seed_cell_block(one_body, operator,
            center_cell + 1; provenance)
        right, right_basis = grow_cells(one_body, operator,
            center_cell+2:atoms, 4, provenance)
        physical_basis = Matrix{Float64}(I, 18, 18)
        rank = left.rank + 36 + right.rank
        occupied = min(atoms ÷ 2, rank)
        density, oracle_energy, oracle_fock, _ = center_oracle(one_body,
            operator, (left, first_center, second_center, right),
            (left_basis, physical_basis, physical_basis, right_basis), occupied)
        prepared = HFDMRG.FastRHFPreparedCenter(rank, operator;
            maximum_rhs=3, maximum_side_rank=max(left.rank, right.rank))
        HFDMRG._prepare_lifecycle_center!(prepared, left, right, center_cell,
            one_body, operator)
        energy = HFDMRG._center_energy_fock!(prepared, density, occupied)
        @test abs(energy.total - oracle_energy.total) <= 2e-11
        @test abs(energy.one_body - oracle_energy.one_body) <= 2e-12
        @test abs(energy.hartree - oracle_energy.hartree) <= 2e-11
        @test abs(energy.exchange - oracle_energy.exchange) <= 2e-11
        @test maximum(abs, prepared.fock[1:rank, 1:rank] - oracle_fock) <= 2e-11
        panel = reshape(collect(range(-0.3, 0.9; length=3rank)), rank, 3)
        panel_output = similar(panel)
        HFDMRG.apply_center_fock!(panel_output, prepared, panel)
        @test panel_output ≈ oracle_fock * panel atol=2e-11 rtol=2e-11
        @test warmed_nonfactor_allocation(prepared, left, right, center_cell,
            one_body, operator, density, occupied, panel_output, panel) == 0
        larger = HFDMRG.FastRHFPreparedCenter(rank + 3, operator;
            maximum_rhs=4, maximum_side_rank=max(left.rank, right.rank))
        inactive_before = reinterpret(UInt64, copy(larger.h1[rank+1:end, :]))
        HFDMRG._prepare_lifecycle_center!(larger, left, right, center_cell,
            one_body, operator)
        @test reinterpret(UInt64, copy(larger.h1[rank+1:end, :])) ==
            inactive_before
        @test HFDMRG.reflect_block(HFDMRG.reflect_block(left, atoms), atoms).cells ==
            left.cells
        if atoms == 10
            first_cell = HFDMRG.seed_cell_block(one_body, operator, 1)
            second_cell = HFDMRG.seed_cell_block(one_body, operator, 2)
            map = deterministic_isometry(36, 3, 7)
            link = HFDMRG.RHFStateLink(18, 18, map,
                HFDMRG.StateSelection(0, 3, 33, 1e-12, :maximum_rank))
            builder = HFDMRG.BlockBuildWorkspace(36, 1)
            HFDMRG.build_outer_block(first_cell, second_cell, link,
                one_body, operator, builder)
            publication = @allocated published = HFDMRG.build_outer_block(
                first_cell, second_cell, link, one_body, operator, builder)
            durable_arrays = sizeof(published.h1) + sizeof(published.pair_field) +
                sizeof(published.left_feature) + sizeof(published.right_feature) +
                sizeof(published.left_edge) + sizeof(published.right_edge)
            @test durable_arrays <= publication <= Base.summarysize(published)
            @test @allocated(HFDMRG._fill_merged!(builder, first_cell,
                second_cell, one_body, operator)) == 0
        end
    end
end

@testset "terminal, exact, reflection, and ownership" begin
    one_body, operator = nucleus_fixture(2)
    first = HFDMRG.seed_cell_block(one_body, operator, 1; provenance=UInt(2))
    second = HFDMRG.seed_cell_block(one_body, operator, 2; provenance=UInt(2))
    density, oracle_energy, oracle_fock, _ = center_oracle(one_body, operator,
        (first, second), (Matrix{Float64}(I, 18, 18),
            Matrix{Float64}(I, 18, 18)), 1)
    empty = HFDMRG._empty_outer_block(operator, UInt(2))
    prepared = HFDMRG.FastRHFPreparedCenter(36, operator; maximum_rhs=2)
    HFDMRG._prepare_lifecycle_center!(prepared, empty, empty, 1,
        one_body, operator)
    energy = HFDMRG._center_energy_fock!(prepared, density, 1)
    @test energy.total ≈ oracle_energy.total atol=2e-11
    @test prepared.fock[1:36, 1:36] ≈ oracle_fock atol=2e-11
    @test reverse(reverse(first)).h1 == first.h1
    @test HFDMRG.reflect_block(HFDMRG.reflect_block(first, 2), 2).cells ==
        first.cells
    durable = sum(sizeof, (first.h1, first.pair_field, first.left_feature,
        first.right_feature, first.left_core, first.right_core,
        first.left_edge, first.right_edge, first.link.close_map,
        first.link.completed_map))
    @test HFDMRG.block_storage_bytes(first) == durable
    @test !any(size(array, 1) == operator.sites && size(array, 2) > 18
        for array in (first.h1, first.pair_field, first.left_feature,
            first.right_feature))
end

@testset "validation, zero rank, aliases, and BLAS" begin
    threads = BLAS.get_num_threads()
    selection = HFDMRG.StateSelection(0, 0, 0, 0.0, :empty)
    link = HFDMRG.RHFStateLink(0, 0, zeros(0, 0), selection)
    @test link.active_rank == 0
    @test_throws ArgumentError HFDMRG.StateSelection(-1, 0, 0, 0.0, :bad)
    @test_throws ArgumentError HFDMRG.RHFRoot([0.8 0.0; 0.0 0.0]; occupied=1)
    one_body, operator = nucleus_fixture(2)
    block = HFDMRG.seed_cell_block(one_body, operator, 1)
    empty = HFDMRG._empty_outer_block(operator, UInt(0))
    prepared = HFDMRG.FastRHFPreparedCenter(54, operator;
        maximum_side_rank=18)
    @test_throws ArgumentError HFDMRG._prepare_lifecycle_center!(prepared,
        block, empty, 1, one_body, operator)
    snapshot = copy(prepared.h1)
    @test_throws ArgumentError HFDMRG._prepare_lifecycle_center!(prepared,
        block, empty, 1, one_body, operator)
    @test prepared.h1 == snapshot
    other = HFDMRG.seed_cell_block(one_body, operator, 2; provenance=UInt(1))
    @test_throws ArgumentError HFDMRG._prepare_lifecycle_center!(prepared,
        empty, other, 1, one_body, operator)
    @test prepared.h1 == snapshot
    vector = ones(36)
    prepared.active_rank = 36
    fill!(@view(prepared.fock[1:36, 1:36]), 0.0)
    @test_throws ArgumentError HFDMRG.apply_center_fock!(vector, prepared,
        vector)
    @test BLAS.get_num_threads() == threads

    residual_operator = HFDMRG.UnitCellInteraction(operator.sites, 1,
        Float64[], reshape([operator.tail_transfer[1]], 1, 1),
        operator.phase_source, operator.phase_sink, operator.local_triangle)
    residual_prepared = HFDMRG.FastRHFPreparedCenter(36, residual_operator)
    residual_empty = HFDMRG._empty_outer_block(residual_operator, UInt(0))
    HFDMRG._prepare_lifecycle_center!(residual_prepared, residual_empty,
        residual_empty, 1, one_body, residual_operator)
    @test @allocated(HFDMRG._prepare_lifecycle_center!(residual_prepared,
        residual_empty, residual_empty, 1, one_body, residual_operator)) == 0
end
