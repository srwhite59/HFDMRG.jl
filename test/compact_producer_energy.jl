using Serialization

function _fixed_result(lifecycle::HFDMRG.RHFFixedLifecycle)
    HFDMRG.HFDMRGResult(true, :fixed_state, 0.0, lifecycle.half_sweeps,
        0, 0, 0, maximum(block.rank for block in lifecycle.collection),
        0, 0.0, Float64[], 0.0, 0.0, 0.0,
        HFDMRG.RHFCompactState(lifecycle))
end

function _fixed_result(lifecycle::HFDMRG.UHFFixedLifecycle)
    HFDMRG.UHFDMRGResult(true, :fixed_state, 0.0, lifecycle.half_sweeps,
        0, 0, 0, 0, 0, 0, 0, 0.0, 0.0, Float64[], 0.0, 0.0, 0.0,
        0.0, 0.0, HFDMRG.UHFCompactState(lifecycle))
end

function _dense_producer(coefficients_alpha, coefficients_beta, h1, v,
        constant)
    Pa = coefficients_alpha * transpose(coefficients_alpha)
    Pb = coefficients_beta === nothing ? Pa :
        coefficients_beta * transpose(coefficients_beta)
    da, db = diag(Pa), diag(Pb)
    ha, hb = dot(h1, Pa), dot(h1, Pb)
    hartree = 0.5dot(da+db, v*(da+db))
    ka = -0.5dot(Pa, v .* Pa)
    kb = -0.5dot(Pb, v .* Pb)
    (one_body_alpha=ha, one_body_beta=hb, hartree=hartree,
        exchange_alpha=ka, exchange_beta=kb,
        total=constant+ha+hb+hartree+ka+kb)
end

function _lifecycle_signature(lifecycle)
    io = IOBuffer()
    Serialization.serialize(io, lifecycle)
    take!(io)
end

function _producer_fixture(one_body, operator; constant=0.37)
    h1 = dense_h1(one_body); interaction = dense_interaction(operator)
    HFDMRG.ProducerHamiltonian(h1, interaction; constant,
        h1_bandwidth=1), h1, interaction
end

function _producer_uhf_solve(atoms)
    one_body, operator = nucleus_fixture(atoms)
    result = HFDMRG.solve_hfdmrg(one_body, operator; spin=:uhf,
        initialization=:atomic_neel, state_maxdim=min(4, atoms÷2),
        state_cutoff=1e-12, energy_tolerance=1e-7,
        maximum_half_sweeps=60, starting_orientation=:physical,
        exact=false)
    result, one_body, operator
end

function _producer_fragment_fixture(remainder)
    sites = 2HFDMRG.UNIT_CELL_WIDTH + remainder
    bands = zeros(2, sites)
    @inbounds for site = 1:sites
        phase = 1 + (site-1) % HFDMRG.UNIT_CELL_WIDTH
        bands[1, site] = 0.001 * (phase-9.5)^2
        site < sites && (bands[2, site] = site % 18 == 0 ? -0.035 : -0.12)
    end
    one_body = HFDMRG.BandedOneBody(bands)
    decay, amplitude, onsite = 0.24, 0.018, 0.08
    transfer = [exp(-18decay)]
    source = reshape([amplitude*exp(decay*phase) for phase = 1:18], 18, 1)
    sink = reshape([exp(-decay*(18+phase)) for phase = 1:18], 1, 18)
    triangle = [first == second ? onsite :
        amplitude*exp(-decay*(second-first))
        for second = 1:18 for first = 1:second]
    operator = HFDMRG.UnitCellInteraction(sites, 1, transfer, zeros(0, 0),
        source, sink, triangle)
    one_body, operator, [1:18, 19:sites]
end

function _fixed_producer_case(atoms, spin, forward, exact)
    one_body, operator = nucleus_fixture(atoms)
    maximum = exact ? (atoms+1)÷2 : min(4, (atoms+1)÷2)
    state = HFDMRG.RHFStateControl(maximum, 1e-12; exact)
    if spin === :rhf
        work = HFDMRG.RHFLifecycleWorkspace(one_body, operator, maximum)
        lifecycle = HFDMRG.initialize_rhf_dimer_lifecycle(one_body,
            operator, state, work; forward)
        for _ = 1:4
            HFDMRG.run_fixed_half_sweep!(lifecycle, one_body, operator,
                state, work)
        end
        result = _fixed_result(lifecycle)
        alpha = zeros(operator.sites, lifecycle.root.total_occupied)
        HFDMRG.reconstruct_terminal!(alpha, lifecycle, work)
        beta = nothing
    else
        control = HFDMRG.UHFStateControl(state)
        work = HFDMRG.UHFLifecycleWorkspace(one_body, operator, maximum,
            maximum)
        lifecycle = HFDMRG.initialize_uhf_neel_lifecycle(one_body, operator,
            control, work; forward)
        for _ = 1:4
            HFDMRG.run_uhf_fixed_half_sweep!(lifecycle, one_body, operator,
                control, work)
        end
        result = _fixed_result(lifecycle)
        alpha = zeros(operator.sites, lifecycle.root.alpha.total_occupied)
        beta = zeros(operator.sites, lifecycle.root.beta.total_occupied)
        HFDMRG.reconstruct_uhf_terminal!(alpha, beta, lifecycle, work)
    end
    result, one_body, operator, alpha, beta
end

@testset "block-streamed producer energy basic parity" begin
    for spin in (:rhf, :uhf), forward in (true, false)
        one_body, operator = nucleus_fixture(2)
        producer, h1, interaction = _producer_fixture(one_body, operator)
        state = HFDMRG.RHFStateControl(1, 1e-12; exact=true)
        if spin === :rhf
            work = HFDMRG.RHFLifecycleWorkspace(one_body, operator, 1)
            lifecycle = HFDMRG.initialize_rhf_dimer_lifecycle(one_body,
                operator, state, work; forward)
            for _ = 1:4
                HFDMRG.run_fixed_half_sweep!(lifecycle, one_body, operator,
                    state, work)
            end
            result = _fixed_result(lifecycle)
            coefficients = zeros(operator.sites, 1)
            HFDMRG.reconstruct_terminal!(coefficients, lifecycle, work)
            dense = _dense_producer(coefficients, nothing, h1, interaction,
                producer.constant)
        else
            control = HFDMRG.UHFStateControl(state)
            work = HFDMRG.UHFLifecycleWorkspace(one_body, operator, 1, 1)
            lifecycle = HFDMRG.initialize_uhf_neel_lifecycle(one_body,
                operator, control, work; forward)
            for _ = 1:4
                HFDMRG.run_uhf_fixed_half_sweep!(lifecycle, one_body,
                    operator, control, work)
            end
            result = _fixed_result(lifecycle)
            alpha = zeros(operator.sites, 1); beta = similar(alpha)
            HFDMRG.reconstruct_uhf_terminal!(alpha, beta, lifecycle, work)
            dense = _dense_producer(alpha, beta, h1, interaction,
                producer.constant)
        end
        workspace = HFDMRG.ProducerEnergyWorkspace(result, producer)
        before = _lifecycle_signature(lifecycle)
        observed = HFDMRG.producer_energy(result, producer, workspace)
        repeated = HFDMRG.producer_energy(result, producer, workspace)
        @test observed.total ≈ dense.total atol=observed.error_envelope
        @test observed.one_body_alpha ≈ dense.one_body_alpha atol=2e-11
        @test observed.one_body_beta ≈ dense.one_body_beta atol=2e-11
        @test observed.hartree ≈ dense.hartree atol=2e-11
        @test observed.exchange_alpha ≈ dense.exchange_alpha atol=2e-11
        @test observed.exchange_beta ≈ dense.exchange_beta atol=2e-11
        @test isequal(observed.total, repeated.total)
        @test isequal(before, _lifecycle_signature(lifecycle))
    end
end

@testset "H10/H20 exact and finite producer decomposition" begin
    for (atoms, exact, forward) in ((10, false, true),
            (10, true, false), (20, false, false), (20, true, true)),
            spin in (:rhf, :uhf)
        result, one_body, operator, alpha, beta = _fixed_producer_case(
            atoms, spin, forward, exact)
        producer, h1, interaction = _producer_fixture(one_body, operator;
            constant=-0.19)
        observed = HFDMRG.producer_energy(result, producer)
        dense = _dense_producer(alpha, beta, h1, interaction,
            producer.constant)
        @test abs(observed.total-dense.total) <= observed.error_envelope
        @test observed.one_body_alpha ≈ dense.one_body_alpha atol=3e-10
        @test observed.one_body_beta ≈ dense.one_body_beta atol=3e-10
        @test observed.hartree ≈ dense.hartree atol=3e-10
        @test observed.exchange_alpha ≈ dense.exchange_alpha atol=3e-10
        @test observed.exchange_beta ≈ dense.exchange_beta atol=3e-10
    end
end

@testset "nonorthogonal UHF, spin swap, and state observation" begin
    result, one_body, operator = _producer_uhf_solve(10)
    producer, h1, interaction = _producer_fixture(one_body, operator;
        constant=0.11)
    lifecycle = result.state.lifecycle
    work = HFDMRG.UHFLifecycleWorkspace(one_body, operator, 4, 4)
    alpha = zeros(operator.sites, lifecycle.root.alpha.total_occupied)
    beta = zeros(operator.sites, lifecycle.root.beta.total_occupied)
    HFDMRG.reconstruct_uhf_terminal!(alpha, beta, lifecycle, work)
    @test norm(transpose(alpha)*beta) > 1e-5
    dense = _dense_producer(alpha, beta, h1, interaction,
        producer.constant)
    workspace = HFDMRG.ProducerEnergyWorkspace(result, producer)
    signature = _lifecycle_signature(lifecycle)
    observed = HFDMRG.producer_energy(result, producer, workspace)
    allocation = @allocated HFDMRG.producer_energy(result, producer,
        workspace)
    @test abs(observed.total-dense.total) <= observed.error_envelope
    @test _lifecycle_signature(lifecycle) == signature
    @test allocation <= 4096
    threads = BLAS.get_num_threads()
    HFDMRG.producer_energy(result, producer, workspace)
    @test BLAS.get_num_threads() == threads

    HFDMRG._spin_swap_uhf_lifecycle!(lifecycle, operator)
    swapped = HFDMRG.producer_energy(result, producer)
    @test swapped.total == observed.total
    @test swapped.one_body_alpha == observed.one_body_beta
    @test swapped.exchange_alpha == observed.exchange_beta
end

@testset "terminal fragments, zero ranks, and hostile workspace" begin
    for spin in (:rhf, :uhf), forward in (true, false)
        one_body, operator, atoms = _producer_fragment_fixture(3)
        state = HFDMRG.RHFStateControl(1, 1e-12; exact=true)
        if spin === :rhf
            work = HFDMRG.RHFLifecycleWorkspace(one_body, operator, 1)
            lifecycle = HFDMRG.initialize_rhf_dimer_lifecycle(one_body,
                operator, state, work; forward, atom_intervals=atoms)
            result = _fixed_result(lifecycle)
            alpha = zeros(operator.sites, 1)
            HFDMRG.reconstruct_terminal!(alpha, lifecycle, work)
            beta = nothing
        else
            control = HFDMRG.UHFStateControl(state)
            work = HFDMRG.UHFLifecycleWorkspace(one_body, operator, 1, 1)
            lifecycle = HFDMRG.initialize_uhf_neel_lifecycle(one_body,
                operator, control, work; forward, atom_intervals=atoms)
            result = _fixed_result(lifecycle)
            alpha = zeros(operator.sites, 1); beta = similar(alpha)
            HFDMRG.reconstruct_uhf_terminal!(alpha, beta, lifecycle, work)
        end
        producer, h1, interaction = _producer_fixture(one_body, operator)
        observed = HFDMRG.producer_energy(result, producer)
        dense = _dense_producer(alpha, beta, h1, interaction,
            producer.constant)
        @test abs(observed.total-dense.total) <= observed.error_envelope
        @test length(lifecycle.intervals[end]) == 3
        @test any(block -> spin === :rhf ? block.rank == 0 :
            min(block.alpha.rank, block.beta.rank) == 0,
            lifecycle.collection)
    end

    result, one_body, operator, _, _ = _fixed_producer_case(2, :rhf,
        true, true)
    producer, _, interaction = _producer_fixture(one_body, operator)
    workspace = HFDMRG.ProducerEnergyWorkspace(result, producer)
    signature = _lifecycle_signature(result.state.lifecycle)
    resize!(workspace.density_alpha, producer.sites-1)
    @test_throws DimensionMismatch HFDMRG.producer_energy(result, producer,
        workspace)
    @test _lifecycle_signature(result.state.lifecycle) == signature
    coordinate_hostile = HFDMRG.ProducerEnergyWorkspace(result, producer)
    coordinate_hostile.covariance = zeros(0, 0)
    @test_throws DimensionMismatch HFDMRG.producer_energy(result, producer,
        coordinate_hostile)
    @test _lifecycle_signature(result.state.lifecycle) == signature
    stale = HFDMRG.ProducerEnergyWorkspace(result, producer)
    stale.replacements += 1
    @test_throws ArgumentError HFDMRG.producer_energy(result, producer, stale)
    bad_interaction = copy(interaction); bad_interaction[1, 1] = NaN
    bad = HFDMRG.ProducerHamiltonian(producer.one_body, bad_interaction;
        constant=producer.constant, h1_bandwidth=producer.h1_bandwidth)
    @test_throws ArgumentError HFDMRG.producer_energy(result, bad)
    @test_throws DimensionMismatch HFDMRG.ProducerHamiltonian(zeros(2, 3),
        zeros(2, 3); h1_bandwidth=1)
end
