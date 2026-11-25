@testset "Historical sweep energy and energy zero" begin
    N = 12
    H = diagm(0 => range(-0.2, 0.3; length=N),
        1 => fill(-0.12, N-1), -1 => fill(-0.12, N-1))
    Hdn = H + Diagonal(range(-0.05, 0.05; length=N))
    V = [5 / (1 + abs(i-j)) for i = 1:N, j = 1:N]
    rng = MersenneTwister(7)
    up = Matrix(qr(randn(rng, N, 2)).Q)
    dn = Matrix(qr(randn(rng, N, 1)).Q)
    layout = HFDMRG.SliceLayout(fill(2, N ÷ 2))
    # The sliced tensor encodes the same density-density Hamiltonian.
    V6 = zeros(2, 2, 2, 2, N ÷ 2, N ÷ 2)
    for t = 1:N÷2, s = 1:N÷2, b = 1:2, a = 1:2
        V6[a, a, b, b, s, t] = V[2s-2+a, 2t-2+b]
    end
    sliced = HFDMRG.SlicedBasisBackend(layout, V6)
    function physical_energy(Hu, Hd, a, b)
        Pa, Pb = a * a', b * b'
        q = diag(Pa + Pb)
        dot(Hu, Pa) + dot(Hd, Pb) + 0.5dot(q, V * q) -
            0.5sum(V .* (Pa.^2 + Pb.^2))
    end
    for (kind, interaction) in ((:density, V), (:sliced, sliced)),
            spin in (:rhf, :common_uhf, :split_uhf)
        @testset "$kind $spin" begin
            partition = kind === :sliced ? layout : nothing
            function run(shift=0.0; kwargs...)
                shifted_up, shifted_dn = H + shift*I, Hdn + shift*I
                kw = (; blocksize=2, block_partition=partition, kwargs...)
                spin === :rhf &&
                    return solve_hfdmrg(shifted_up, interaction, up; kw...)
                spin === :common_uhf &&
                    return solve_hfdmrg(shifted_up, interaction, up, dn; kw...)
                solve_hfdmrg(shifted_up, shifted_dn, interaction, up, dn; kw...)
            end
            Hd = spin === :split_uhf ? Hdn : H
            @test_throws ArgumentError run(; maxiter=0)
            @test_throws ArgumentError run(; maxiter=-1)

            diagnostics = HFDMRG._RunDiagnostics()
            budget = run(; maxiter=1, cutoff=0.0, scf_cutoff=0.0,
                _diagnostics=diagnostics)
            @test budget[3] ≈ physical_energy(H, Hd, budget[1], budget[2]) atol=2e-12
            if kind === :density && spin === :rhf
                @test abs(last(diagnostics.windows).local_energies[end] -
                    budget[3]) > 1e-3
            end

            seen = HFDMRG.SweepInfo[]
            observer = info -> begin
                push!(seen, info)
                @test info.energy ≈ physical_energy(H, Hd, info.psiup,
                    info.psidn) atol=2e-12
                info.sweep == 2
            end
            stopped = run(; maxiter=4, cutoff=0.0, scf_cutoff=0.0, observer)
            @test length(seen) == 2
            @test stopped == (seen[end].psiup, seen[end].psidn, seen[end].energy)

            # Use one local update here to isolate the energy-zero test from
            # roundoff-scale rise decisions in a multi-update local SCF.
            first = run(; maxiter=1, cutoff=0.0, scf_cutoff=Inf)
            second = run(; maxiter=2, cutoff=0.0, scf_cutoff=Inf)
            electrons = size(first[1], 2) + size(first[2], 2)
            shift = (10000.0 - first[3]) / electrons
            shifted_sweeps = HFDMRG.SweepInfo[]
            shifted_observer = info -> (push!(shifted_sweeps, info); false)
            shifted = run(shift; maxiter=2, cutoff=1e-8, scf_cutoff=Inf,
                observer=shifted_observer)
            @test length(shifted_sweeps) == 2
            @test !shifted_sweeps[1].converged
            @test shifted_sweeps[1].energy ≈ 10000.0 atol=1e-8 rtol=0
            @test shifted[3] - electrons*shift ≈ second[3] atol=2e-8 rtol=0
            @test shifted[1] * shifted[1]' ≈ second[1] * second[1]' atol=2e-8
            @test shifted[2] * shifted[2]' ≈ second[2] * second[2]' atol=2e-8

            if kind === :density && spin === :rhf
                large = HFDMRG._RunDiagnostics()
                large_shift = 3.0e9
                run(large_shift; maxiter=1, scf_cutoff=Inf, _diagnostics=large)
                @test all(w -> only(w.local_energies) > 1e10 &&
                    isempty(w.rise_updates) && w.damping_after == 1,
                    large.windows)
            end
        end
    end
end
