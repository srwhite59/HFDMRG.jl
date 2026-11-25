function compact_hopping_fixture(i, j)
    sites = 72
    bands = zeros(j - i + 1, sites)
    bands[end, i] = -0.75
    dense = zeros(sites, sites)
    dense[i, j] = dense[j, i] = -0.75
    operator = UnitCellInteraction(sites, 1, [0.0], zeros(0, 0),
        zeros(18, 1), zeros(1, 18), zeros(171))
    BandedOneBody(bands), operator, dense
end

@testset "compact projected one-body hopping" begin
    for (i, j) in ((18, 55), (18, 19))
        @testset "hop $i to $j" begin
            one_body, operator, dense = compact_hopping_fixture(i, j)
            provenance = UInt(1)
            left_map, right_map = zeros(18, 1), zeros(18, 1)
            left_map[end, 1] = right_map[1, 1] = 1.0
            selection = HFDMRG.StateSelection(0, 1, 17, 0.0, :exact)
            left_link = HFDMRG.RHFStateLink(0, 18, left_map, selection)
            right_link = HFDMRG.RHFStateLink(0, 18, right_map, selection)
            left_rows, right_rows = 1:18, 55:72
            basis = Matrix{Float64}(I, 72, 72)[:, 18:55]
            projected = basis' * dense * basis
            bonding = (basis[i, :] + basis[j, :]) / sqrt(2)
            density = bonding * bonding'
            rank = size(basis, 2)

            empty = HFDMRG._empty_outer_block(operator, provenance)
            seed_left = HFDMRG.seed_cell_block(one_body, operator, 1; provenance)
            seed_right = HFDMRG.seed_cell_block(one_body, operator, 4; provenance)
            builder = HFDMRG.BlockBuildWorkspace(18, 1)
            left = HFDMRG.build_outer_block(empty, seed_left, left_link,
                one_body, operator, builder)
            right = HFDMRG.build_outer_block(seed_right, empty, right_link,
                one_body, operator, builder)
            rhf = HFDMRG.FastRHFPreparedCenter(rank, operator)
            HFDMRG._prepare_lifecycle_center!(rhf, left, right, 2,
                19:36, 37:54, one_body, operator)
            @test rhf.h1[1:rank, 1:rank] == projected
            @test HFDMRG._center_energy_fock!(rhf, density, 1).total ≈
                -1.5 atol=1e-12 rtol=0

            spin_builder = HFDMRG._IndependentSpinEntryBuildWorkspace(
                1, 1, 18, operator)
            spin_left = HFDMRG._build_independent_spin_entry!(spin_builder,
                nothing, left_link, left_link, left_rows, one_body, operator,
                true, provenance)
            spin_right = HFDMRG._build_independent_spin_entry!(spin_builder,
                nothing, right_link, right_link, right_rows, one_body,
                operator, false, provenance)
            uhf = HFDMRG.FastUHFPreparedCenter(rank, rank, operator)
            HFDMRG.prepare_uhf_lifecycle_center!(uhf, spin_left, spin_right,
                2, 19:36, 37:54, one_body, operator)
            @test uhf.alpha.h1[1:rank, 1:rank] == projected
            @test uhf.beta.h1[1:rank, 1:rank] == projected
            @test HFDMRG.uhf_center_energy_fock!(uhf, density, density,
                1, 1).total ≈ -1.5 atol=1e-12 rtol=0
        end
    end
end

@testset "compact hopping energy on the returned determinant" begin
    for (i, j) in ((18, 55), (18, 19)), spin in (:rhf, :uhf),
            starting_orientation in (:physical, :reflected)
        @testset "hop $i to $j, $spin, $starting_orientation" begin
            one_body, operator, dense = compact_hopping_fixture(i, j)
            # The first atom contains both hopping sites. In particular,
            # the UHF Néel initializer occupies the bonding alpha orbital.
            result = solve_hfdmrg(one_body, operator; spin, starting_orientation,
                atom_intervals=[1:55, 56:72], exact=true, state_maxdim=2,
                state_cutoff=0.0, energy_tolerance=1e-11, maximum_half_sweeps=4)
            lifecycle = result.state.lifecycle
            alpha = zeros(operator.sites, 1)
            beta = similar(alpha)
            if spin === :rhf
                workspace = HFDMRG.RHFLifecycleWorkspace(one_body, operator, 2)
                HFDMRG.reconstruct_terminal!(alpha, lifecycle, workspace)
                copyto!(beta, alpha)
            else
                workspace = HFDMRG.UHFLifecycleWorkspace(one_body, operator, 2, 2)
                HFDMRG.reconstruct_uhf_terminal!(alpha, beta, lifecycle, workspace)
            end
            reference = dot(alpha, dense * alpha) + dot(beta, dense * beta)
            # Require at least one bond so an all-zero state energy cannot
            # pass the consistency check. RHF puts both spins in that bond.
            if spin === :rhf
                @test reference ≈ -1.5 atol=1e-12 rtol=0
            else
                @test reference <= -0.75 + 1e-12
            end
            @test result.represented_energy ≈ reference atol=1e-12 rtol=0
        end
    end
end
