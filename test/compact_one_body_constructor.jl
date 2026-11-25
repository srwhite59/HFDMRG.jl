@testset "BandedOneBody validates and owns matrix and view inputs" begin
    for input in (identity, x -> view(x, :, :))
        @test_throws DimensionMismatch BandedOneBody(input(zeros(0, 4)))
        @test_throws DimensionMismatch BandedOneBody(input(zeros(3, 2)))
        @test_throws ArgumentError BandedOneBody(input(fill(NaN, 1, 4)))
        @test_throws ArgumentError BandedOneBody(input(fill(Inf, 1, 4)))
        @test_throws ArgumentError BandedOneBody(input([1.0 2.0; -0.5 0.25]))

        data = [1.0 2.0 3.0; -0.5 0.25 0.0; -0.125 0.0 -0.0]
        original = copy(data)
        one_body = BandedOneBody(input(data))
        fill!(data, 99.0)
        @test one_body.bands == original
        @test HFDMRG.h1_entry(one_body, 1, 3) == -0.125
        @test HFDMRG.h1_entry(one_body, 3, 3) == 3.0
    end
end
