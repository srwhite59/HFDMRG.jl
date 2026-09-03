@testset "independent hierarchical timing utility" begin
    timing = HFDMRG.TimeG
    timing.set_timing!(false)
    timing.set_timing_live!(false)
    timing.reset_timing_report!()
    @test !timing.timing_enabled()
    @test isempty(timing.current_timing_report().roots)

    timing.set_timing_thresholds!(expand=0.0, drop=0.0)
    timing.set_timing!(true)
    value = try
        HFDMRG.@timeg "outer" begin
            HFDMRG.@timeg "inner" sum(1:10)
        end
    finally
        timing.set_timing!(false)
    end
    @test value == 55
    report = timing.current_timing_report()
    @test length(report.roots) == 1
    @test report.roots[1].label == "outer"
    @test report.roots[1].call_count == 1
    @test only(report.roots[1].children).label == "inner"
    @test occursin("HFDMRG timing report", timing.timing_report(report))
    timing.reset_timing_report!()
    @test isempty(timing.current_timing_report().roots)
end
