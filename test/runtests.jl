using Test
using EntrapmentAnalyses

@testset "EntrapmentAnalyses.jl Tests" begin
    # Dummy test to verify testing infrastructure
    @testset "Basic Setup" begin
        @test true
        @test 1 + 1 == 2
        
        # Verify package can be loaded
        @test isdefined(Main, :EntrapmentAnalyses)
    end
end