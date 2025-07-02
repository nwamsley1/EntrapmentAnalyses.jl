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
    
    # Unit tests
    @testset "Unit Tests" begin
        # IO tests
        include("unit/io/test_data_loaders.jl")
        
        # Core tests
        include("unit/core/test_pairing.jl")
    end
end