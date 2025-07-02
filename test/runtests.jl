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
        include("unit/core/test_efdr_methods.jl")
        include("unit/core/test_scoring.jl")
        
        # Analysis tests
        include("unit/analysis/test_qvalue_calculation.jl")
        include("unit/analysis/test_protein_analysis.jl")
        
        # Plotting tests
        include("unit/plotting/test_visualization.jl")
    end
    
    # Integration tests
    @testset "Integration Tests" begin
        include("integration/test_api.jl")
    end
end