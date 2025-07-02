using Test
using EntrapmentAnalyses

@testset "Scoring Tests" begin
    
    @testset "monotonize! (in-place)" begin
        # Test basic monotonization
        @testset "Basic monotonization" begin
            values = Float32[0.1, 0.3, 0.2, 0.4, 0.5]
            result = EntrapmentAnalyses.monotonize!(values)
            
            @test result === values  # Should modify in-place
            @test values == Float32[0.1, 0.2, 0.2, 0.4, 0.5]
            @test issorted(values)  # Should be monotonic
        end
        
        # Test already monotonic
        @testset "Already monotonic" begin
            values = Float32[0.1, 0.2, 0.3, 0.4, 0.5]
            original = copy(values)
            EntrapmentAnalyses.monotonize!(values)
            
            @test values == original  # Should not change
        end
        
        # Test decreasing values
        @testset "Decreasing values" begin
            values = Float32[0.5, 0.4, 0.3, 0.2, 0.1]
            EntrapmentAnalyses.monotonize!(values)
            
            @test all(values .== 0.1f0)  # All should be set to minimum
        end
        
        # Test with duplicates
        @testset "With duplicates" begin
            values = Float32[0.1, 0.3, 0.3, 0.2, 0.4]
            EntrapmentAnalyses.monotonize!(values)
            
            @test values == Float32[0.1, 0.2, 0.2, 0.2, 0.4]
        end
        
        # Test single element
        @testset "Single element" begin
            values = Float32[0.5]
            EntrapmentAnalyses.monotonize!(values)
            
            @test values == Float32[0.5]
        end
        
        # Test empty vector
        @testset "Empty vector" begin
            values = Float32[]
            EntrapmentAnalyses.monotonize!(values)
            
            @test isempty(values)
        end
        
        # Test with Float64
        @testset "Float64 type" begin
            values = Float64[0.1, 0.3, 0.2, 0.4]
            EntrapmentAnalyses.monotonize!(values)
            
            @test values == Float64[0.1, 0.2, 0.2, 0.4]
            @test eltype(values) == Float64
        end
    end
    
    @testset "monotonize! with missing values" begin
        # Test with missing at end
        @testset "Missing at end" begin
            values = Union{Missing, Float32}[0.1, 0.3, 0.2, missing]
            EntrapmentAnalyses.monotonize!(values)
            
            @test values[1:3] == Float32[0.1, 0.2, 0.2]
            @test ismissing(values[4])
        end
        
        # Test with missing at beginning
        @testset "Missing at beginning" begin
            values = Union{Missing, Float32}[missing, 0.3, 0.2, 0.4]
            EntrapmentAnalyses.monotonize!(values)
            
            @test ismissing(values[1])
            @test values[2:4] == Float32[0.2, 0.2, 0.4]
        end
        
        # Test with multiple missing
        @testset "Multiple missing" begin
            values = Union{Missing, Float32}[0.1, missing, 0.3, missing, 0.2]
            EntrapmentAnalyses.monotonize!(values)
            
            @test values[1] == 0.1f0
            @test ismissing(values[2])
            @test values[3] == 0.2f0
            @test ismissing(values[4])
            @test values[5] == 0.2f0
        end
        
        # Test all missing
        @testset "All missing" begin
            values = Union{Missing, Float32}[missing, missing, missing]
            EntrapmentAnalyses.monotonize!(values)
            
            @test all(ismissing, values)
        end
    end
    
    @testset "monotonize (non-mutating)" begin
        # Test that original is not modified
        @testset "Original unchanged" begin
            original = Float32[0.1, 0.3, 0.2, 0.4]
            result = EntrapmentAnalyses.monotonize(original)
            
            @test original == Float32[0.1, 0.3, 0.2, 0.4]  # Original unchanged
            @test result == Float32[0.1, 0.2, 0.2, 0.4]    # Result is monotonic
            @test result !== original  # Different arrays
        end
        
        # Test with different types
        @testset "Different types" begin
            original = [0.1, 0.3, 0.2, 0.4]  # Float64
            result = EntrapmentAnalyses.monotonize(original)
            
            @test result == [0.1, 0.2, 0.2, 0.4]
            @test eltype(result) == Float64
        end
    end
    
    @testset "Edge cases and special values" begin
        # Test with values near 1.0
        @testset "Values near 1.0" begin
            values = Float32[0.8, 1.1, 0.9, 1.0]  # 1.1 should be capped at 1.0
            EntrapmentAnalyses.monotonize!(values)
            
            @test values == Float32[0.8, 0.9, 0.9, 1.0]
        end
        
        # Test with very small differences
        @testset "Small differences" begin
            values = Float32[0.10000, 0.10002, 0.10001, 0.10003]
            EntrapmentAnalyses.monotonize!(values)
            
            @test issorted(values)
            @test values[3] ≤ values[4]  # Monotonic
        end
    end
end