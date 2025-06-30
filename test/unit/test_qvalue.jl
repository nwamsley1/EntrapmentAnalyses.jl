using Test
using EntrapmentAnalyses
using DataFrames

@testset "Q-value Calculations" begin
    @testset "monotonize!" begin
        # Test basic monotonization
        values = Float32[0.01, 0.005, 0.02, 0.015, 0.03]
        monotonize!(values)
        @test issorted(values)
        @test values == Float32[0.005, 0.005, 0.015, 0.015, 0.03]
        
        # Test with missing values
        values_missing = Union{Missing, Float32}[0.01, missing, 0.02, 0.015, missing, 0.03]
        monotonize!(values_missing)
        @test values_missing[1] == 0.01
        @test ismissing(values_missing[2])
        @test values_missing[3] == 0.015  # Monotonized
        @test values_missing[4] == 0.015
        @test ismissing(values_missing[5])
        @test values_missing[6] == 0.03
        
        # Test non-mutating version
        original = Float32[0.01, 0.005, 0.02]
        result = EntrapmentAnalyses.monotonize(original)
        @test original == Float32[0.01, 0.005, 0.02]  # Original unchanged
        @test result == Float32[0.005, 0.005, 0.02]
    end
    
    @testset "calculate_qvalues with vectors" begin
        # Simple test case
        scores = Float32[0.9, 0.8, 0.7, 0.6, 0.5, 0.4]
        is_decoy = [false, true, false, false, true, false]
        
        qvals = EntrapmentAnalyses.calculate_qvalues(scores, is_decoy)
        
        @test length(qvals) == length(scores)
        @test all(0 .<= qvals .<= 1)
        
        # Check specific values
        # Position 1 (score 0.9): 0 decoys, 1 target → 0/1 = 0
        # Position 2 (score 0.8): 1 decoy, 1 target → 1/1 = 1
        # Position 3 (score 0.7): 1 decoy, 2 targets → 1/2 = 0.5
        # etc.
        
        # Verify monotonicity
        sorted_indices = sortperm(scores, rev=true)
        sorted_qvals = qvals[sorted_indices]
        @test issorted(sorted_qvals)
    end
    
    @testset "calculate_qvalues! with DataFrame" begin
        df = DataFrame(
            PredVal = Float32[0.9, 0.8, 0.7, 0.6, 0.5],
            decoy = [false, true, false, true, false]
        )
        
        calculate_qvalues!(df)
        
        @test hasproperty(df, :local_qvalue)
        @test all(0 .<= df.local_qvalue .<= 1)
        @test df.local_qvalue isa Vector{Float32}
        
        # Test error handling
        bad_df = DataFrame(score = [0.9, 0.8])  # Wrong column name
        @test_throws ErrorException calculate_qvalues!(bad_df)
    end
    
    @testset "calculate_global_qvalues!" begin
        # Create test data with multiple PSMs per precursor
        df = DataFrame(
            stripped_seq = ["PEPTIDE", "PEPTIDE", "SEQUENCE", "SEQUENCE", "ANOTHER"],
            z = [2, 2, 3, 3, 2],
            PredVal = Float32[0.9, 0.7, 0.8, 0.6, 0.5],
            decoy = [false, false, true, true, false],
            local_qvalue = Float32[0.1, 0.2, 0.3, 0.4, 0.5]
        )
        
        calculate_global_qvalues!(df)
        
        @test hasproperty(df, :global_qvalue)
        @test all(0 .<= df.global_qvalue .<= 1)
        
        # Check that same precursor gets same global q-value
        @test df.global_qvalue[1] == df.global_qvalue[2]  # Same PEPTIDE z=2
        @test df.global_qvalue[3] == df.global_qvalue[4]  # Same SEQUENCE z=3
    end
    
    @testset "Sort order validation" begin
        # Test validate_sort_order function
        @test EntrapmentAnalyses.validate_sort_order([1, 2, 3, 4], true)
        @test EntrapmentAnalyses.validate_sort_order([4, 3, 2, 1], false)
        @test !EntrapmentAnalyses.validate_sort_order([1, 3, 2, 4], true)
        @test !EntrapmentAnalyses.validate_sort_order([4, 2, 3, 1], false)
    end
end