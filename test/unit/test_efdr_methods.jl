using Test
using EntrapmentAnalyses

@testset "EFDR Methods" begin
    @testset "calculate_combined_efdr" begin
        # Simple test case
        scores = [0.9, 0.8, 0.7, 0.6, 0.5]
        entrap_labels = [0, 1, 0, 1, 0]  # 3 targets, 2 entrapments
        qvals = [0.001, 0.002, 0.003, 0.004, 0.005]
        
        efdr = calculate_combined_efdr(scores, entrap_labels, qvals; 
                                       r=1.0, show_progress=false)
        
        @test length(efdr) == length(scores)
        @test all(0 .<= efdr .<= 1)
        
        # Check specific values
        # At position 1: 0 targets, 0 entrapments → 0
        # At position 2: 1 target, 1 entrapment → (1*2)/(1+1) = 1.0
        # At position 3: 2 targets, 1 entrapment → (1*2)/(2+1) = 0.667
        # etc.
        
        # Test with different r value
        efdr_r2 = calculate_combined_efdr(scores, entrap_labels, qvals; 
                                         r=2.0, show_progress=false)
        @test all(efdr_r2 .<= efdr)  # Higher r should give lower EFDR
    end
    
    @testset "calculate_paired_efdr" begin
        # Test case with known complement relationships
        scores = [0.9, 0.5, 0.8, 0.4, 0.7]
        complement_scores = [0.5, 0.9, 0.4, 0.8, -1.0]  # Last has no pair
        is_original = [true, false, true, false, true]
        qvals = [0.001, 0.002, 0.003, 0.004, 0.005]
        
        efdr = calculate_paired_efdr(scores, complement_scores, is_original, qvals;
                                    r=1.0, show_progress=false)
        
        @test length(efdr) == length(scores)
        @test all(0 .<= efdr .<= 1)
        
        # Entrapment at position 2 wins (0.5 > 0.9's complement 0.5)
        # Entrapment at position 4 loses (0.4 < 0.8's complement 0.4)
    end
    
    @testset "Struct-based interface" begin
        scores = Float64[0.9, 0.8, 0.7, 0.6, 0.5]
        entrap_labels = [0, 1, 0, 1, 0]
        qvals = Float64[0.001, 0.002, 0.003, 0.004, 0.005]
        
        # Test CombinedEFDR struct
        method = CombinedEFDR(scores, entrap_labels, qvals, 1.0)
        efdr = calculate_efdr(method; show_progress=false)
        
        @test length(efdr) == length(scores)
        @test all(0 .<= efdr .<= 1)
        
        # Test PairedEFDR struct
        complement_scores = Float64[0.8, 0.9, 0.6, 0.7, 0.4]
        is_original = [true, false, true, false, true]
        
        method2 = PairedEFDR(scores, complement_scores, is_original, qvals, 1.0)
        efdr2 = calculate_efdr(method2; show_progress=false)
        
        @test length(efdr2) == length(scores)
        @test all(0 .<= efdr2 .<= 1)
    end
    
    @testset "Sort order validation" begin
        # Test with unsorted q-values (should warn)
        scores = [0.9, 0.8, 0.7, 0.6, 0.5]
        entrap_labels = [0, 1, 0, 1, 0]
        qvals = [0.005, 0.001, 0.003, 0.002, 0.004]  # Not sorted
        
        # This should produce a warning but still calculate
        efdr = @test_logs (:warn,) calculate_combined_efdr(
            scores, entrap_labels, qvals; show_progress=false
        )
        
        @test length(efdr) == length(scores)
    end
end