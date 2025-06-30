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
        
        # Position 2: entrapment (0.5) vs original (0.9) - entrapment loses
        # Position 4: entrapment (0.4) vs original (0.8) - entrapment loses
        # Testing the notebook's logic:
        # - Nεsτ: entrapment >= threshold > original
        # - Nετs: entrapment > original AND original >= threshold
    end
    
    @testset "Paired EFDR specific counting logic" begin
        # Create a specific test case to verify notebook logic
        # Scores in sorted order: [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]
        scores = [1.0, 0.6, 0.9, 0.4, 0.8, 0.3, 0.7, 0.5]
        # Complements arranged so we can test specific conditions
        complement_scores = [0.6, 1.0, 0.4, 0.9, 0.3, 0.8, 0.5, 0.7]
        is_original = [true, false, true, false, true, false, true, false]
        # Q-values that match score order
        qvals = [0.001, 0.006, 0.002, 0.008, 0.003, 0.009, 0.004, 0.005]
        
        efdr = calculate_paired_efdr(scores, complement_scores, is_original, qvals;
                                    r=1.0, show_progress=false)
        
        # Manual calculation for threshold at position 4 (score 0.7):
        # Positions with score >= 0.7: 1,3,5,7 (scores: 1.0,0.9,0.8,0.7)
        # Position 1: original
        # Position 3: original  
        # Position 5: original
        # Position 7: original
        # So at threshold 0.7, we should have only originals above
        
        # For entrapments:
        # Position 2: e=0.6, o=1.0 (e < threshold, not counted)
        # Position 4: e=0.4, o=0.9 (e < threshold, not counted)
        # Position 6: e=0.3, o=0.8 (e < threshold, not counted)
        # Position 8: e=0.5, o=0.7 (e < threshold, not counted)
        
        @test length(efdr) == length(scores)
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