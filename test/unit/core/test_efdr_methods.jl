using Test
using DataFrames
using EntrapmentAnalyses

@testset "EFDR Methods Tests" begin
    
    @testset "calculate_combined_efdr" begin
        # Test basic functionality
        @testset "Basic calculation" begin
            scores = Float32[0.9, 0.8, 0.7, 0.6, 0.5]
            entrap_labels = [0, 1, 0, 1, 0]  # 3 originals, 2 entrapments
            qvals = Float32[0.01, 0.02, 0.03, 0.04, 0.05]
            
            efdr = EntrapmentAnalyses.calculate_combined_efdr(
                scores, entrap_labels, qvals;
                show_progress = false
            )
            
            @test length(efdr) == 5
            @test all(0 .<= efdr .<= 1)
            
            # Check specific values - as we go through sorted order:
            # i=1: score=0.9, entrap=0 -> Nτ=1, Nε=0 -> efdr=0
            # i=2: score=0.8, entrap=1 -> Nτ=1, Nε=1 -> efdr=(1*(1+1))/2=1
            # i=3: score=0.7, entrap=0 -> Nτ=2, Nε=1 -> efdr=(1*2)/3=0.667
            # etc.
            @test efdr[1] ≈ 0.0  # First item (highest score, lowest qval)
            @test efdr[2] ≈ 1 
            @test efdr[3] ≈ 2/3 
            @test efdr[4] ≈ 1 
            @test efdr[5] ≈ 4/5
        end
        
        # Test with all originals
        @testset "All originals" begin
            scores = Float32[0.9, 0.8, 0.7]
            entrap_labels = [0, 0, 0]  # All originals
            qvals = Float32[0.01, 0.02, 0.03]
            
            efdr = EntrapmentAnalyses.calculate_combined_efdr(
                scores, entrap_labels, qvals;
                show_progress = false
            )
            
            @test all(efdr .== 0.0)  # No entrapments means EFDR = 0
        end
        
        # Test with all entrapments
        @testset "All entrapments" begin
            scores = Float32[0.9, 0.8, 0.7]
            entrap_labels = [1, 1, 1]  # All entrapments
            qvals = Float32[0.01, 0.02, 0.03]
            
            efdr = EntrapmentAnalyses.calculate_combined_efdr(
                scores, entrap_labels, qvals;
                show_progress = false
            )
            
            @test all(efdr .== 1.0)  # All entrapments means EFDR = 1
        end
        
        # Test with different r value
        @testset "Different r value" begin
            scores = Float32[0.9, 0.8, 0.7, 0.6]
            entrap_labels = [0, 1, 0, 1]
            qvals = Float32[0.01, 0.02, 0.03, 0.04]
            
            efdr_r1 = EntrapmentAnalyses.calculate_combined_efdr(
                scores, entrap_labels, qvals;
                r = 1.0f0, show_progress = false
            )
            
            efdr_r2 = EntrapmentAnalyses.calculate_combined_efdr(
                scores, entrap_labels, qvals;
                r = 2.0f0, show_progress = false
            )
            
            # Different r values should give different results
            @test !all(efdr_r1 .≈ efdr_r2)
        end
        
        # Test input validation
        @testset "Input validation" begin
            scores = Float32[0.9, 0.8]
            entrap_labels = [0, 1, 0]  # Wrong length
            qvals = Float32[0.01, 0.02]
            
            @test_throws ErrorException EntrapmentAnalyses.calculate_combined_efdr(
                scores, entrap_labels, qvals;
                show_progress = false
            )
        end
    end
    
    @testset "calculate_paired_efdr" begin
        # Test basic functionality
        @testset "Basic calculation" begin
            # Original-entrapment pairs
            scores = Float32[0.9, 0.8, 0.7, 0.6]
            complement_scores = Float32[0.8, 0.9, 0.6, 0.7]  # Swapped for pairs
            is_original = [true, false, true, false]
            qvals = Float32[0.01, 0.02, 0.03, 0.04]
            
            efdr = EntrapmentAnalyses.calculate_paired_efdr(
                scores, complement_scores, is_original, qvals;
                show_progress = false
            )
            
            @test length(efdr) == 4
            @test all(0 .<= efdr .<= 1)
            @test efdr ≈ [0, 1/2, 1/3, 1/2]  # Expected EFDR values
        end
        
        # Test unpaired peptides
        @testset "Unpaired peptides" begin
            scores = Float32[0.9, 0.8, 0.7]
            complement_scores = Float32[-1.0, 0.9, -1.0]  # First and third have no complement
            is_original = [true, false, true]
            qvals = Float32[0.01, 0.02, 0.03]
            
            efdr = EntrapmentAnalyses.calculate_paired_efdr(
                scores, complement_scores, is_original, qvals;
                show_progress = false
            )
            
            @test length(efdr) == 3
            # Should handle unpaired peptides gracefully
        end
        
        # Test mutually exclusive conditions
        @testset "Mutually exclusive conditions" begin
            # Case 1: e_score >= s && s > o_score (Nεsτ case)
            scores = Float32[0.8, 0.7, 0.6]  # s = 0.7 for middle
            complement_scores = Float32[0.0, 0.0, 0.5]  # o_score = 0.5 for last
            is_original = [true, true, false]  # Last is entrapment
            qvals = Float32[0.01, 0.02, 0.03]
            
            efdr = EntrapmentAnalyses.calculate_paired_efdr(
                scores, complement_scores, is_original, qvals;
                show_progress = false
            )
            
            # For the entrapment: e_score=0.6, o_score=0.5
            # When s=0.7: neither condition is met (e_score < s)
            # When s=0.6: e_score >= s && s > o_score is true
            @test length(efdr) == 3
            @test efdr ≈ [0.0, 0.0, 2/3]  # First two originals, last is entrapment
        end
        
        @testset "Mutually exclusive conditions Nets and Nest" begin
            # Case 1: e_score >= s && s > o_score (Nεsτ case)
            scores = Float32[0.8, 0.7, 0.6, 0.5]  # s = 0.7 for middle
            complement_scores = Float32[0.0, 0.0, 0.5, 0.6]  # o_score = 0.5 for last
            is_original = [true, true, false, true]  # Last is entrapment
            qvals = Float32[0.01, 0.02, 0.03, 0.04]
            
            efdr = EntrapmentAnalyses.calculate_paired_efdr(
                scores, complement_scores, is_original, qvals;
                show_progress = false
            )
            
            # For the entrapment: e_score=0.6, o_score=0.5
            # When s=0.7: neither condition is met (e_score < s)
            # When s=0.6: e_score >= s && s > o_score is true
            @test length(efdr) == 4
            @test efdr ≈ [0.0, 0.0, 2/3, 3/4]  # First two originals, last is entrapment
        end

        # Test with all originals
        @testset "All originals" begin
            scores = Float32[0.9, 0.8, 0.7]
            complement_scores = Float32[-1.0, -1.0, -1.0]
            is_original = [true, true, true]
            qvals = Float32[0.01, 0.02, 0.03]
            
            efdr = EntrapmentAnalyses.calculate_paired_efdr(
                scores, complement_scores, is_original, qvals;
                show_progress = false
            )
            
            @test all(efdr .== 0.0)  # No entrapments
        end
        
        # Test input validation
        @testset "Input validation" begin
            scores = Float32[0.9, 0.8]
            complement_scores = Float32[0.8]  # Wrong length
            is_original = [true, false]
            qvals = Float32[0.01, 0.02]
            
            @test_throws ErrorException EntrapmentAnalyses.calculate_paired_efdr(
                scores, complement_scores, is_original, qvals;
                show_progress = false
            )
        end
    end
    
    @testset "EFDR Method structs" begin
        # Test CombinedEFDR
        @testset "CombinedEFDR struct" begin
            scores = [0.9, 0.8, 0.7]
            entrap_labels = [0, 1, 0]
            qvals = [0.01, 0.02, 0.03]
            
            method = EntrapmentAnalyses.CombinedEFDR(
                scores, entrap_labels, qvals, 1.0
            )
            
            efdr = EntrapmentAnalyses.calculate_efdr(method; show_progress=false)
            @test length(efdr) == 3
            @test all(0 .<= efdr .<= 1)
        end
        
        # Test PairedEFDR
        @testset "PairedEFDR struct" begin
            scores = [0.9, 0.8, 0.7, 0.6]
            complement_scores = [0.8, 0.9, 0.6, 0.7]
            is_original = [true, false, true, false]
            qvals = [0.01, 0.02, 0.03, 0.04]
            
            method = EntrapmentAnalyses.PairedEFDR(
                scores, complement_scores, is_original, qvals, 1.0
            )
            
            efdr = EntrapmentAnalyses.calculate_efdr(method; show_progress=false)
            @test length(efdr) == 4
            @test all(0 .<= efdr .<= 1)
        end
    end
    
    @testset "DataFrame calculate_paired_efdr" begin
        # Create test DataFrame with pairing info
        df = DataFrame(
            PredVal = Float32[0.9, 0.8, 0.7, 0.6],
            local_qvalue = Float32[0.01, 0.02, 0.03, 0.04]
        )
        
        pairing_info = (
            is_original = [true, false, true, false],
            complement_indices = [2, 1, 4, 3],
            pair_ids = [100, 100, 200, 200],
            entrap_labels = [0, 1, 0, 1]
        )
        
        # Test the DataFrame interface with deprecation warning
        @test_logs (:warn,) match_mode=:any begin
            efdr = EntrapmentAnalyses.calculate_paired_efdr(
                df;
                pairing_info = pairing_info,
                show_progress = false
            )
            
            @test length(efdr) == 4
            @test all(0 .<= efdr .<= 1)
        end
    end
end