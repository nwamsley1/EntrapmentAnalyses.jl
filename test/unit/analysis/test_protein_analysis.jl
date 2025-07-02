using Test
using DataFrames
using EntrapmentAnalyses

@testset "Protein Analysis Tests" begin
    
    @testset "add_entrapment_group_column!" begin
        # Create test data
        @testset "Basic functionality" begin
            results_df = DataFrame(
                stripped_seq = ["PEPTIDE1", "PEPTIDE2", "PEPTIDE3", "PEPTIDE4"],
                score = [0.9, 0.8, 0.7, 0.6]
            )
            
            target_seqs = Set(["PEPTIDE1", "PEPTIDE3"])
            entrap_seqs = Set(["PEPTIDE2", "PEPTIDE4"])
            
            EntrapmentAnalyses.add_entrapment_group_column!(results_df, target_seqs, entrap_seqs)
            
            @test hasproperty(results_df, :entrapment_group)
            @test results_df.entrapment_group == [false, true, false, true]
        end
        
        # Test with overlapping sets (should not happen in practice)
        @testset "Overlapping sets" begin
            results_df = DataFrame(
                stripped_seq = ["PEPTIDE1", "PEPTIDE2"]
            )
            
            target_seqs = Set(["PEPTIDE1", "PEPTIDE2"])
            entrap_seqs = Set(["PEPTIDE2"])  # PEPTIDE2 in both sets
            
            EntrapmentAnalyses.add_entrapment_group_column!(results_df, target_seqs, entrap_seqs)
            
            # Should prioritize entrapment
            @test results_df.entrapment_group == [false, true]
        end
        
        # Test with empty DataFrame
        @testset "Empty DataFrame" begin
            results_df = DataFrame(stripped_seq = String[])
            target_seqs = Set{String}()
            entrap_seqs = Set{String}()
            
            EntrapmentAnalyses.add_entrapment_group_column!(results_df, target_seqs, entrap_seqs)
            
            @test hasproperty(results_df, :entrapment_group)
            @test isempty(results_df.entrapment_group)
        end
    end
    
    @testset "calculate_protein_qvalues_per_run!" begin
        # Test basic per-run q-value calculation
        @testset "Multiple runs" begin
            protein_df = DataFrame(
                file_name = ["run1", "run1", "run1", "run2", "run2", "run2"],
                decoy = [false, true, false, false, false, true],
                PredVal = Float32[0.9, 0.8, 0.7, 0.85, 0.75, 0.65]
            )
            
            # Sort by score within runs (required for q-value calculation)
            sort!(protein_df, [:file_name, :PredVal], rev=[false, true])
            
            EntrapmentAnalyses.calculate_protein_qvalues_per_run!(protein_df)
            
            @test hasproperty(protein_df, :Protein_Qvalue)
            @test all(0 .<= protein_df.Protein_Qvalue .<= 1)
            
            # Check that q-values are calculated per run
            run1_df = protein_df[protein_df.file_name .== "run1", :]
            run2_df = protein_df[protein_df.file_name .== "run2", :]
            
            @test run1_df.Protein_Qvalue ≈ [0.0, 1/2, 1/2]
            @test run2_df.Protein_Qvalue ≈ [0.0, 0, 1/2]
            # Each run should have its own q-value calculation
            @test maximum(run1_df.Protein_Qvalue) ≤ 1.0
            @test maximum(run2_df.Protein_Qvalue) ≤ 1.0
        end
        
        # Test single run
        @testset "Single run" begin
            protein_df = DataFrame(
                file_name = ["run1", "run1", "run1"],
                decoy = [false, true, false],
                PredVal = Float32[0.9, 0.8, 0.7]
            )
            
            sort!(protein_df, :PredVal, rev=true)
            
            EntrapmentAnalyses.calculate_protein_qvalues_per_run!(protein_df)
            
            @test hasproperty(protein_df, :Protein_Qvalue)
            @test protein_df.Protein_Qvalue ≈ [0.0, 1/2, 1/2]
        end
        
        # Test all decoys
        @testset "All decoys" begin
            protein_df = DataFrame(
                file_name = ["run1", "run1"],
                decoy = [true, true]
            )
            
            EntrapmentAnalyses.calculate_protein_qvalues_per_run!(protein_df)
            
            # Should handle division by zero gracefully
            @test hasproperty(protein_df, :Protein_Qvalue)
            @test all(protein_df.Protein_Qvalue .≈ 1)
        end
    end
    
    @testset "prepare_protein_analysis" begin
        # Create test library
        library_df = DataFrame(
            PeptideSequence = ["PEPTIDE1", "PEPTIDE2", "PEPTIDE3", "PEPTIDE4"],
            EntrapmentGroupId = [0, 1, 0, 1]  # 1,3 are targets; 2,4 are entrapments
        )
        
        # Create test results
        results_df = DataFrame(
            stripped_seq = ["PEPTIDE1", "PEPTIDE2", "PEPTIDE1", "PEPTIDE3", "PEPTIDE4"],
            PredVal = Float32[0.9, 0.8, 0.85, 0.7, 0.6],
            protein = ["PROT1", "PROT2", "PROT1", "PROT3", "PROT4"],
            file_name = ["run1", "run1", "run1", "run1", "run1"],
            channel = UInt8[0, 0, 0, 0, 0],
            decoy = [false, false, false, false, false]
        )
        
        protein_df = EntrapmentAnalyses.prepare_protein_analysis(results_df, library_df)
        
        @test hasproperty(protein_df, :entrapment_group)
        @test hasproperty(protein_df, :Protein_Qvalue)
        
        # Should have best PSM per protein group
        # PROT1 appears twice with scores 0.9 and 0.85, should keep 0.9
        prot1_rows = protein_df[protein_df.protein .== "PROT1", :]
        @test nrow(prot1_rows) == 1
        @test prot1_rows.PredVal[1] ≈ 0.9f0
        
        # Check that entrapment groups are correctly assigned
        @info display(results_df)
        @test results_df[results_df[!,:stripped_seq].=="PEPTIDE2",:entrapment_group] == [true]
        @test results_df[results_df[!,:stripped_seq].=="PEPTIDE1",:entrapment_group] == [false,false]
        @test results_df[results_df[!,:stripped_seq].=="PEPTIDE3",:entrapment_group] == [false]
        @test results_df[results_df[!,:stripped_seq].=="PEPTIDE4",:entrapment_group] == [true]
    end
    
    @testset "rollup_to_protein_groups (vector-based)" begin
        # Test basic rollup
        @testset "Basic rollup" begin
            scores = Float32[0.9, 0.8, 0.85, 0.7, 0.6]
            proteins = ["PROT1", "PROT2", "PROT1", "PROT3", "PROT3"]
            file_names = ["run1", "run1", "run1", "run1", "run1"]
            channels = [0, 0, 0, 0, 0]
            is_decoy = [false, false, false, false, true]
            is_entrapment = [false, true, false, false, true]
            
            result = EntrapmentAnalyses.rollup_to_protein_groups(
                scores, proteins, file_names, channels, is_decoy, is_entrapment
            )
            
            # Should have 4 unique groups:
            # - PROT1, run1, channel=0, decoy=false, entrap=false (best: 0.9)
            # - PROT2, run1, channel=0, decoy=false, entrap=true (0.8)
            # - PROT3, run1, channel=0, decoy=false, entrap=false (0.7)
            # - PROT3, run1, channel=0, decoy=true, entrap=true (0.6)
            @test length(result.protein_scores) == 4
            
            # Find PROT1 entry
            prot1_idx = findfirst(result.protein_names .== "PROT1")
            @test !isnothing(prot1_idx)
            @test result.protein_scores[prot1_idx] ≈ 0.9f0
        end
        
        # Test with multiple files and channels
        @testset "Multiple files and channels" begin
            scores = Float32[0.9, 0.8, 0.85]
            proteins = ["PROT1", "PROT1", "PROT1"]
            file_names = ["run1", "run2", "run1"]
            channels = [0, 0, 1]
            is_decoy = [false, false, false]
            is_entrapment = [false, false, false]
            
            result = EntrapmentAnalyses.rollup_to_protein_groups(
                scores, proteins, file_names, channels, is_decoy, is_entrapment
            )
            
            # Should have 3 groups (different file/channel combinations)
            @test length(result.protein_scores) == 3
        end
        
        # Test input validation
        @testset "Input validation" begin
            scores = Float32[0.9, 0.8]
            proteins = ["PROT1"]  # Wrong length
            file_names = ["run1", "run1"]
            channels = [0, 0]
            is_decoy = [false, false]
            is_entrapment = [false, false]
            
            @test_throws ErrorException EntrapmentAnalyses.rollup_to_protein_groups(
                scores, proteins, file_names, channels, is_decoy, is_entrapment
            )
        end
    end
    
    @testset "rollup_to_protein_groups (DataFrame)" begin
        df = DataFrame(
            PredVal = Float32[0.9, 0.8, 0.85, 0.7],
            protein = ["PROT1", "PROT2", "PROT1", "PROT2"],
            file_name = ["run1", "run1", "run1", "run1"],
            channel = UInt8[0, 0, 0, 0],
            decoy = [false, false, false, true],
            entrapment_group = [false, true, false, true]
        )
        
        result = EntrapmentAnalyses.rollup_to_protein_groups(df)
        
        # Should have 3 unique groups:
        # - PROT1, decoy=false, entrap=false (best: 0.9)
        # - PROT2, decoy=false, entrap=true (0.8)
        # - PROT2, decoy=true, entrap=true (0.7)
        @test length(result.protein_scores) == 3
        
        # Verify all expected fields are present
        @test haskey(result, :protein_scores)
        @test haskey(result, :protein_names)
        @test haskey(result, :protein_files)
        @test haskey(result, :protein_channels)
        @test haskey(result, :protein_is_decoy)
        @test haskey(result, :protein_is_entrapment)
        @test haskey(result, :original_indices)
    end
end