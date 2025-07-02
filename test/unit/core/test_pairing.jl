using Test
using DataFrames
using EntrapmentAnalyses

@testset "Pairing Tests" begin
    
    @testset "PeptideKey struct" begin
        # Test construction and equality
        key1 = EntrapmentAnalyses.PeptideKey("PEPTIDE", UInt8(2))
        key2 = EntrapmentAnalyses.PeptideKey("PEPTIDE", UInt8(2))
        key3 = EntrapmentAnalyses.PeptideKey("PEPTIDE", UInt8(3))
        key4 = EntrapmentAnalyses.PeptideKey("DIFFERENT", UInt8(2))
        
        @test key1 == key2
        @test key1 != key3
        @test key1 != key4
        
        # Test hashing
        @test hash(key1) == hash(key2)
        @test hash(key1) != hash(key3)
    end
    
    @testset "PlexPairKey struct" begin
        # Test construction and equality
        key1 = EntrapmentAnalyses.PlexPairKey(UInt8(1), 100)
        key2 = EntrapmentAnalyses.PlexPairKey(UInt8(1), 100)
        key3 = EntrapmentAnalyses.PlexPairKey(UInt8(2), 100)
        key4 = EntrapmentAnalyses.PlexPairKey(UInt8(1), 200)
        
        @test key1 == key2
        @test key1 != key3
        @test key1 != key4
        
        # Test hashing
        @test hash(key1) == hash(key2)
        @test hash(key1) != hash(key3)
    end
    
    @testset "ScorePair struct" begin
        pair = EntrapmentAnalyses.ScorePair(0.9f0, 0.8f0)
        @test pair.original_score == 0.9f0
        @test pair.entrapment_score == 0.8f0
    end
    
    @testset "init_entrapment_pairs_dict" begin
        # Create test library
        lib_df = DataFrame(
            PeptideSequence = ["PEPTIDE1", "PEPTIDE1", "PEPTIDE2", "PEPTIDE2"],
            PrecursorCharge = [2, 2, 3, 3],
            EntrapmentGroupId = [0, 1, 0, 1],  # First is original, second is entrapment
            PrecursorIdx = [100, 100, 200, 200]  # Same pair ID for pairs
        )
        
        pair_dict, is_original_dict = EntrapmentAnalyses.init_entrapment_pairs_dict(lib_df)
        
        # Test pair_dict
        key1 = EntrapmentAnalyses.PeptideKey("PEPTIDE1", UInt8(2))
        key2 = EntrapmentAnalyses.PeptideKey("PEPTIDE2", UInt8(3))
        
        @test haskey(pair_dict, key1)
        @test haskey(pair_dict, key2)
        @test pair_dict[key1] == 100
        @test pair_dict[key2] == 200
        
        # Test is_original_dict
        @test is_original_dict[key1] == true  # First occurrence is original
        @test is_original_dict[key2] == true  # First occurrence is original
    end
    
    @testset "add_plex_complement_scores!" begin
        # Create test library
        lib_df = DataFrame(
            PeptideSequence = ["PEPTIDE1", "EDITPEP1", "PEPTIDE2", "EDITPEP2"],
            PrecursorCharge = [2, 2, 3, 3],
            EntrapmentGroupId = [0, 1, 0, 1],
            PrecursorIdx = [100, 100, 200, 200]
        )
        
        # Create test results
        results_df = DataFrame(
            stripped_seq = ["PEPTIDE1", "EDITPEP1", "PEPTIDE2", "EDITPEP2"],
            z = [2, 2, 3, 3],
            PredVal = Float32[0.9, 0.8, 0.7, 0.6],
            channel = [0, 0, 0, 0],  # Same plex
            file_name = ["file1", "file1", "file1", "file1"]
        )
        
        # Initialize dictionaries
        pair_dict, is_original_dict = EntrapmentAnalyses.init_entrapment_pairs_dict(lib_df)
        
        # Add complement scores
        EntrapmentAnalyses.add_plex_complement_scores!(
            results_df, pair_dict, is_original_dict;
            show_progress = false
        )
        
        # Check columns were added
        @test hasproperty(results_df, :complement_score)
        @test hasproperty(results_df, :is_original)
        @test hasproperty(results_df, :pair_id)
        
        # Check values
        @test results_df.is_original == [true, false, true, false]
        @test results_df.pair_id == [100, 100, 200, 200]
        @test results_df.complement_score ≈ [0.8f0, 0.9f0, 0.6f0, 0.7f0]
    end
    
    @testset "add_plex_complement_scores! - multi-plex" begin
        # Test that same peptide in different plexes can have different complement scores
        lib_df = DataFrame(
            PeptideSequence = ["PEPTIDE1", "EDITPEP1"],
            PrecursorCharge = [2, 2],
            EntrapmentGroupId = [0, 1],
            PrecursorIdx = [100, 100]
        )
        
        # Same peptide pair in two different plexes
        results_df = DataFrame(
            stripped_seq = ["PEPTIDE1", "EDITPEP1", "PEPTIDE1", "EDITPEP1"],
            z = [2, 2, 2, 2],
            PredVal = Float32[0.9, 0.8, 0.95, 0.85],  # Different scores in different plexes
            channel = [0, 0, 1, 1],  # Two different plexes
            file_name = ["file1", "file1", "file1", "file1"]
        )
        
        pair_dict, is_original_dict = EntrapmentAnalyses.init_entrapment_pairs_dict(lib_df)
        
        EntrapmentAnalyses.add_plex_complement_scores!(
            results_df, pair_dict, is_original_dict;
            show_progress = false
        )
        
        # Check plex-specific complement scores
        # Plex 0: original=0.9, entrapment=0.8
        # Plex 1: original=0.95, entrapment=0.85
        @test results_df.complement_score[1] ≈ 0.8f0   # Original in plex 0 gets entrapment score
        @test results_df.complement_score[2] ≈ 0.9f0   # Entrapment in plex 0 gets original score
        @test results_df.complement_score[3] ≈ 0.85f0  # Original in plex 1 gets entrapment score
        @test results_df.complement_score[4] ≈ 0.95f0  # Entrapment in plex 1 gets original score
    end
    
    @testset "compute_pairing_vectors!" begin
        # Create test library
        lib_df = DataFrame(
            PeptideSequence = ["PEPTIDE1", "EDITPEP1", "PEPTIDE2", "EDITPEP2"],
            PrecursorCharge = [2, 2, 3, 3],
            EntrapmentGroupId = [0, 1, 0, 1],
            PrecursorIdx = [100, 100, 200, 200]
        )
        
        # Create test results
        results_df = DataFrame(
            stripped_seq = ["PEPTIDE1", "EDITPEP1", "PEPTIDE2", "EDITPEP2"],
            z = [2, 2, 3, 3],
            PredVal = Float32[0.9, 0.8, 0.7, 0.6],
            channel = [0, 0, 0, 0],
            file_name = ["file1", "file1", "file1", "file1"]
        )
        
        # Run the main function
        EntrapmentAnalyses.compute_pairing_vectors!(
            lib_df, results_df;
            show_progress = false
        )
        
        # Check all columns were added
        @test hasproperty(results_df, :is_original)
        @test hasproperty(results_df, :pair_id)
        @test hasproperty(results_df, :entrap_label)
        @test hasproperty(results_df, :complement_score)
        @test hasproperty(results_df, :complement_indices)
        
        # Check values
        @test results_df.is_original == [true, false, true, false]
        @test results_df.entrap_label == [0, 1, 0, 1]
        @test results_df.complement_score ≈ [0.8f0, 0.9f0, 0.6f0, 0.7f0]
        @test all(results_df.complement_indices .== -1)  # Deprecated column
    end
    
    @testset "compute_pairing_vectors! - missing peptide" begin
        # Library missing a peptide that's in results
        lib_df = DataFrame(
            PeptideSequence = ["PEPTIDE1", "PEPTIDE1"],
            PrecursorCharge = [2, 2],
            EntrapmentGroupId = [0, 1],
            PrecursorIdx = [100, 100]
        )
        
        results_df = DataFrame(
            stripped_seq = ["PEPTIDE1", "PEPTIDE2"],  # PEPTIDE2 not in library
            z = [2, 3],
            PredVal = Float32[0.9, 0.8],
            channel = [0, 0],
            file_name = ["file1", "file1"]
        )
        
        # Should throw error for missing peptide
        @test_throws ErrorException EntrapmentAnalyses.compute_pairing_vectors!(
            lib_df, results_df;
            show_progress = false
        )
    end
    
    @testset "get_complement_scores (deprecated)" begin
        scores = Float32[0.9, 0.8, 0.7, 0.6]
        complement_indices = [2, 1, 4, 3]
        
        # Test deprecated function still works
        @test_logs (:warn,) match_mode=:any begin
            comp_scores = EntrapmentAnalyses.get_complement_scores(scores, complement_indices)
            @test comp_scores ≈ Float32[0.8, 0.9, 0.6, 0.7]
        end
        
        # Test with invalid indices
        complement_indices_invalid = [2, -1, 5, 3]  # -1 and 5 are invalid
        @test_logs (:warn,) match_mode=:any begin
            comp_scores = EntrapmentAnalyses.get_complement_scores(scores, complement_indices_invalid)
            @test comp_scores[2] == -1.0f0  # No complement
            @test comp_scores[3] == -1.0f0  # Out of bounds
        end
    end
end