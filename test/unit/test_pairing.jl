using Test
using EntrapmentAnalyses
using DataFrames

@testset "Pairing Functions" begin
    @testset "compute_pairing_vectors with vectors" begin
        # Create test data
        lib_sequences = ["PEPTIDE", "EPPTIDE", "SEQUENCE", "SEQEUNCE", "ANOTHER", "ANOHTER"]
        lib_charges = [2, 2, 3, 3, 2, 2]
        lib_entrap_groups = [0, 1, 0, 1, 0, 1]  # 0 = original, 1 = entrapment
        lib_pair_ids = [1, 1, 2, 2, 3, 3]
        
        results_sequences = ["PEPTIDE", "EPPTIDE", "SEQUENCE", "ANOTHER"]
        results_charges = [2, 2, 3, 2]
        
        pairing = compute_pairing_vectors(
            lib_sequences, lib_charges, lib_entrap_groups, lib_pair_ids,
            results_sequences, results_charges;
            show_progress = false
        )
        
        @test pairing.is_original == [true, false, true, true]
        @test pairing.pair_indices == [1, 1, 2, 3]
        @test pairing.entrap_labels == [0, 1, 0, 0]
        @test pairing.complement_indices == [2, 1, -1, -1]  # SEQUENCE has no pair in results
    end
    
    @testset "compute_pairing_vectors with DataFrames" begin
        # Create test library
        library_df = DataFrame(
            PeptideSequence = ["PEPTIDE", "EPPTIDE", "SEQUENCE", "SEQEUNCE"],
            PrecursorCharge = [2, 2, 3, 3],
            EntrapmentGroupId = [0, 1, 0, 1],
            PrecursorIdx = [1, 1, 2, 2]
        )
        
        results_df = DataFrame(
            stripped_seq = ["PEPTIDE", "EPPTIDE"],
            z = [2, 2]
        )
        
        pairing = compute_pairing_vectors(library_df, results_df; show_progress=false)
        
        @test pairing.is_original == [true, false]
        @test pairing.pair_indices == [1, 1]
        @test pairing.complement_indices == [2, 1]
    end
    
    @testset "get_complement_scores" begin
        scores = [0.9, 0.8, 0.7, 0.6, 0.5]
        complement_indices = [2, 1, 5, -1, 3]
        
        comp_scores = EntrapmentAnalyses.get_complement_scores(scores, complement_indices)
        
        @test comp_scores[1] ≈ 0.8  # scores[2]
        @test comp_scores[2] ≈ 0.9  # scores[1]
        @test comp_scores[3] ≈ 0.5  # scores[5]
        @test comp_scores[4] == -1.0  # No complement
        @test comp_scores[5] ≈ 0.7  # scores[3]
    end
    
    @testset "Error handling" begin
        # Mismatched lengths
        @test_throws ErrorException compute_pairing_vectors(
            ["A", "B"], [1], [0, 1], [1, 1], ["A"], [1];
            show_progress = false
        )
        
        # Sequence not in library
        @test_throws ErrorException compute_pairing_vectors(
            ["A"], [1], [0], [1], ["B"], [1];
            show_progress = false
        )
    end
end