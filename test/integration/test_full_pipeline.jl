using Test
using EntrapmentAnalyses
using DataFrames
using CSV

@testset "Full Pipeline Integration" begin
    @testset "Create test data" begin
        # Create a small test dataset
        
        # Create test library
        library_df = DataFrame(
            PeptideSequence = repeat(["PEPTIDEA", "EPTIDPEA", "PEPTIDEB", "EPTIDPEB", 
                                     "PEPTIDEC", "EPTIDPEC"], 2),
            PrecursorCharge = repeat([2, 2, 2, 2, 3, 3], 2),
            EntrapmentGroupId = repeat([0, 1, 0, 1, 0, 1], 2),
            PrecursorIdx = repeat([1, 1, 2, 2, 3, 3], 2)
        )
        
        # Create test results
        results_df = DataFrame(
            stripped_seq = ["PEPTIDEA", "EPTIDPEA", "PEPTIDEB", "EPTIDPEB", 
                           "PEPTIDEC", "EPTIDPEC", "PEPTIDEA", "PEPTIDEB"],
            z = [2, 2, 2, 2, 3, 3, 2, 2],
            PredVal = Float32[0.95, 0.85, 0.90, 0.70, 0.88, 0.92, 0.80, 0.75],
            decoy = [false, false, false, false, false, false, true, true],
            channel = UInt8[1, 1, 1, 1, 1, 1, 1, 1],
            file_name = fill("test_run1", 8),
            protein = ["ProtA", "ProtA", "ProtB", "ProtB", "ProtC", "ProtC", "ProtA", "ProtB"]
        )
        
        # Save test files
        test_dir = mktempdir()
        library_path = joinpath(test_dir, "test_library.tsv")
        results_path = joinpath(test_dir, "test_results.arrow")
        
        CSV.write(library_path, library_df, delim='\t')
        
        # For this test, we'll work directly with the DataFrame
        # In real usage, this would be saved as Parquet/Arrow
        
        @test isfile(library_path)
        @test nrow(library_df) == 12
        @test nrow(results_df) == 8
        
        # Test pairing
        pairing_info = compute_pairing_vectors(library_df, results_df; show_progress=false)
        
        @test length(pairing_info.is_original) == nrow(results_df)
        @test sum(pairing_info.is_original) == 4  # 4 original sequences
        @test all(pairing_info.complement_indices[1:6] .> 0)  # All have pairs
        
        # Clean up
        rm(test_dir, recursive=true)
    end
    
    @testset "Test vector-based calculations" begin
        # Test with simple known data
        scores = Float64[0.9, 0.5, 0.8, 0.4]
        is_decoy = [false, false, true, false]
        
        qvals = EntrapmentAnalyses.calculate_qvalues(scores, is_decoy)
        
        # At score 0.9: 0 decoys, 1 target → 0.0
        # At score 0.8: 1 decoy, 1 target → 1.0
        # At score 0.5: 1 decoy, 2 targets → 0.5
        # At score 0.4: 1 decoy, 3 targets → 0.333...
        
        @test qvals[1] ≈ 0.0
        @test qvals[3] ≈ 1.0
        @test qvals[2] ≈ 0.5
        @test qvals[4] ≈ 1/3
        
        # After monotonization, values should be non-decreasing in score order
        sorted_indices = sortperm(scores, rev=true)
        @test issorted(qvals[sorted_indices])
    end
end