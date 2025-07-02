using Test
using DataFrames
using EntrapmentAnalyses
using CSV

@testset "Data Loaders Tests" begin
    
    @testset "handle_missing_values!" begin
        # Test with missing values
        @testset "Replace missing values" begin
            df = DataFrame(
                col1 = [1, missing, 3, missing, 5],
                col2 = ["a", missing, "c", "d", missing]
            )
            
            # Test numeric column
            n_missing = EntrapmentAnalyses.handle_missing_values!(df, :col1, 0, Int, "0")
            @test n_missing == 2
            @test all(df.col1 .== [1, 0, 3, 0, 5])
            @test eltype(df.col1) == Int
            
            # Test string column
            n_missing = EntrapmentAnalyses.handle_missing_values!(df, :col2, "", String, "empty string")
            @test n_missing == 2
            @test all(df.col2 .== ["a", "", "c", "d", ""])
            @test eltype(df.col2) == String
        end
        
        # Test with no missing values
        @testset "No missing values" begin
            df = DataFrame(col1 = [1, 2, 3, 4, 5])
            n_missing = EntrapmentAnalyses.handle_missing_values!(df, :col1, 0, Int, "0")
            @test n_missing == 0
            @test all(df.col1 .== [1, 2, 3, 4, 5])
        end
        
        # Test type conversion
        @testset "Type conversion" begin
            df = DataFrame(col1 = [1.0, 2.0, 3.0])
            EntrapmentAnalyses.handle_missing_values!(df, :col1, 0, Int, "0")
            @test eltype(df.col1) == Int
            @test all(df.col1 .== [1, 2, 3])
        end
        
        # Test error on missing column
        @testset "Missing column error" begin
            df = DataFrame(col1 = [1, 2, 3])
            @test_throws ErrorException EntrapmentAnalyses.handle_missing_values!(df, :col2, 0, Int, "0")
        end
    end
    
    @testset "load_parquet" begin
        # Create a test DataFrame
        test_df = DataFrame(
            stripped_seq = ["PEPTIDE1", "PEPTIDE2", "PEPTIDE3"],
            z = UInt8[2, 3, 2],
            PredVal = Float32[0.9, 0.8, 0.7],
            decoy = [false, true, false],
            file_name = ["file1", "file1", "file1"],
            protein = ["PROT1", "PROT2", "PROT3"]
        )
        
        # Save as temporary parquet file
        temp_file = tempname() * ".parquet"
        
        # For now, we'll skip this test as it requires Parquet writing capability
        @test_skip "Parquet file loading"
        
        # Test non-existent file error
        @test_throws ErrorException EntrapmentAnalyses.load_parquet("nonexistent.parquet")
    end
    
    @testset "load_spectral_library" begin
        # Create test library data
        test_lib = DataFrame(
            PeptideSequence = ["PEPTIDE1", "PEPTIDE2", "PEPTIDE3", "PEPTIDE4"],
            PrecursorCharge = UInt8[2, 2, 3, 3],
            EntrapmentGroupId = [0, 1, 0, 1],
            PrecursorIdx = [1, 1, 2, 2]
        )
        
        # Save as temporary TSV file
        temp_lib = tempname() * ".tsv"
        CSV.write(temp_lib, test_lib; delim='\t')
        
        # Test loading
        @testset "Valid library loading" begin
            lib_df = EntrapmentAnalyses.load_spectral_library(temp_lib)
            @test nrow(lib_df) == 4
            @test all(propertynames(lib_df) .== [:PeptideSequence, :PrecursorCharge, :EntrapmentGroupId, :PrecursorIdx])
            @test sum(lib_df.EntrapmentGroupId .== 0) == 2  # 2 targets
            @test sum(lib_df.EntrapmentGroupId .> 0) == 2   # 2 entrapments
        end
        
        # Clean up
        rm(temp_lib)
        
        # Test missing file error
        @test_throws ErrorException EntrapmentAnalyses.load_spectral_library("nonexistent.tsv")
        
        # Test missing columns error
        @testset "Missing columns error" begin
            incomplete_lib = DataFrame(
                PeptideSequence = ["PEPTIDE1", "PEPTIDE2"],
                PrecursorCharge = [2, 3]
                # Missing EntrapmentGroupId and PrecursorIdx
            )
            temp_incomplete = tempname() * ".tsv"
            CSV.write(temp_incomplete, incomplete_lib; delim='\t')
            
            @test_throws ErrorException EntrapmentAnalyses.load_spectral_library(temp_incomplete)
            
            rm(temp_incomplete)
        end
    end
    
    @testset "load_spectral_library with missing values" begin
        # Create library with missing values
        test_lib = DataFrame(
            PeptideSequence = ["PEPTIDE1", missing, "PEPTIDE3"],
            PrecursorCharge = [2, missing, 3],
            EntrapmentGroupId = [0, 1, missing],
            PrecursorIdx = [1, missing, 2]
        )
        
        temp_lib = tempname() * ".tsv"
        CSV.write(temp_lib, test_lib; delim='\t')
        
        # Test that missing values are handled
        lib_df = EntrapmentAnalyses.load_spectral_library(temp_lib)
        
        @test !any(ismissing, lib_df.PeptideSequence)
        @test !any(ismissing, lib_df.PrecursorCharge)
        @test !any(ismissing, lib_df.EntrapmentGroupId)
        @test !any(ismissing, lib_df.PrecursorIdx)
        
        # Check default values
        @test lib_df.PeptideSequence[2] == ""
        @test lib_df.PrecursorCharge[2] == 0
        @test lib_df.EntrapmentGroupId[3] == 0
        @test lib_df.PrecursorIdx[2] == 0
        
        rm(temp_lib)
    end
end