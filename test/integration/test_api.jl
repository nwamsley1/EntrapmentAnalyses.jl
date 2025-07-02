using Test
using DataFrames
using CSV
using EntrapmentAnalyses

@testset "API Integration Tests" begin
    
    # Helper function to create test data files
    function create_test_files()
        temp_dir = mktempdir()
        
        # Create test library file
        library_df = DataFrame(
            PeptideSequence = ["PEPTIDE1", "EDITPEP1", "PEPTIDE2", "EDITPEP2", "PEPTIDE3", "EDITPEP3"],
            PrecursorCharge = [2, 2, 3, 3, 2, 2],
            EntrapmentGroupId = [0, 1, 0, 1, 0, 1],
            PrecursorIdx = [100, 100, 200, 200, 300, 300]
        )
        library_path = joinpath(temp_dir, "test_library.tsv")
        CSV.write(library_path, library_df, delim='\t')
        
        # Create test results file (would be Parquet in real use)
        # For testing, we'll create a CSV and pretend it's loaded from Parquet
        results_df = DataFrame(
            stripped_seq = ["PEPTIDE1", "EDITPEP1", "PEPTIDE2", "EDITPEP2", "PEPTIDE3", "EDITPEP3",
                           "PEPTIDE1", "EDITPEP1", "PEPTIDE2"],
            z = UInt8[2, 2, 3, 3, 2, 2, 2, 2, 3],
            PredVal = Float32[0.9, 0.8, 0.7, 0.6, 0.85, 0.75, 0.88, 0.78, 0.68],
            decoy = [false, false, false, false, false, false, true, true, true],
            channel = UInt8[0, 0, 0, 0, 0, 0, 0, 0, 0],
            file_name = ["run1", "run1", "run1", "run1", "run1", "run1", "run1", "run1", "run1"],
            protein = ["PROT1", "PROT1", "PROT2", "PROT2", "PROT3", "PROT3", "PROT1", "PROT1", "PROT2"]
        )
        
        return temp_dir, library_path, results_df
    end
    
    @testset "run_efdr_analysis" begin
        temp_dir, library_path, results_df = create_test_files()
        output_dir = joinpath(temp_dir, "efdr_output")
        
        # Mock the load_parquet_results function for testing
        EntrapmentAnalyses.load_parquet_results(files::Vector{String}) = results_df
        
        @testset "Basic precursor analysis" begin
            # Run analysis
            result_df = EntrapmentAnalyses.run_efdr_analysis(
                ["dummy.parquet"],  # Will use mocked data
                library_path;
                output_dir = output_dir,
                local_qval_threshold = 1.0,  # No filtering for test
                show_progress = false
            )
            
            # Check output DataFrame
            @test isa(result_df, DataFrame)
            @test hasproperty(result_df, :combined_entrapment_fdr)
            @test hasproperty(result_df, :precursor_entrapment_fdr)
            @test hasproperty(result_df, :is_original)
            @test hasproperty(result_df, :pair_id)
            @test hasproperty(result_df, :complement_score)
            
            # Check no decoys in results
            @test !any(result_df.decoy)
            
            # Check output files
            @test isfile(joinpath(output_dir, "precursor_entrapment_results.tsv"))
            @test isfile(joinpath(output_dir, "precursor_efdr_comparison.pdf"))
            @test isfile(joinpath(output_dir, "analysis_report.md"))
            
            # Check EFDR values are reasonable
            @test all(0 .<= result_df.combined_entrapment_fdr .<= 1)
            @test all(0 .<= result_df.precursor_entrapment_fdr .<= 1)
        end
        
        @testset "With filtering" begin
            # Run with q-value filtering
            result_df = EntrapmentAnalyses.run_efdr_analysis(
                ["dummy.parquet"],
                library_path;
                output_dir = joinpath(output_dir, "filtered"),
                local_qval_threshold = 0.05,
                show_progress = false
            )
            
            # Should have fewer results due to filtering
            @test nrow(result_df) < nrow(results_df[.!results_df.decoy, :])
        end
        
        @testset "Single file convenience" begin
            # Test single file interface
            result_df = EntrapmentAnalyses.run_efdr_analysis(
                "dummy.parquet",  # Single string instead of vector
                library_path;
                output_dir = joinpath(output_dir, "single"),
                show_progress = false
            )
            
            @test isa(result_df, DataFrame)
        end
        
        # Clean up
        rm(temp_dir, recursive=true, force=true)
    end
    
    @testset "run_protein_efdr_analysis" begin
        temp_dir, library_path, results_df = create_test_files()
        output_dir = joinpath(temp_dir, "protein_output")
        
        # Mock the load_parquet_results function for testing
        EntrapmentAnalyses.load_parquet_results(files::Vector{String}) = results_df
        
        @testset "Basic protein analysis" begin
            # Run analysis
            result_df = EntrapmentAnalyses.run_protein_efdr_analysis(
                ["dummy.parquet"],
                library_path;
                output_dir = output_dir,
                local_qval_threshold = 1.0,  # No filtering
                show_progress = false
            )
            
            # Check output DataFrame
            @test isa(result_df, DataFrame)
            @test hasproperty(result_df, :combined_protein_fdr)
            @test hasproperty(result_df, :protein_group_entrapment_fdr)
            @test hasproperty(result_df, :Protein_Qvalue)
            @test hasproperty(result_df, :protein)
            
            # Check protein rollup (should have fewer rows than PSMs)
            @test nrow(result_df) < nrow(results_df[.!results_df.decoy, :])
            
            # Check output files
            @test isfile(joinpath(output_dir, "protein_group_entrapment.tsv"))
            @test isfile(joinpath(output_dir, "protein_efdr_comparison.pdf"))
            @test isfile(joinpath(output_dir, "analysis_report.md"))
            
            # Check EFDR values
            @test all(0 .<= result_df.combined_protein_fdr .<= 1)
            @test all(0 .<= result_df.protein_group_entrapment_fdr .<= 1)
        end
        
        # Clean up
        rm(temp_dir, recursive=true, force=true)
    end
    
    @testset "Edge cases" begin
        temp_dir = mktempdir()
        
        @testset "Empty results" begin
            # Create library but empty results
            library_df = DataFrame(
                PeptideSequence = ["PEPTIDE1", "EDITPEP1"],
                PrecursorCharge = [2, 2],
                EntrapmentGroupId = [0, 1],
                PrecursorIdx = [100, 100]
            )
            library_path = joinpath(temp_dir, "library.tsv")
            CSV.write(library_path, library_df, delim='\t')
            
            # Empty results
            empty_results = DataFrame(
                stripped_seq = String[],
                z = UInt8[],
                PredVal = Float32[],
                decoy = Bool[],
                channel = UInt8[],
                file_name = String[],
                protein = String[]
            )
            
            EntrapmentAnalyses.load_parquet_results(files::Vector{String}) = empty_results
            
            # Should handle empty data gracefully
            result_df = EntrapmentAnalyses.run_efdr_analysis(
                ["empty.parquet"],
                library_path;
                output_dir = joinpath(temp_dir, "empty_output"),
                show_progress = false
            )
            
            @test nrow(result_df) == 0
        end
        
        @testset "All decoys" begin
            # Create results with only decoys
            library_df = DataFrame(
                PeptideSequence = ["PEPTIDE1", "EDITPEP1"],
                PrecursorCharge = [2, 2],
                EntrapmentGroupId = [0, 1],
                PrecursorIdx = [100, 100]
            )
            library_path = joinpath(temp_dir, "library_decoys.tsv")
            CSV.write(library_path, library_df, delim='\t')
            
            decoy_results = DataFrame(
                stripped_seq = ["PEPTIDE1", "EDITPEP1"],
                z = UInt8[2, 2],
                PredVal = Float32[0.9, 0.8],
                decoy = [true, true],  # All decoys
                channel = UInt8[0, 0],
                file_name = ["run1", "run1"],
                protein = ["PROT1", "PROT1"]
            )
            
            EntrapmentAnalyses.load_parquet_results(files::Vector{String}) = decoy_results
            
            result_df = EntrapmentAnalyses.run_efdr_analysis(
                ["decoys.parquet"],
                library_path;
                output_dir = joinpath(temp_dir, "decoy_output"),
                show_progress = false
            )
            
            # Should have no results after removing decoys
            @test nrow(result_df) == 0
        end
        
        # Clean up
        rm(temp_dir, recursive=true, force=true)
    end
    
    @testset "Multiple files" begin
        temp_dir, library_path, _ = create_test_files()
        
        # Create multi-file results
        multi_results = DataFrame(
            stripped_seq = ["PEPTIDE1", "EDITPEP1", "PEPTIDE2", "EDITPEP2",
                           "PEPTIDE1", "EDITPEP1", "PEPTIDE3", "EDITPEP3"],
            z = UInt8[2, 2, 3, 3, 2, 2, 2, 2],
            PredVal = Float32[0.9, 0.8, 0.7, 0.6, 0.85, 0.75, 0.82, 0.72],
            decoy = [false, false, false, false, false, false, false, false],
            channel = UInt8[0, 0, 0, 0, 0, 0, 0, 0],
            file_name = ["run1", "run1", "run1", "run1", "run2", "run2", "run2", "run2"],
            protein = ["PROT1", "PROT1", "PROT2", "PROT2", "PROT1", "PROT1", "PROT3", "PROT3"]
        )
        
        EntrapmentAnalyses.load_parquet_results(files::Vector{String}) = multi_results
        
        result_df = EntrapmentAnalyses.run_efdr_analysis(
            ["file1.parquet", "file2.parquet"],  # Multiple files
            library_path;
            output_dir = joinpath(temp_dir, "multi_output"),
            show_progress = false
        )
        
        # Check that both runs are present
        @test "run1" in result_df.file_name
        @test "run2" in result_df.file_name
        
        # Check EFDR calculated per run
        @test hasproperty(result_df, :combined_entrapment_fdr)
        @test hasproperty(result_df, :precursor_entrapment_fdr)
        
        # Clean up
        rm(temp_dir, recursive=true, force=true)
    end
end