using Test
using DataFrames
using EntrapmentAnalyses

@testset "Q-value Calculation Tests" begin
    
    @testset "validate_sort_order" begin
        # Test ascending order
        @test EntrapmentAnalyses.validate_sort_order([1, 2, 3, 4], true) == true
        @test EntrapmentAnalyses.validate_sort_order([4, 3, 2, 1], true) == false
        
        # Test descending order
        @test EntrapmentAnalyses.validate_sort_order([4, 3, 2, 1], false) == true
        @test EntrapmentAnalyses.validate_sort_order([1, 2, 3, 4], false) == false
        
        # Test with duplicates
        @test EntrapmentAnalyses.validate_sort_order([1, 2, 2, 3], true) == true
        @test EntrapmentAnalyses.validate_sort_order([3, 2, 2, 1], false) == true
        
        # Test empty and single element
        @test EntrapmentAnalyses.validate_sort_order([], true) == true
        @test EntrapmentAnalyses.validate_sort_order([5], false) == true
    end
    
    @testset "calculate_qvalues" begin
        # Basic test with mixed targets and decoys
        @testset "Basic calculation" begin
            scores = Float32[0.9, 0.8, 0.7, 0.6, 0.5]
            is_decoy = [false, true, false, true, false]  # 3 targets, 2 decoys
            
            qvals = EntrapmentAnalyses.calculate_qvalues(scores, is_decoy)
            
            @test length(qvals) == 5
            @test all(0 .<= qvals .<= 1)
            @test issorted(qvals[sortperm(scores, rev=true)])  # Q-values should be monotonic
        end
        
        # Test with all targets
        @testset "All targets" begin
            scores = Float32[0.9, 0.8, 0.7]
            is_decoy = [false, false, false]
            
            qvals = EntrapmentAnalyses.calculate_qvalues(scores, is_decoy)
            
            @test all(qvals .== 0.0)  # No decoys means q-value = 0
        end
        
        # Test with all decoys
        @testset "All decoys" begin
            scores = Float32[0.9, 0.8, 0.7]
            is_decoy = [true, true, true]
            
            qvals = EntrapmentAnalyses.calculate_qvalues(scores, is_decoy)
            
            @test all(qvals .== 0.0)  # No targets means q-value = 0 (special case)
        end
        
        # Test monotonicity
        @testset "Monotonicity" begin
            # Create a case where raw q-values would not be monotonic
            scores = Float32[0.9, 0.85, 0.8, 0.75, 0.7, 0.65]
            is_decoy = [false, true, false, false, true, false]
            
            qvals = EntrapmentAnalyses.calculate_qvalues(scores, is_decoy)
            
            # Check that q-values are monotonic when sorted by score
            sorted_indices = sortperm(scores, rev=true)
            sorted_qvals = qvals[sorted_indices]
            @test issorted(sorted_qvals)
        end
        
        # Test input validation
        @testset "Input validation" begin
            scores = Float32[0.9, 0.8]
            is_decoy = [true, true, false]  # Wrong length
            
            @test_throws ErrorException EntrapmentAnalyses.calculate_qvalues(scores, is_decoy)
        end
    end
    
    @testset "calculate_qvalues!" begin
        # Create test DataFrame
        @testset "Basic DataFrame q-value calculation" begin
            df = DataFrame(
                PredVal = Float32[0.9, 0.8, 0.7, 0.6, 0.5],
                decoy = [false, true, false, true, false],
                stripped_seq = ["PEP1", "PEP2", "PEP3", "PEP4", "PEP5"],
                z = [2, 2, 3, 3, 2]
            )
            
            # Calculate q-values
            EntrapmentAnalyses.calculate_qvalues!(df)
            
            @test hasproperty(df, :local_qvalue)
            @test all(0 .<= df.local_qvalue .<= 1)
            @test issorted(df.PredVal, rev=true)  # Should be sorted by score
        end
        
        # Test without sorting
        @testset "Without sorting" begin
            df = DataFrame(
                PredVal = Float32[0.5, 0.9, 0.7],  # Not sorted
                decoy = [false, false, true]
            )
            
            EntrapmentAnalyses.calculate_qvalues!(df; sort_df=false)
            
            @test df.PredVal == Float32[0.5, 0.9, 0.7]  # Order unchanged
            @test hasproperty(df, :local_qvalue)
        end
        
        # Test missing columns
        @testset "Missing columns" begin
            df = DataFrame(
                score = Float32[0.9, 0.8],  # Wrong column name
                decoy = [false, true]
            )
            
            @test_throws ErrorException EntrapmentAnalyses.calculate_qvalues!(df)
        end
    end
    
    @testset "calculate_global_qvalues!" begin
        # Create test DataFrame with multiple PSMs per precursor
        @testset "Global q-value calculation" begin
            df = DataFrame(
                PredVal = Float32[0.9, 0.8, 0.85, 0.7, 0.75, 0.6],
                decoy = [false, false, false, true, true, true],
                channel = [0, 0, 0, 0, 0, 0],  # Dummy channel
                stripped_seq = ["PEP1", "PEP1", "PEP2", "PEP1", "PEP1", "PEP2"],
                z = [2, 2, 3, 2, 2, 3]
            )
            
            EntrapmentAnalyses.calculate_global_qvalues!(df)
            
            @test hasproperty(df, :global_qvalue)
            @test all(0 .<= df.global_qvalue .<= 1)
            
            # Check that peptides with same (seq, z, decoy) have same global q-value
            gdf = groupby(df, [:decoy, :stripped_seq, :z])
            for group in gdf
                @test length(unique(group.global_qvalue)) == 1
            end
        end
        
        # Test with single PSM per precursor
        @testset "Single PSM per precursor" begin
            df = DataFrame(
                PredVal = Float32[0.9, 0.8, 0.7],
                decoy = [false, true, false],
                channel = [0, 0, 0],  # Dummy channel
                stripped_seq = ["PEP1", "PEP2", "PEP3"],
                z = [2, 2, 3]
            )
            
            EntrapmentAnalyses.calculate_global_qvalues!(df)
            
            @test hasproperty(df, :global_qvalue)
            # Should essentially match local q-values in this case
        end
    end
    
    @testset "calculate_qvalues_per_file!" begin
        # Test with multiple files
        @testset "Multiple files" begin
            df = DataFrame(
                PredVal = Float32[0.9, 0.8, 0.7, 0.85, 0.75, 0.65],
                decoy = [false, true, false, false, true, false],
                channel = [0, 0, 0, 0, 0, 0],  # Dummy channel
                stripped_seq = ["PEP1", "PEP2", "PEP3", "PEP4", "PEP5", "PEP6"],
                z = [2, 2, 3, 2, 3, 2],
                file_name = ["file1", "file1", "file1", "file2", "file2", "file2"]
            )
            
            EntrapmentAnalyses.calculate_qvalues_per_file!(df)
            
            @test hasproperty(df, :local_qvalue)
            @test hasproperty(df, :global_qvalue)
            @test all(0 .<= df.local_qvalue .<= 1)
            @test all(0 .<= df.global_qvalue .<= 1)
            
            # Check that q-values are calculated per file
            # File 1: 2 targets, 1 decoy
            # File 2: 2 targets, 1 decoy
            file1_df = df[df.file_name .== "file1", :]
            file2_df = df[df.file_name .== "file2", :]
            
            # Each file should have its own q-value calculation
            @test maximum(file1_df.local_qvalue) ≤ 1.0
            @test maximum(file2_df.local_qvalue) ≤ 1.0
        end
        
        # Test without file column (single file)
        @testset "Single file (no file column)" begin
            df = DataFrame(
                PredVal = Float32[0.9, 0.8, 0.7],
                decoy = [false, true, false],
                channel = [0, 0, 0],  # Dummy channel
                stripped_seq = ["PEP1", "PEP2", "PEP3"],
                z = [2, 2, 3]
            )
            
            EntrapmentAnalyses.calculate_qvalues_per_file!(df)
            
            @test hasproperty(df, :local_qvalue)
            @test hasproperty(df, :global_qvalue)
            @test all(0 .<= df.local_qvalue .<= 1)
            @test all(0 .<= df.global_qvalue .<= 1)
        end
    end
    
end