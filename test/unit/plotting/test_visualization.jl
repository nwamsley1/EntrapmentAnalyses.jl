using Test
using DataFrames
using EntrapmentAnalyses
using Plots

@testset "Visualization Tests" begin
    
    # Create test data
    function create_test_data()
        return DataFrame(
            local_qvalue = Float32[0.01, 0.02, 0.03, 0.04, 0.05, 0.06],
            combined_entrapment_fdr = Float32[0.008, 0.018, 0.025, 0.032, 0.04, 0.045],
            paired_entrapment_fdr = Float32[0.009, 0.019, 0.028, 0.035, 0.042, 0.048],
            precursor_entrapment_fdr = Float32[0.008, 0.017, 0.024, 0.031, 0.038, 0.044],
            protein_group_entrapment_fdr = Float32[0.01, 0.02, 0.029, 0.036, 0.043, 0.049],
            Protein_Qvalue = Float32[0.012, 0.022, 0.032, 0.042, 0.052, 0.062],
            decoy = [false, false, true, false, true, false],
            entrapment_group = [false, true, false, true, false, true],
            file_name = ["run1", "run1", "run1", "run2", "run2", "run2"]
        )
    end
    
    @testset "plot_precursor_efdr_comparison" begin
        df = create_test_data()
        
        # Test basic plotting
        @testset "Basic plot" begin
            temp_file = tempname() * ".pdf"
            p = EntrapmentAnalyses.plot_precursor_efdr_comparison(df; output_path=temp_file)
            
            @test isa(p, Plots.Plot)
            @test isfile(temp_file)
            
            # Clean up
            rm(temp_file, force=true)
        end
        
        # Test with custom xlim
        @testset "Custom xlim" begin
            temp_file = tempname() * ".pdf"
            p = EntrapmentAnalyses.plot_precursor_efdr_comparison(
                df; 
                output_path=temp_file,
                xlim=(0, 0.1)
            )
            
            @test isa(p, Plots.Plot)
            @test isfile(temp_file)
            
            rm(temp_file, force=true)
        end
        
        # Test single file (no file_name column)
        @testset "Single file" begin
            df_single = select(df, Not(:file_name))
            temp_file = tempname() * ".pdf"
            
            p = EntrapmentAnalyses.plot_precursor_efdr_comparison(
                df_single; 
                output_path=temp_file
            )
            
            @test isa(p, Plots.Plot)
            @test isfile(temp_file)
            
            rm(temp_file, force=true)
        end
    end
    
    @testset "plot_protein_efdr_comparison" begin
        df = create_test_data()
        
        # Test basic plotting
        @testset "Basic plot" begin
            temp_file = tempname() * ".pdf"
            p = EntrapmentAnalyses.plot_protein_efdr_comparison(df; output_path=temp_file)
            
            @test isa(p, Plots.Plot)
            @test isfile(temp_file)
            
            rm(temp_file, force=true)
        end
    end
    
    @testset "create_efdr_comparison_plot" begin
        # Test vector-based plotting
        fdr_values = Float32[0.01, 0.02, 0.03, 0.04, 0.05]
        efdr_values = Float32[0.008, 0.018, 0.025, 0.032, 0.04]
        
        p = EntrapmentAnalyses.create_efdr_comparison_plot(
            fdr_values, 
            efdr_values;
            title="Test Plot"
        )
        
        @test isa(p, Plots.Plot)
    end
    
    @testset "plot_combined_efdr" begin
        df = create_test_data()
        temp_file = tempname() * ".pdf"
        
        p = EntrapmentAnalyses.plot_combined_efdr(df; output_path=temp_file)
        
        @test isa(p, Plots.Plot)
        @test isfile(temp_file)
        @test isfile(replace(temp_file, ".pdf" => ".png"))  # Should also save PNG
        
        # Clean up
        rm(temp_file, force=true)
        rm(replace(temp_file, ".pdf" => ".png"), force=true)
    end
    
    @testset "plot_paired_efdr" begin
        df = create_test_data()
        temp_file = tempname() * ".pdf"
        
        p = EntrapmentAnalyses.plot_paired_efdr(df; output_path=temp_file)
        
        @test isa(p, Plots.Plot)
        @test isfile(temp_file)
        @test isfile(replace(temp_file, ".pdf" => ".png"))
        
        # Clean up
        rm(temp_file, force=true)
        rm(replace(temp_file, ".pdf" => ".png"), force=true)
    end
    
    @testset "plot_efdr_comparison_both_methods" begin
        df = create_test_data()
        temp_file = tempname() * ".pdf"
        
        p = EntrapmentAnalyses.plot_efdr_comparison_both_methods(df; output_path=temp_file)
        
        @test isa(p, Plots.Plot)
        @test isfile(temp_file)
        @test isfile(replace(temp_file, ".pdf" => ".png"))
        
        # Clean up
        rm(temp_file, force=true)
        rm(replace(temp_file, ".pdf" => ".png"), force=true)
    end
    
    @testset "generate_analysis_report" begin
        df = create_test_data()
        temp_dir = mktempdir()
        
        # Test report generation
        report_path = EntrapmentAnalyses.generate_analysis_report(df, temp_dir)
        
        @test isfile(report_path)
        @test isfile(joinpath(temp_dir, "combined_efdr.pdf"))
        @test isfile(joinpath(temp_dir, "combined_efdr.png"))
        @test isfile(joinpath(temp_dir, "paired_efdr.pdf"))
        @test isfile(joinpath(temp_dir, "paired_efdr.png"))
        @test isfile(joinpath(temp_dir, "comparison_both_methods.pdf"))
        @test isfile(joinpath(temp_dir, "comparison_both_methods.png"))
        
        # Check report content
        report_content = read(report_path, String)
        @test occursin("# Entrapment FDR Analysis Report", report_content)
        @test occursin("## Data Summary", report_content)
        @test occursin("## Methods", report_content)
        @test occursin("## Results at Common FDR Thresholds", report_content)
        @test occursin("## Visualizations", report_content)
        
        # Test with only combined EFDR
        df_combined_only = select(df, Not(:paired_entrapment_fdr))
        report_path2 = EntrapmentAnalyses.generate_analysis_report(
            df_combined_only, 
            temp_dir;
            paired_efdr_col = :nonexistent_col
        )
        @test isfile(report_path2)
        
        # Clean up
        rm(temp_dir, recursive=true, force=true)
    end
    
    @testset "Plot styling verification" begin
        # Test that plots have correct styling
        df = create_test_data()
        p = EntrapmentAnalyses.plot_combined_efdr(
            df; 
            output_path=tempname() * ".pdf"
        )
        
        # Verify plot attributes (size, fonts, etc.)
        @test p.attr[:size] == (600, 450)
        @test p.attr[:titlefontsize] == 16
        @test p.attr[:xguidefontsize] == 16
        @test p.attr[:yguidefontsize] == 16
        @test p.attr[:tickfontsize] == 12
    end
end