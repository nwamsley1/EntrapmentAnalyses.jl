# Visualization functions with notebook styling

"""
    plot_precursor_efdr_comparison(df::DataFrame; kwargs...)

Create EFDR comparison plot for precursors using exact notebook styling.
"""
function plot_precursor_efdr_comparison(
    df::DataFrame;
    output_path::String = "precursor_efdr_comparison.pdf",
    fdr_col::Symbol = :local_qvalue,
    efdr_col::Symbol = :paired_entrapment_fdr,
    title::String = "Entrapment Analysis Precursors",
    xlim::Tuple{Real, Real} = (0, 0.05)
)
    # Extract typed vectors
    fdr_values = df[!, fdr_col]::AbstractVector{Float32}
    efdr_values = df[!, efdr_col]::AbstractVector{Float32}
    
    # Determine y-axis range from data
    mask = fdr_values .<= xlim[2]
    if any(mask)
        max_efdr = maximum(efdr_values[mask])
        ylim_max = max(max_efdr * 1.1, xlim[2])  # Add 10% padding or at least match x-axis
    else
        ylim_max = xlim[2]
    end
    
    # Plot styling from notebook with dynamic limits
    p = plot(
        size = (400*1.5, 300*1.5),  # 600x450 pixels
        xlim = xlim,
        ylim = (0, ylim_max),
        xlabel = "FDR",
        ylabel = "Entrapment FDR",
        title = title,
        titlefontsize = 16,
        xguidefontsize = 16,
        yguidefontsize = 16,
        tickfontsize = 12,
        legendfontsize = 12
    )
    
    # Diagonal reference line (exact style)
    plot!(p, [0, xlim[2]], [0, xlim[2]], 
          color = :black, 
          lw = 3, 
          label = nothing, 
          linestyle = :dash)
    
    # Plot each run with notebook's exact color
    if hasproperty(df, :file_name)
        file_names = df.file_name::AbstractVector{String}
        for (key, group) in pairs(groupby(DataFrame(fdr=fdr_values, efdr=efdr_values, file=file_names), :file))
            plot!(p,
                  group.fdr,
                  group.efdr,
                  lw = 3,
                  label = nothing,  # No labels as per notebook
                  color = RGB(0.39215686, 0.58431373, 0.92941176),  # Exact color
                  alpha = 0.75)
        end
    else
        # Single run case
        plot!(p,
              fdr_values,
              efdr_values,
              lw = 3,
              label = nothing,
              color = RGB(0.39215686, 0.58431373, 0.92941176),
              alpha = 0.75)
    end
    
    savefig(p, output_path)
    println("Plot saved to: $output_path")
    
    return p
end

"""
    plot_protein_efdr_comparison(df::DataFrame; kwargs...)

Create EFDR comparison plot for protein groups using exact notebook styling.
Supports both combined and paired EFDR columns.
"""
function plot_protein_efdr_comparison(
    df::DataFrame;
    output_path::String = "protein_efdr_comparison.pdf",
    fdr_col::Symbol = :Protein_Qvalue,
    efdr_col::Symbol = :protein_group_entrapment_fdr,
    combined_efdr_col::Union{Symbol, Nothing} = :combined_protein_fdr,
    title::String = "Entrapment Analysis Protein Groups",
    xlim::Tuple{Real, Real} = (0, 0.05)
)
    # Extract typed vectors
    fdr_values = df[!, fdr_col]::AbstractVector{Float32}
    efdr_values = df[!, efdr_col]::AbstractVector{Float32}
    
    # Determine y-axis range from data
    mask = fdr_values .<= xlim[2]
    if any(mask)
        max_efdr = maximum(efdr_values[mask])
        ylim_max = max(max_efdr * 1.1, xlim[2])  # Add 10% padding or at least match x-axis
    else
        ylim_max = xlim[2]
    end
    
    # Plot styling with dynamic limits
    p = plot(
        size = (400*1.5, 300*1.5),  # 600x450 pixels
        xlim = xlim,
        ylim = (0, ylim_max),
        xlabel = "FDR",
        ylabel = "Entrapment FDR",
        title = title
    )
    
    # Note: titlefontsize and other font sizes are set after the diagonal line in notebook
    plot!(p, [0, xlim[2]], [0, xlim[2]], 
          color = :black, 
          lw = 3, 
          label = nothing, 
          linestyle = :dash,
          titlefontsize = 16,
          xguidefontsize = 16,
          yguidefontsize = 16,
          tickfontsize = 12,
          legendfontsize = 12)
    
    # Plot each run
    if hasproperty(df, :file_name)
        file_names = df.file_name::AbstractVector{String}
        for (key, group) in pairs(groupby(DataFrame(fdr=fdr_values, efdr=efdr_values, file=file_names), :file))
            plot!(p,
                  group.fdr,
                  group.efdr,
                  lw = 3,
                  label = nothing,
                  color = RGB(0.39215686, 0.58431373, 0.92941176),
                  alpha = 0.75)
        end
    else
        # Single run case
        plot!(p,
              fdr_values,
              efdr_values,
              lw = 3,
              label = nothing,
              color = RGB(0.39215686, 0.58431373, 0.92941176),
              alpha = 0.75)
    end
    
    savefig(p, output_path)
    println("Plot saved to: $output_path")
    
    return p
end

"""
    create_efdr_comparison_plot(fdr_values::AbstractVector{T}, efdr_values::AbstractVector{T};
                               title::String, xlim=(0, 0.05), ylim=(0, 0.05)) where T<:Real

Low-level vector-based plotting function.
"""
function create_efdr_comparison_plot(
    fdr_values::AbstractVector{T},
    efdr_values::AbstractVector{T};
    title::String = "EFDR Comparison",
    xlim = (0, 0.05),
    ylim = (0, 0.05)
) where T<:Real
    
    p = plot(
        size = (400*1.5, 300*1.5),
        xlim = xlim,
        ylim = ylim,
        xlabel = "FDR",
        ylabel = "Entrapment FDR",
        title = title,
        titlefontsize = 16,
        xguidefontsize = 16,
        yguidefontsize = 16,
        tickfontsize = 12,
        legendfontsize = 12
    )
    
    # Diagonal reference
    plot!(p, [xlim[1], xlim[2]], [ylim[1], ylim[2]], 
          color = :black, 
          lw = 3, 
          label = nothing, 
          linestyle = :dash)
    
    # Data
    plot!(p,
          fdr_values,
          efdr_values,
          lw = 3,
          label = nothing,
          color = RGB(0.39215686, 0.58431373, 0.92941176),
          alpha = 0.75)
    
    return p
end

"""
    plot_combined_efdr(df::DataFrame; kwargs...)

Create a plot showing only the combined EFDR results.

# Arguments
- `df`: DataFrame with EFDR results
- `output_path`: Output file path (default: "combined_efdr.pdf")
- `fdr_col`: FDR column name (default: :local_qvalue)
- `efdr_col`: EFDR column name (default: :combined_entrapment_fdr)
"""
function plot_combined_efdr(
    df::DataFrame;
    output_path::String = "combined_efdr.pdf",
    fdr_col::Symbol = :local_qvalue,
    efdr_col::Symbol = :combined_entrapment_fdr,
    title::String = "Combined EFDR Method",
    xlim::Tuple{Real, Real} = (0, 0.05)
)
    # Extract typed vectors
    fdr_values = df[!, fdr_col]::AbstractVector{Float32}
    efdr_values = df[!, efdr_col]::AbstractVector{Float32}
    
    # Determine y-axis range from data
    mask = fdr_values .<= xlim[2]
    if any(mask)
        max_efdr = maximum(efdr_values[mask])
        ylim_max = max(max_efdr * 1.1, xlim[2])
    else
        ylim_max = xlim[2]
    end
    
    # Create plot with notebook styling
    p = plot(
        size = (400*1.5, 300*1.5),  # 600x450 pixels
        xlim = xlim,
        ylim = (0, ylim_max),
        xlabel = "FDR",
        ylabel = "Entrapment FDR",
        title = title,
        titlefontsize = 16,
        xguidefontsize = 16,
        yguidefontsize = 16,
        tickfontsize = 12,
        legendfontsize = 12
    )
    
    # Diagonal reference line
    plot!(p, [0, xlim[2]], [0, xlim[2]], 
          color = :black, 
          lw = 3, 
          label = nothing, 
          linestyle = :dash)
    
    # Plot data
    if hasproperty(df, :file_name)
        file_names = df.file_name::AbstractVector{String}
        for (key, group) in pairs(groupby(DataFrame(fdr=fdr_values, efdr=efdr_values, file=file_names), :file))
            plot!(p,
                  group.fdr,
                  group.efdr,
                  lw = 3,
                  label = nothing,
                  color = RGB(0.39215686, 0.58431373, 0.92941176),
                  alpha = 0.75)
        end
    else
        plot!(p,
              fdr_values,
              efdr_values,
              lw = 3,
              label = nothing,
              color = RGB(0.39215686, 0.58431373, 0.92941176),
              alpha = 0.75)
    end
    
    # Save as both PDF and PNG
    savefig(p, output_path)
    png_path = replace(output_path, ".pdf" => ".png")
    savefig(p, png_path)
    println("Plot saved to: $output_path and $png_path")
    
    return p
end

"""
    plot_paired_efdr(df::DataFrame; kwargs...)

Create a plot showing only the paired EFDR results.

# Arguments
- `df`: DataFrame with EFDR results
- `output_path`: Output file path (default: "paired_efdr.pdf")
- `fdr_col`: FDR column name (default: :local_qvalue)
- `efdr_col`: EFDR column name (default: :paired_entrapment_fdr)
"""
function plot_paired_efdr(
    df::DataFrame;
    output_path::String = "paired_efdr.pdf",
    fdr_col::Symbol = :local_qvalue,
    efdr_col::Symbol = :paired_entrapment_fdr,
    title::String = "Paired EFDR Method",
    xlim::Tuple{Real, Real} = (0, 0.05)
)
    # Extract typed vectors
    fdr_values = df[!, fdr_col]::AbstractVector{Float32}
    efdr_values = df[!, efdr_col]::AbstractVector{Float32}
    
    # Determine y-axis range from data
    mask = fdr_values .<= xlim[2]
    if any(mask)
        max_efdr = maximum(efdr_values[mask])
        ylim_max = max(max_efdr * 1.1, xlim[2])
    else
        ylim_max = xlim[2]
    end
    
    # Create plot with notebook styling
    p = plot(
        size = (400*1.5, 300*1.5),  # 600x450 pixels
        xlim = xlim,
        ylim = (0, ylim_max),
        xlabel = "FDR",
        ylabel = "Entrapment FDR",
        title = title,
        titlefontsize = 16,
        xguidefontsize = 16,
        yguidefontsize = 16,
        tickfontsize = 12,
        legendfontsize = 12
    )
    
    # Diagonal reference line
    plot!(p, [0, xlim[2]], [0, xlim[2]], 
          color = :black, 
          lw = 3, 
          label = nothing, 
          linestyle = :dash)
    
    # Plot data
    if hasproperty(df, :file_name)
        file_names = df.file_name::AbstractVector{String}
        for (key, group) in pairs(groupby(DataFrame(fdr=fdr_values, efdr=efdr_values, file=file_names), :file))
            plot!(p,
                  group.fdr,
                  group.efdr,
                  lw = 3,
                  label = nothing,
                  color = RGB(0.39215686, 0.58431373, 0.92941176),
                  alpha = 0.75)
        end
    else
        plot!(p,
              fdr_values,
              efdr_values,
              lw = 3,
              label = nothing,
              color = RGB(0.39215686, 0.58431373, 0.92941176),
              alpha = 0.75)
    end
    
    # Save as both PDF and PNG
    savefig(p, output_path)
    png_path = replace(output_path, ".pdf" => ".png")
    savefig(p, png_path)
    println("Plot saved to: $output_path and $png_path")
    
    return p
end

"""
    plot_efdr_comparison_both_methods(df::DataFrame; kwargs...)

Create a comparison plot showing both combined and paired EFDR methods.

# Arguments
- `df`: DataFrame with both combined and paired EFDR results
- `output_path`: Output file path (default: "efdr_comparison_both_methods.pdf")
- `fdr_col`: FDR column name (default: :local_qvalue)
- `combined_efdr_col`: Combined EFDR column name (default: :combined_entrapment_fdr)
- `paired_efdr_col`: Paired EFDR column name (default: :paired_entrapment_fdr)
"""
function plot_efdr_comparison_both_methods(
    df::DataFrame;
    output_path::String = "efdr_comparison_both_methods.pdf",
    fdr_col::Symbol = :local_qvalue,
    combined_efdr_col::Symbol = :combined_entrapment_fdr,
    paired_efdr_col::Symbol = :paired_entrapment_fdr,
    title::String = "EFDR Methods Comparison",
    xlim::Tuple{Real, Real} = (0, 0.05)
)
    # Extract typed vectors
    fdr_values = df[!, fdr_col]::AbstractVector{Float32}
    combined_efdr = df[!, combined_efdr_col]::AbstractVector{Float32}
    paired_efdr = df[!, paired_efdr_col]::AbstractVector{Float32}
    
    # Determine y-axis range from data
    mask = fdr_values .<= xlim[2]
    if any(mask)
        max_efdr = max(maximum(combined_efdr[mask]), maximum(paired_efdr[mask]))
        ylim_max = max(max_efdr * 1.1, xlim[2])
    else
        ylim_max = xlim[2]
    end
    
    # Create plot with notebook styling
    p = plot(
        size = (400*1.5, 300*1.5),  # 600x450 pixels
        xlim = xlim,
        ylim = (0, ylim_max),
        xlabel = "FDR",
        ylabel = "Entrapment FDR",
        title = title,
        titlefontsize = 16,
        xguidefontsize = 16,
        yguidefontsize = 16,
        tickfontsize = 12,
        legendfontsize = 12,
        legend = :bottomright
    )
    
    # Diagonal reference line
    plot!(p, [0, xlim[2]], [0, xlim[2]], 
          color = :black, 
          lw = 3, 
          label = "y = x", 
          linestyle = :dash)
    
    # Plot both methods
    if hasproperty(df, :file_name)
        file_names = df.file_name::AbstractVector{String}
        
        # Combined method
        for (i, (key, group)) in enumerate(pairs(groupby(DataFrame(fdr=fdr_values, efdr=combined_efdr, file=file_names), :file)))
            plot!(p,
                  group.fdr,
                  group.efdr,
                  lw = 3,
                  label = i == 1 ? "Combined" : nothing,
                  color = RGB(0.39215686, 0.58431373, 0.92941176),
                  alpha = 0.75)
        end
        
        # Paired method
        for (i, (key, group)) in enumerate(pairs(groupby(DataFrame(fdr=fdr_values, efdr=paired_efdr, file=file_names), :file)))
            plot!(p,
                  group.fdr,
                  group.efdr,
                  lw = 3,
                  label = i == 1 ? "Paired" : nothing,
                  color = RGB(0.92941176, 0.58431373, 0.39215686),  # Complementary color
                  alpha = 0.75)
        end
    else
        # Single file case
        plot!(p,
              fdr_values,
              combined_efdr,
              lw = 3,
              label = "Combined",
              color = RGB(0.39215686, 0.58431373, 0.92941176),
              alpha = 0.75)
        
        plot!(p,
              fdr_values,
              paired_efdr,
              lw = 3,
              label = "Paired",
              color = RGB(0.92941176, 0.58431373, 0.39215686),
              alpha = 0.75)
    end
    
    # Save as both PDF and PNG
    savefig(p, output_path)
    png_path = replace(output_path, ".pdf" => ".png")
    savefig(p, png_path)
    println("Plot saved to: $output_path and $png_path")
    
    return p
end

"""
    generate_analysis_report(df::DataFrame, output_dir::String; kwargs...)

Generate a markdown report with embedded plots and analysis summary.

# Arguments
- `df`: DataFrame with EFDR results
- `output_dir`: Directory to save report and plots
- `paired_df`: Optional DataFrame with paired EFDR results (if calculated separately)
"""
function generate_analysis_report(
    df::DataFrame,
    output_dir::String;
    paired_df::Union{DataFrame, Nothing} = nothing,
    combined_efdr_col::Symbol = :combined_entrapment_fdr,
    paired_efdr_col::Symbol = :paired_entrapment_fdr,
    fdr_col::Symbol = :local_qvalue,
    xlim::Tuple{Real, Real} = (0, 0.05)
)
    # Ensure output directory exists
    mkpath(output_dir)
    
    # Generate plots
    println("Generating plots for report...")
    
    # Determine if we have both methods in one dataframe
    has_combined = hasproperty(df, combined_efdr_col)
    has_paired = hasproperty(df, paired_efdr_col)
    
    if has_combined
        plot_combined_efdr(df; 
            output_path = joinpath(output_dir, "combined_efdr.pdf"),
            efdr_col = combined_efdr_col,
            xlim = xlim)
    end
    
    if has_paired
        plot_paired_efdr(df; 
            output_path = joinpath(output_dir, "paired_efdr.pdf"),
            efdr_col = paired_efdr_col,
            xlim = xlim)
    end
    
    if has_combined && has_paired
        plot_efdr_comparison_both_methods(df; 
            output_path = joinpath(output_dir, "comparison_both_methods.pdf"),
            combined_efdr_col = combined_efdr_col,
            paired_efdr_col = paired_efdr_col,
            xlim = xlim)
    end
    
    # Generate report
    report_path = joinpath(output_dir, "analysis_report.md")
    println("Generating markdown report...")
    
    open(report_path, "w") do io
        # Header
        write(io, "# Entrapment FDR Analysis Report\n\n")
        write(io, "Generated on: $(Dates.now())\n\n")
        
        # Data summary
        write(io, "## Data Summary\n\n")
        n_psms = nrow(df)
        n_targets = sum(.!df.decoy)
        n_entrapments = sum(df.entrapment_group)
        write(io, "- Total PSMs: $(n_psms)\n")
        write(io, "- Target PSMs: $(n_targets)\n")
        write(io, "- Entrapment PSMs: $(n_entrapments)\n")
        
        if hasproperty(df, :file_name)
            n_files = length(unique(df.file_name))
            write(io, "- Number of files: $(n_files)\n")
        end
        write(io, "\n")
        
        # Method descriptions
        write(io, "## Methods\n\n")
        
        if has_combined
            write(io, "### Combined EFDR\n")
            write(io, "The combined empirical FDR method calculates:\n")
            write(io, "```\n")
            write(io, "EFDR = (Nε × (1 + 1/r)) / (Nτ + Nε)\n")
            write(io, "```\n")
            write(io, "where Nε is the number of entrapments and Nτ is the number of targets above the score threshold.\n\n")
        end
        
        if has_paired
            write(io, "### Paired EFDR\n")
            write(io, "The paired empirical FDR method calculates:\n")
            write(io, "```\n")
            write(io, "EFDR = (Nε + Nεsτ + 2×Nετs) / (Nτ + Nε)\n")
            write(io, "```\n")
            write(io, "where:\n")
            write(io, "- Nεsτ: entrapments winning with score ≥ threshold > paired target\n")
            write(io, "- Nετs: entrapments winning with both ≥ threshold\n\n")
        end
        
        # Results at common FDR thresholds
        write(io, "## Results at Common FDR Thresholds\n\n")
        
        thresholds = [0.01, 0.02, 0.05]
        write(io, "| FDR Threshold | ")
        if has_combined
            write(io, "Combined EFDR | ")
        end
        if has_paired
            write(io, "Paired EFDR | ")
        end
        write(io, "PSMs |\n")
        
        write(io, "|")
        write(io, "---------------|")
        if has_combined
            write(io, "--------------|")
        end
        if has_paired
            write(io, "------------|")
        end
        write(io, "-----|\n")
        
        for thresh in thresholds
            mask = df[!, fdr_col] .<= thresh
            n_at_thresh = sum(mask)
            write(io, "| $(thresh) | ")
            
            if has_combined && n_at_thresh > 0
                efdr_at_thresh = maximum(df[mask, combined_efdr_col])
                write(io, "$(round(efdr_at_thresh, digits=4)) | ")
            elseif has_combined
                write(io, "N/A | ")
            end
            
            if has_paired && n_at_thresh > 0
                efdr_at_thresh = maximum(df[mask, paired_efdr_col])
                write(io, "$(round(efdr_at_thresh, digits=4)) | ")
            elseif has_paired
                write(io, "N/A | ")
            end
            
            write(io, "$(n_at_thresh) |\n")
        end
        write(io, "\n")
        
        # Plots
        write(io, "## Visualizations\n\n")
        
        if has_combined
            write(io, "### Combined EFDR\n")
            write(io, "![Combined EFDR](combined_efdr.png)\n\n")
        end
        
        if has_paired
            write(io, "### Paired EFDR\n")
            write(io, "![Paired EFDR](paired_efdr.png)\n\n")
        end
        
        if has_combined && has_paired
            write(io, "### Method Comparison\n")
            write(io, "![EFDR Methods Comparison](comparison_both_methods.png)\n\n")
        end
        
        # Footer
        write(io, "## Notes\n\n")
        write(io, "- Q-values were calculated separately for each file\n")
        write(io, "- Plots use the exact styling from the reference notebook\n")
        write(io, "- FDR values have been monotonized to ensure proper ordering\n")
    end
    
    println("Report saved to: $report_path")
    return report_path
end