# Visualization functions with notebook styling

"""
    plot_precursor_efdr_comparison(df::DataFrame; kwargs...)

Create EFDR comparison plot for precursors using exact notebook styling.
"""
function plot_precursor_efdr_comparison(
    df::DataFrame;
    output_path::String = "precursor_efdr_comparison.pdf",
    fdr_col::Symbol = :local_qvalue,
    efdr_col::Symbol = :precursor_entrapment_fdr,
    title::String = "Entrapment Analysis Precursors"
)
    # Extract typed vectors
    fdr_values = df[!, fdr_col]::Vector{Float32}
    efdr_values = df[!, efdr_col]::Vector{Float32}
    
    # Exact plot styling from notebook
    p = plot(
        size = (400*1.5, 300*1.5),  # 600x450 pixels
        xlim = (0, 0.05),
        ylim = (0, 0.05),
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
    plot!(p, [0, 0.05], [0, 0.05], 
          color = :black, 
          lw = 3, 
          label = nothing, 
          linestyle = :dash)
    
    # Plot each run with notebook's exact color
    if hasproperty(df, :file_name)
        file_names = df.file_name::Vector{String}
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
"""
function plot_protein_efdr_comparison(
    df::DataFrame;
    output_path::String = "protein_efdr_comparison.pdf",
    fdr_col::Symbol = :Protein_Qvalue,
    efdr_col::Symbol = :protein_group_entrapment_fdr,
    title::String = "Entrapment Analysis Protein Groups"
)
    # Extract typed vectors
    fdr_values = df[!, fdr_col]::Vector{Float32}
    efdr_values = df[!, efdr_col]::Vector{Float32}
    
    # Exact plot styling from notebook
    p = plot(
        size = (400*1.5, 300*1.5),  # 600x450 pixels
        xlim = (0, 0.05),
        ylim = (0, 0.05),
        xlabel = "FDR",
        ylabel = "Entrapment FDR",
        title = title
    )
    
    # Note: titlefontsize and other font sizes are set after the diagonal line in notebook
    plot!(p, [0, 0.05], [0, 0.05], 
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
        file_names = df.file_name::Vector{String}
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