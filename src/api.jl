# Main API functions


"""
    run_efdr_analysis(parquet_files::Vector{String}, library_path::String; kwargs...)

Run empirical FDR analysis at the precursor level.

# Arguments
- `parquet_files`: Vector of paths to Parquet files with PSM results
- `library_path`: Path to TSV file with spectral library

# Keyword Arguments
- `output_dir`: Output directory (default: "efdr_output")
- `score_col`: Score column name (default: :PredVal)
- `global_qval_threshold`: Global q-value filter threshold (default: 1.0 = no filter)
- `local_qval_threshold`: Local q-value filter threshold (default: 0.05)
- `r_lib`: Library to real entrapment ratio (default: 1.0)
- `show_progress`: Show progress bars (default: true)

# Returns
- DataFrame with EFDR columns added
"""
function run_efdr_analysis(
    parquet_files::Vector{String},
    library_path::String;
    output_dir::String = "efdr_output",
    score_col::Symbol = :PredVal,
    global_qval_threshold::Float64 = 1.0,
    local_qval_threshold::Float64 = 0.01,
    r_lib::Float64 = 1.0,
    show_progress::Bool = true
)
    # Create output directory
    mkpath(output_dir)
    
    # Step 1: Load data
    results_df = load_parquet_results(parquet_files)
    library_df = load_spectral_library(library_path)
    
    # Step 2: Add entrapment labels
    target_seqs = Set{String}()
    entrap_seqs = Set{String}()
    
    entrap_groups = library_df.EntrapmentGroupId::AbstractVector{Int}
    peptide_seqs = library_df.PeptideSequence::AbstractVector{String}
    
    for i in eachindex(entrap_groups)
        if entrap_groups[i] == 0
            push!(target_seqs, peptide_seqs[i])
        else
            push!(entrap_seqs, peptide_seqs[i])
        end
    end
    
    # Type-stable column addition
    stripped_seqs = results_df.stripped_seq::AbstractVector{String}
    entrapment_group = [seq ∈ entrap_seqs for seq in stripped_seqs]
    results_df[!, :entrapment_group] = entrapment_group
    
    # Step 3: Calculate q-values per file
    calculate_qvalues_per_file!(results_df; score_col=score_col, file_col=:file_name)
    
    # Step 4: Filter by global q-value (optional)
    if global_qval_threshold < 1.0
        n_before = nrow(results_df)
        global_qvals = results_df.global_qvalue::AbstractVector{Float32}
        keep_mask = global_qvals .<= global_qval_threshold
        results_df = results_df[keep_mask, :]
        n_after = nrow(results_df)
        println("Filtered to $(n_after)/$(n_before) PSMs at $(global_qval_threshold) global FDR")
    end
    
    # Step 5: Filter by local q-value to reduce computation time
    if local_qval_threshold < 1.0
        n_before = nrow(results_df)
        local_qvals = results_df.local_qvalue::AbstractVector{Float32}
        keep_mask = local_qvals .<= local_qval_threshold
        results_df = results_df[keep_mask, :]
        n_after = nrow(results_df)
        println("Filtered to $(n_after)/$(n_before) PSMs at $(local_qval_threshold) local FDR")
    end
    
    # Step 6: Remove decoys for EFDR calculation
    decoy_mask = results_df.decoy::AbstractVector{Bool}
    results_no_decoys = results_df[.!decoy_mask, :]
    println("Analyzing $(nrow(results_no_decoys)) target PSMs")
    
    # Step 7: Compute pairing vectors
    pairing_info = compute_pairing_vectors(library_df, results_no_decoys; 
                                          show_progress=show_progress)
    
    # Step 8: Initialize EFDR columns
    results_no_decoys[!, :combined_entrapment_fdr] = zeros(Float32, nrow(results_no_decoys))
    results_no_decoys[!, :precursor_entrapment_fdr] = zeros(Float32, nrow(results_no_decoys))
    
    # Step 9: Get entrapment labels from pairing info
    entrap_labels = [orig ? 0 : 1 for orig in pairing_info.is_original]
    
    # Step 10: Group by run and calculate both EFDR methods
    if hasproperty(results_no_decoys, :file_name)
        # Add row indices for tracking
        results_no_decoys[!, :_row_idx] = 1:nrow(results_no_decoys)
        gdf = groupby(results_no_decoys, :file_name)
        
        for (key, run_df) in pairs(gdf)
            println("\nProcessing run: $(key.file_name)")
            
            # Get indices for this run in the original results_no_decoys
            run_indices = run_df._row_idx
            
            # Extract data for this run
            scores = run_df[!, score_col]::AbstractVector{Float32}
            qvals = run_df.local_qvalue::AbstractVector{Float32}
            run_entrap_labels = entrap_labels[run_indices]
            
            # Calculate combined EFDR for this run
            combined_efdr_values = calculate_combined_efdr(
                Float64.(scores),
                run_entrap_labels,
                Float64.(qvals);
                r = r_lib,
                show_progress = show_progress
            )
            monotonize!(combined_efdr_values)
            
            # Debug: Check combined EFDR
            max_combined = maximum(combined_efdr_values)
            println("  Max combined EFDR for this run: $max_combined")
            
            # Store combined EFDR results
            results_no_decoys[run_indices, :combined_entrapment_fdr] = Float32.(combined_efdr_values)
            
            # Get complement scores for paired EFDR
            all_scores = results_no_decoys[!, score_col]::AbstractVector{Float32}
            complement_scores_all = get_complement_scores(all_scores, pairing_info.complement_indices)
            complement_scores = complement_scores_all[run_indices]
            
            # Calculate paired EFDR
            paired_efdr_values = calculate_paired_efdr(
                Float64.(scores),
                Float64.(complement_scores),
                pairing_info.is_original[run_indices],
                Float64.(qvals);
                r = r_lib,
                show_progress = show_progress
            )
            
            monotonize!(paired_efdr_values)
            
            # Debug: Check paired EFDR
            max_paired = maximum(paired_efdr_values)
            println("  Max paired EFDR for this run: $max_paired")
            
            # Store paired EFDR results
            results_no_decoys[run_indices, :precursor_entrapment_fdr] = Float32.(paired_efdr_values)
        end
        
        # Remove temporary column
        select!(results_no_decoys, Not(:_row_idx))
    else
        # Single file case
        println("Processing single file...")
        
        scores = results_no_decoys[!, score_col]::AbstractVector{Float32}
        qvals = results_no_decoys.local_qvalue::AbstractVector{Float32}
        
        # Calculate combined EFDR
        combined_efdr_values = calculate_combined_efdr(
            Float64.(scores),
            entrap_labels,
            Float64.(qvals);
            r = r_lib,
            show_progress = show_progress
        )
        monotonize!(combined_efdr_values)
        results_no_decoys[!, :combined_entrapment_fdr] = Float32.(combined_efdr_values)
        
        # Calculate paired EFDR
        complement_scores = get_complement_scores(scores, pairing_info.complement_indices)
        paired_efdr_values = calculate_paired_efdr(
            Float64.(scores),
            Float64.(complement_scores),
            pairing_info.is_original,
            Float64.(qvals);
            r = r_lib,
            show_progress = show_progress
        )
        
        monotonize!(paired_efdr_values)
        results_no_decoys[!, :precursor_entrapment_fdr] = Float32.(paired_efdr_values)
    end
    
    # Step 11: Generate outputs
    output_file = joinpath(output_dir, "precursor_entrapment_results.tsv")
    CSV.write(output_file, results_no_decoys, delim='\t')
    println("\nResults saved to: $output_file")
    
    # Step 12: Create visualizations and report
    # Set x-axis limit to minimum of the two thresholds
    xlim_max = min(local_qval_threshold, global_qval_threshold)
    
    # Generate paired EFDR plot (for backward compatibility)
    plot_precursor_efdr_comparison(
        results_no_decoys;
        output_path = joinpath(output_dir, "precursor_efdr_comparison.pdf"),
        title = "Entrapment Analysis Precursors",
        xlim = (0, xlim_max)
    )
    
    # Generate comprehensive report with all plots
    generate_analysis_report(
        results_no_decoys,
        output_dir;
        combined_efdr_col = :combined_entrapment_fdr,
        paired_efdr_col = :precursor_entrapment_fdr,
        fdr_col = :local_qvalue,
        xlim = (0, xlim_max)
    )
    
    return results_no_decoys
end

"""
    run_protein_efdr_analysis(parquet_files::Vector{String}, library_path::String; kwargs...)

Run empirical FDR analysis at the protein level following the notebook implementation.

# Arguments
- `parquet_files`: Vector of paths to Parquet files with PSM results
- `library_path`: Path to TSV file with spectral library

# Keyword Arguments
- `output_dir`: Output directory (default: "efdr_output")
- `score_col`: Score column name (default: :PredVal)
- `global_qval_threshold`: Global q-value filter threshold (default: 1.0 = no filter)
- `local_qval_threshold`: Local q-value filter threshold (default: 0.05)
- `r_lib`: Library to real entrapment ratio (default: 1.0)
- `show_progress`: Show progress bars (default: true)

# Returns
- DataFrame with protein-level EFDR analysis
"""
function run_protein_efdr_analysis(
    parquet_files::Vector{String},
    library_path::String;
    output_dir::String = "efdr_output",
    score_col::Symbol = :PredVal,
    global_qval_threshold::Float64 = 1.0,
    local_qval_threshold::Float64 = 0.05,
    r_lib::Float64 = 1.0,
    show_progress::Bool = true
)
    # Create output directory
    mkpath(output_dir)
    
    # Step 1: Load data
    results_df = load_parquet_results(parquet_files)
    library_df = load_spectral_library(library_path)
    
    # Step 2: Add entrapment labels
    target_seqs = Set{String}()
    entrap_seqs = Set{String}()
    
    entrap_groups = library_df.EntrapmentGroupId::AbstractVector{Int}
    peptide_seqs = library_df.PeptideSequence::AbstractVector{String}
    
    for i in eachindex(entrap_groups)
        if entrap_groups[i] == 0
            push!(target_seqs, peptide_seqs[i])
        else
            push!(entrap_seqs, peptide_seqs[i])
        end
    end
    
    # Type-stable column addition
    stripped_seqs = results_df.stripped_seq::AbstractVector{String}
    entrapment_group = [seq ∈ entrap_seqs for seq in stripped_seqs]
    results_df[!, :entrapment_group] = entrapment_group
    
    # Step 3: Calculate q-values per file at PSM level
    calculate_qvalues_per_file!(results_df; score_col=score_col, file_col=:file_name)
    
    # Step 4: Filter by global q-value (optional)
    if global_qval_threshold < 1.0
        n_before = nrow(results_df)
        global_qvals = results_df.global_qvalue::AbstractVector{Float32}
        keep_mask = global_qvals .<= global_qval_threshold
        results_df = results_df[keep_mask, :]
        n_after = nrow(results_df)
        println("Filtered to $(n_after)/$(n_before) PSMs at $(global_qval_threshold) global FDR")
    end
    
    # Step 5: Filter by local q-value to reduce computation time
    if local_qval_threshold < 1.0
        n_before = nrow(results_df)
        local_qvals = results_df.local_qvalue::AbstractVector{Float32}
        keep_mask = local_qvals .<= local_qval_threshold
        results_df = results_df[keep_mask, :]
        n_after = nrow(results_df)
        println("Filtered to $(n_after)/$(n_before) PSMs at $(local_qval_threshold) local FDR")
    end
    
    # Step 6: Prepare protein-level data (rollup to best PSM per protein)
    protein_df = prepare_protein_analysis(results_df, library_df; score_col=score_col)
    
    # Step 7: Remove decoys for EFDR calculation
    decoy_mask = protein_df.decoy::AbstractVector{Bool}
    protein_no_decoys = protein_df[.!decoy_mask, :]
    println("Analyzing $(nrow(protein_no_decoys)) target protein groups")
    
    # Step 8: Compute pairing vectors for proteins
    pairing_info = compute_pairing_vectors(library_df, protein_no_decoys; 
                                          show_progress=show_progress)
    
    # Step 9: Initialize EFDR columns
    protein_no_decoys[!, :combined_protein_fdr] = zeros(Float32, nrow(protein_no_decoys))
    protein_no_decoys[!, :protein_group_entrapment_fdr] = zeros(Float32, nrow(protein_no_decoys))
    
    # Step 10: Get entrapment labels from pairing info
    entrap_labels = [orig ? 0 : 1 for orig in pairing_info.is_original]
    
    # Step 11: Group by run and calculate both EFDR methods
    if hasproperty(protein_no_decoys, :file_name)
        # Add row indices for tracking
        protein_no_decoys[!, :_row_idx] = 1:nrow(protein_no_decoys)
        gdf = groupby(protein_no_decoys, :file_name)
        
        for (key, run_df) in pairs(gdf)
            println("\nProcessing protein run: $(key.file_name)")
            
            # Get indices for this run in the original protein_no_decoys
            run_indices = run_df._row_idx
            
            # Extract data for this run
            scores = run_df[!, score_col]::AbstractVector{Float32}
            qvals = run_df.Protein_Qvalue::AbstractVector{Float32}
            run_entrap_labels = entrap_labels[run_indices]
            
            # Calculate combined EFDR for this run
            combined_efdr_values = calculate_combined_efdr(
                Float64.(scores),
                run_entrap_labels,
                Float64.(qvals);
                r = r_lib,
                show_progress = show_progress
            )
            monotonize!(combined_efdr_values)
            
            # Debug: Check combined EFDR
            max_combined = maximum(combined_efdr_values)
            println("  Max combined protein EFDR for this run: $max_combined")
            
            # Store combined EFDR results
            protein_no_decoys[run_indices, :combined_protein_fdr] = Float32.(combined_efdr_values)
            
            # Get complement scores for paired EFDR
            all_scores = protein_no_decoys[!, score_col]::AbstractVector{Float32}
            complement_scores_all = get_complement_scores(all_scores, pairing_info.complement_indices)
            complement_scores = complement_scores_all[run_indices]
            
            # Calculate paired EFDR
            paired_efdr_values = calculate_paired_efdr(
                Float64.(scores),
                Float64.(complement_scores),
                pairing_info.is_original[run_indices],
                Float64.(qvals);
                r = r_lib,
                show_progress = show_progress
            )
            
            monotonize!(paired_efdr_values)
            
            # Debug: Check paired EFDR
            max_paired = maximum(paired_efdr_values)
            println("  Max paired protein EFDR for this run: $max_paired")
            
            # Store paired EFDR results
            protein_no_decoys[run_indices, :protein_group_entrapment_fdr] = Float32.(paired_efdr_values)
        end
        
        # Remove temporary column
        select!(protein_no_decoys, Not(:_row_idx))
    else
        # Single file case
        println("Processing single protein file...")
        
        scores = protein_no_decoys[!, score_col]::AbstractVector{Float32}
        qvals = protein_no_decoys.Protein_Qvalue::AbstractVector{Float32}
        
        # Calculate combined EFDR
        combined_efdr_values = calculate_combined_efdr(
            Float64.(scores),
            entrap_labels,
            Float64.(qvals);
            r = r_lib,
            show_progress = show_progress
        )
        monotonize!(combined_efdr_values)
        protein_no_decoys[!, :combined_protein_fdr] = Float32.(combined_efdr_values)
        
        # Calculate paired EFDR
        complement_scores = get_complement_scores(scores, pairing_info.complement_indices)
        paired_efdr_values = calculate_paired_efdr(
            Float64.(scores),
            Float64.(complement_scores),
            pairing_info.is_original,
            Float64.(qvals);
            r = r_lib,
            show_progress = show_progress
        )
        
        monotonize!(paired_efdr_values)
        protein_no_decoys[!, :protein_group_entrapment_fdr] = Float32.(paired_efdr_values)
    end
    
    # Step 12: Generate outputs
    output_file = joinpath(output_dir, "protein_group_entrapment.tsv")
    CSV.write(output_file, protein_no_decoys, delim='\t')
    println("\nProtein results saved to: $output_file")
    
    # Step 13: Create visualizations and report
    # Set x-axis limit to minimum of the two thresholds
    xlim_max = min(local_qval_threshold, global_qval_threshold)
    
    # Generate paired EFDR plot (for backward compatibility)
    plot_protein_efdr_comparison(
        protein_no_decoys;
        output_path = joinpath(output_dir, "protein_efdr_comparison.pdf"),
        title = "Entrapment Analysis Protein Groups",
        xlim = (0, xlim_max)
    )
    
    # Generate comprehensive report with all plots
    generate_analysis_report(
        protein_no_decoys,
        output_dir;
        combined_efdr_col = :combined_protein_fdr,
        paired_efdr_col = :protein_group_entrapment_fdr,
        fdr_col = :Protein_Qvalue,
        xlim = (0, xlim_max)
    )
    
    return protein_no_decoys
end

# Single file convenience functions
run_efdr_analysis(parquet_file::String, library_path::String; kwargs...) = 
    run_efdr_analysis([parquet_file], library_path; kwargs...)

run_protein_efdr_analysis(parquet_file::String, library_path::String; kwargs...) = 
    run_protein_efdr_analysis([parquet_file], library_path; kwargs...)