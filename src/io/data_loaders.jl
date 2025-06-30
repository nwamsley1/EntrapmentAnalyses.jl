# Data loading functions for Parquet and TSV files

"""
    load_parquet_results(filepaths::Vector{String})

Load and combine multiple Parquet files containing PSM results.

# Expected columns:
- stripped_seq: String - Peptide sequence without modifications
- z: Integer - Charge state
- PredVal: Float - Prediction score (higher is better)
- decoy: Boolean - Decoy status
- channel: Integer - Multiplex channel (optional, will add dummy if missing)
- file_name: String - Source file identifier
- protein: String - Protein group (for protein-level analysis)

# Returns
- DataFrame with combined results
"""
function load_parquet_results(filepaths::Vector{String})
    # Validate files exist
    for filepath in filepaths
        if !isfile(filepath)
            error("File not found: $filepath")
        end
    end
    
    # Load all files
    dfs = DataFrame[]
    
    for (i, filepath) in enumerate(filepaths)
        println("Loading file $i/$(length(filepaths)): $(basename(filepath))")
        df = DataFrame(Dataset(filepath); copycols=true)
        push!(dfs, df)
    end
    
    # Combine dataframes
    results_df = vcat(dfs...)
    
    # Validate required columns
    required_cols = [:stripped_seq, :z, :PredVal, :decoy]
    missing_cols = setdiff(required_cols, propertynames(results_df))
    
    if !isempty(missing_cols)
        error("Missing required columns: $(join(missing_cols, ", "))")
    end
    
    # Add dummy channel if missing
    if !hasproperty(results_df, :channel)
        println("No channel column detected, adding dummy channel 0")
        results_df[!, :channel] = zero(UInt8)
    end
    
    println("Loaded $(nrow(results_df)) total PSMs")
    
    return results_df
end

"""
    load_spectral_library(filepath::String)

Load spectral library from TSV file.

# Expected columns:
- PeptideSequence: String - Modified peptide sequence
- PrecursorCharge: Integer - Charge state
- EntrapmentGroupId: Integer - 0 for originals, >0 for entrapments
- PrecursorIdx: Integer - Unique pair identifier

# Returns
- DataFrame with library information
"""
function load_spectral_library(filepath::String)
    if !isfile(filepath)
        error("Library file not found: $filepath")
    end
    
    println("Loading spectral library from: $(basename(filepath))")
    library_df = DataFrame(CSV.File(filepath))
    
    # Validate required columns
    required_cols = [:PeptideSequence, :PrecursorCharge, :EntrapmentGroupId, :PrecursorIdx]
    missing_cols = setdiff(required_cols, propertynames(library_df))
    
    if !isempty(missing_cols)
        error("Missing required columns in library: $(join(missing_cols, ", "))")
    end
    
    # Summary statistics
    n_targets = sum(library_df.EntrapmentGroupId .== 0)
    n_entrapments = sum(library_df.EntrapmentGroupId .> 0)
    
    println("Library loaded: $n_targets targets, $n_entrapments entrapments")
    
    return library_df
end