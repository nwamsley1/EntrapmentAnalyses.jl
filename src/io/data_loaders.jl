# Data loading functions for Parquet and TSV files

# Column specifications for data validation and missing value handling
const PARQUET_COLUMN_SPECS = [
    (col=:stripped_seq, default="", type=String, desc="empty string"),
    (col=:decoy, default=false, type=Bool, desc="false"),
    (col=:z, default=0, type=UInt8, desc="0"),
    (col=:PredVal, default=0.0f0, type=Float32, desc="0.0"),
    (col=:file_name, default="", type=String, desc="empty string")
]

const LIBRARY_COLUMN_SPECS = [
    (col=:PeptideSequence, default="", type=AbstractString, desc="empty string"),
    (col=:PrecursorCharge, default=0, type=UInt8, desc="0"),
    (col=:EntrapmentGroupId, default=0, type=Int, desc="0"),
    (col=:PrecursorIdx, default=0, type=Int, desc="0")
]

"""
    handle_missing_values!(df::DataFrame, col::Symbol, default_value, target_type::Type, description::String)

Replace missing values in a DataFrame column with a default value and ensure correct type.

# Arguments
- `df`: DataFrame to modify
- `col`: Column symbol to check
- `default_value`: Value to use for missing entries
- `target_type`: Type to convert to
- `description`: Human-readable description of the default value

# Returns
- Number of missing values that were replaced
"""
function handle_missing_values!(df::DataFrame, col::Symbol, default_value, target_type::Type, description::String)
    if !hasproperty(df, col)
        error("Column $col not found in DataFrame")
    end
    
    n_missing = count(ismissing, df[!, col])
    
    if Missing <: eltype(df[!, col])
        @warn "Found $n_missing missing values in column '$col', replacing with $description"
        df[!, col] = [target_type(coalesce(x, default_value)) for x in df[!, col]]
    end
    
    # Validate final type
    actual_type = eltype(df[!, col])
    if actual_type != target_type && !(actual_type <: target_type)
        df[!, col] = convert(Vector{target_type}, df[!, col])
        #error("Column $col has type $actual_type, expected $target_type")
    end
    
    return n_missing
end

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
        if !endswith(filepath, ".parquet")
            continue
        end
        println("Loading file $i/$(length(filepaths)): $(basename(filepath))")
        df = DataFrame(Dataset(filepath); copycols=true)
        push!(dfs, df)
    end
    
    # Combine dataframes
    results_df = vcat(dfs...)
    
    # Validate required columns
    required_cols = [:stripped_seq, :z, :PredVal, :decoy, :file_name]
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
    
    # Handle missing values for all required columns
    println("\nChecking for missing values...")
    total_missing = 0
    for spec in PARQUET_COLUMN_SPECS
        n_missing = handle_missing_values!(
            results_df, 
            spec.col, 
            spec.default, 
            spec.type, 
            spec.desc
        )
        total_missing += n_missing
    end
    
    if total_missing > 0
        println("Replaced $total_missing total missing values across all columns")
    end
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

"""
    load_parquet(filepath::String)

Load a single Parquet file containing PSM results.

# Arguments
- `filepath`: Path to the Parquet file

# Returns
- DataFrame with PSM results
"""
function load_parquet(filepath::String)
    if !isfile(filepath)
        error("File not found: $filepath")
    end
    
    println("Loading Parquet file: $(basename(filepath))")
    df = DataFrame(Dataset(filepath); copycols=true)
    println("Loaded $(nrow(df)) PSMs")
    
    return df
end

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
    
    # Handle missing values for all required columns
    println("\nChecking for missing values in library...")
    total_missing = 0
    for spec in LIBRARY_COLUMN_SPECS
        n_missing = handle_missing_values!(
            library_df, 
            spec.col, 
            spec.default, 
            spec.type, 
            spec.desc
        )
        total_missing += n_missing
    end
    
    if total_missing > 0
        println("Replaced $total_missing total missing values in library")
    end

    return library_df
end