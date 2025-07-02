# Plan to Fix Data Loader Warnings and Code Duplication

## Current Issues in `load_parquet_results`

1. **Incorrect warning messages**: All warnings say "Found missing decoy values" regardless of the column
2. **Inconsistent replacement descriptions**: All say "replacing with empty string" even when using different types
3. **Code duplication**: Same pattern repeated for each column
4. **Missing validations**: Only `stripped_seq` has a type check at the end

## Current Issues in `load_spectral_library`

1. **Same incorrect warning**: Says "Found missing decoy values" when checking `PrecursorCharge`
2. **No checks for other columns**: Only `PrecursorCharge` is checked, missing `PeptideSequence`, `EntrapmentGroupId`, `PrecursorIdx`

## Proposed Solution

### Step 1: Create a Helper Function

```julia
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
    
    if n_missing > 0
        @warn "Found $n_missing missing values in column '$col', replacing with $description"
        df[!, col] = [target_type(coalesce(x, default_value)) for x in df[!, col]]
    end
    
    # Validate final type
    actual_type = eltype(df[!, col])
    if actual_type != target_type && !(actual_type <: target_type)
        error("Column $col has type $actual_type, expected $target_type")
    end
    
    return n_missing
end
```

### Step 2: Define Column Specifications

```julia
# For load_parquet_results
const PARQUET_COLUMN_SPECS = [
    (col=:stripped_seq, default="", type=String, desc="empty string"),
    (col=:decoy, default=false, type=Bool, desc="false"),
    (col=:z, default=0, type=UInt8, desc="0"),
    (col=:PredVal, default=0.0f0, type=Float32, desc="0.0"),
    (col=:file_name, default="", type=String, desc="empty string")
]

# For load_spectral_library
const LIBRARY_COLUMN_SPECS = [
    (col=:PeptideSequence, default="", type=String, desc="empty string"),
    (col=:PrecursorCharge, default=0, type=UInt8, desc="0"),
    (col=:EntrapmentGroupId, default=0, type=Int, desc="0"),
    (col=:PrecursorIdx, default=0, type=Int, desc="0")
]
```

### Step 3: Update `load_parquet_results`

Replace the repetitive code with:

```julia
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
```

### Step 4: Update `load_spectral_library`

Add similar handling after loading the library:

```julia
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
```

## Alternative Approach: More Flexible Configuration

```julia
struct ColumnSpec
    col::Symbol
    default::Any
    type::Type
    description::String
    required::Bool  # Whether column must exist
    
    ColumnSpec(col, default, type, desc; required=true) = 
        new(col, default, type, desc, required)
end

# Then use it like:
const PARQUET_SPECS = [
    ColumnSpec(:stripped_seq, "", String, "empty string"),
    ColumnSpec(:decoy, false, Bool, "false"),
    ColumnSpec(:z, UInt8(0), UInt8, "0"),
    ColumnSpec(:PredVal, 0.0f0, Float32, "0.0"),
    ColumnSpec(:file_name, "", String, "empty string"),
    ColumnSpec(:protein, "", String, "empty string", required=false)  # Optional column
]
```

## Benefits

1. **DRY (Don't Repeat Yourself)**: Single function handles all columns
2. **Accurate warnings**: Each column gets its own specific warning
3. **Type safety**: Validates types after replacement
4. **Maintainability**: Easy to add/modify columns
5. **Consistency**: Same approach for both loaders
6. **Better logging**: Shows count of missing values

## Implementation Steps

1. Add the helper function to `data_loaders.jl`
2. Define column specifications as constants
3. Replace repetitive code in `load_parquet_results`
4. Add missing value handling to `load_spectral_library`
5. Test with data containing missing values
6. Update documentation if needed

## Testing Strategy

Create test data with:
- All columns present, no missing values
- Some missing values in each column
- Entire columns missing (for optional columns)
- Wrong types that need conversion

Verify:
- Correct warning messages
- Proper type conversion
- Missing value counts
- Error handling for required columns