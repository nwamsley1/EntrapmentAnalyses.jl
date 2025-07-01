# Plan to Fix Safety Checks in Data Loaders

## Current Issues

### In `load_parquet_results` (lines 56-83)

1. **Incorrect warning messages**: All warnings say "Found missing decoy values" regardless of which column is being checked
2. **Inaccurate replacement descriptions**: All warnings say "replacing with empty string" even when:
   - `decoy` → `false` (boolean)
   - `z` → `0` (UInt8)
   - `PredVal` → `0` (Float32)
3. **Code duplication**: Repetitive pattern for each column
4. **Inconsistent type checking**: Only `stripped_seq` has a type check at the end

### In `load_spectral_library` (lines 145-148)

1. **Same incorrect warning message**: Says "Found missing decoy values" when checking `PrecursorCharge`
2. **Wrong replacement description**: Says "replacing with empty string" when actually replacing with `0` (UInt8)
3. **No safety checks for other columns**: Only `PrecursorCharge` is checked for missing values, but `PeptideSequence`, `EntrapmentGroupId`, and `PrecursorIdx` could also have missing values

## Proposed Solution

### Option 1: Fix Individual Warnings (Minimal Change)
- Correct each warning message to reference the actual column name
- Fix replacement descriptions to match actual values
- Add type checks for all columns that undergo transformation

### Option 2: Refactor with Helper Function (Recommended)
Create a helper function to handle missing value replacement:

```julia
function replace_missing!(df::DataFrame, col::Symbol, default_value, target_type::Type, description::String)
    if any(ismissing.(df[!, col]))
        @warn "Found missing $col values, replacing with $description"
        df[!, col] = [target_type(coalesce(x, default_value)) for x in df[!, col]]
    end
end
```

Then use it for each column in `load_parquet_results`:
```julia
replace_missing!(results_df, :stripped_seq, "", String, "empty string")
replace_missing!(results_df, :decoy, false, Bool, "false")
replace_missing!(results_df, :z, 0, UInt8, "0")
replace_missing!(results_df, :PredVal, 0, Float32, "0.0")
replace_missing!(results_df, :file_name, "", String, "empty string")
```

And in `load_spectral_library`:
```julia
replace_missing!(library_df, :PeptideSequence, "", String, "empty string")
replace_missing!(library_df, :PrecursorCharge, 0, UInt8, "0")
replace_missing!(library_df, :EntrapmentGroupId, 0, Int, "0")
replace_missing!(library_df, :PrecursorIdx, 0, Int, "0")
```

### Option 3: Configuration-Based Approach
Define column specifications and iterate:

For `load_parquet_results`:
```julia
const PARQUET_COLUMN_SPECS = [
    (col=:stripped_seq, default="", type=String, desc="empty string"),
    (col=:decoy, default=false, type=Bool, desc="false"),
    (col=:z, default=0, type=UInt8, desc="0"),
    (col=:PredVal, default=0, type=Float32, desc="0.0"),
    (col=:file_name, default="", type=String, desc="empty string")
]

for spec in PARQUET_COLUMN_SPECS
    if any(ismissing.(results_df[!, spec.col]))
        @warn "Found missing $(spec.col) values, replacing with $(spec.desc)"
        results_df[!, spec.col] = [spec.type(coalesce(x, spec.default)) for x in results_df[!, spec.col]]
    end
end
```

For `load_spectral_library`:
```julia
const LIBRARY_COLUMN_SPECS = [
    (col=:PeptideSequence, default="", type=String, desc="empty string"),
    (col=:PrecursorCharge, default=0, type=UInt8, desc="0"),
    (col=:EntrapmentGroupId, default=0, type=Int, desc="0"),
    (col=:PrecursorIdx, default=0, type=Int, desc="0")
]

for spec in LIBRARY_COLUMN_SPECS
    if any(ismissing.(library_df[!, spec.col]))
        @warn "Found missing $(spec.col) values, replacing with $(spec.desc)"
        library_df[!, spec.col] = [spec.type(coalesce(x, spec.default)) for x in library_df[!, spec.col]]
    end
end
```

## Implementation Steps

1. **Choose approach** - I recommend Option 2 for balance of clarity and DRY principle
2. **Implement the fix** in `data_loaders.jl`:
   - Add the `replace_missing!` helper function at the top of the file
   - Replace the repetitive code in both `load_parquet_results` and `load_spectral_library`
3. **Add comprehensive type checks** for all processed columns in both functions
4. **Document the behavior** in both function docstrings about missing value handling

## Additional Improvements

1. **Consider making defaults configurable** via function parameters
2. **Add option to error on missing values** instead of replacing them
3. **Log statistics** about how many missing values were replaced per column
4. **Validate column types** before processing to catch issues early

## Benefits of This Refactoring

1. **Eliminates code duplication** between the two functions
2. **Fixes all incorrect warning messages** to accurately report which column has missing values
3. **Provides consistent handling** of missing values across all data loading functions
4. **Makes the code more maintainable** - adding new columns or changing defaults is easier
5. **Improves error visibility** - users will know exactly which columns had missing values
6. **Ensures type safety** - all columns will have consistent types after loading