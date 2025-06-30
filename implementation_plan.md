# Implementation Plan: Adapt EntrapmentAnalyses Module for Notebook Data Format

## Executive Summary
This plan outlines the adaptation of the EntrapmentAnalyses.jl module to accept input data in the format used by the Jupyter notebook implementation (Parquet files for PSM results and TSV for spectral library), while maintaining the superior module structure and algorithms from Pioneer.jl.

## Goals
1. **Primary**: Enable the module to process the same data format as the notebook (Parquet + TSV)
2. **Secondary**: Maintain the module's architectural advantages (testing, documentation, extensibility)
3. **Tertiary**: Provide seamless migration path for notebook users

## Phase 1: Core Module Structure (Day 1)

### 1.1 Module Setup
```julia
# src/EntrapmentAnalyses.jl
module EntrapmentAnalyses

using DataFrames
using CSV
using Parquet2
using Plots
using Dictionaries
using Printf
using Statistics
using Dates

# Include all source files
include("io/data_loaders.jl")
include("core/notebook_pairing.jl")
include("core/efdr_methods.jl")
include("core/scoring.jl")
include("analysis/efdr_analysis.jl")
include("analysis/qvalue_calculation.jl")
include("plotting/visualization.jl")
include("api.jl")

# Export main API
export run_notebook_efdr_analysis

# Export types and methods
export NotebookEFDRMethod, CombinedEFDR, PairedEFDR
export monotonize!, getPairedEmpiricalFdr
export initEntrapPairsDict

end
```

### 1.2 Dependencies Update
Update Project.toml to include:
- Parquet2 (for reading Parquet files)
- Dictionaries (for efficient pair mappings)
- Statistics (for analysis functions)

## Phase 2: Data Loading and Format Conversion (Day 1-2)

### 2.1 Data Loaders
Create `src/io/data_loaders.jl`:

```julia
"""
Load PSM results from Parquet files
"""
function load_parquet_results(filepaths::Vector{String})
    # Implementation details:
    # 1. Load each Parquet file using Parquet2
    # 2. Concatenate DataFrames
    # 3. Return combined DataFrame
end

"""
Load spectral library from TSV file
"""
function load_tsv_library(filepath::String)
    # Implementation details:
    # 1. Read TSV using CSV.File
    # 2. Validate required columns
    # 3. Return DataFrame
end

"""
Convert notebook format to module format
"""
function convert_notebook_to_module_format!(df::DataFrame; 
                                          is_results::Bool=true)
    # Column mappings:
    # - decoy → target (invert boolean)
    # - stripped_seq → sequence
    # - PredVal → score
    # - z → charge
    # - EntrapmentGroupId → entrapment_group_id
    # - PeptideSequence → stripped_seq (for library)
    # - PrecursorCharge → prec_charge
end
```

### 2.2 Format Validation
```julia
function validate_notebook_format(df::DataFrame; is_results::Bool=true)
    required_cols = is_results ? 
        [:stripped_seq, :z, :PredVal, :decoy, :channel, :file_name] :
        [:PeptideSequence, :PrecursorCharge, :EntrapmentGroupId, :PrecursorIdx]
    
    # Check for required columns
    # Provide helpful error messages
end
```

## Phase 3: Pairing System Implementation (Day 2-3)

### 3.1 Notebook-Compatible Pairing
Create `src/core/notebook_pairing.jl`:

```julia
"""
Initialize entrapment pairs dictionaries (notebook compatible)
"""
function initEntrapPairsDict(
    lib_df::DataFrame,
    mod_seq_col::Symbol,
    charge_col::Symbol,
    entrap_group_col::Symbol,
    pair_column::Symbol
)
    # Direct port from notebook implementation
    # Returns: pair_dict, is_original_dict
end

"""
Convert notebook pairing to module format
"""
function convert_pairing_to_module_format!(
    results_df::DataFrame,
    library_df::DataFrame,
    pair_dict::Dictionary,
    is_original_dict::Dictionary
)
    # Add entrap_pair_id column
    # Map notebook pairs to module's pair system
end
```

## Phase 4: EFDR Methods (Day 3-4)

### 4.1 Notebook-Compatible EFDR
Update `src/core/efdr_methods.jl`:

```julia
"""
Direct port of notebook's getPairedEmpiricalFdr
"""
function getPairedEmpiricalFdr(
    df::AbstractDataFrame,
    score_col::Symbol,    
    channel_id_col::Symbol,
    mod_seq_col::Symbol,
    charge_col::Symbol,
    is_original_dict::Dictionary,
    pair_id_dict::Dictionary
)
    # Direct implementation from notebook
    # Includes plex-specific pairing
    # Returns empirical FDR vector
end

"""
Monotonize q-values to ensure proper FDR
"""
function monotonize!(values::AbstractVector{Float32})
    # Direct port from notebook
end
```

### 4.2 Integration with Module Structure
```julia
# Adapter to use notebook EFDR with module's struct system
struct NotebookPairedEFDR <: EFDRMethod
    df::DataFrame
    score_col::Symbol
    channel_col::Symbol
    is_original_dict::Dictionary
    pair_id_dict::Dictionary
end

function calculate_efdr(method::NotebookPairedEFDR)
    # Call getPairedEmpiricalFdr internally
end
```

## Phase 5: Analysis Pipeline (Day 4-5)

### 5.1 Q-value Calculation
Create `src/analysis/qvalue_calculation.jl`:

```julia
"""
Calculate local q-values (target-decoy approach)
"""
function calculate_local_qvalues!(df::DataFrame)
    # Sort by score (descending), targets before decoys
    # Calculate running FDR
    # Monotonize
end

"""
Calculate global q-values (best per precursor)
"""
function calculate_global_qvalues!(df::DataFrame)
    # Group by (decoy, stripped_seq, z)
    # Take best scoring per group
    # Calculate q-values
    # Map back to original dataframe
end
```

### 5.2 Main Analysis Function
Create `src/api.jl`:

```julia
"""
Run EFDR analysis using notebook data format
"""
function run_notebook_efdr_analysis(
    parquet_files::Vector{String},
    library_tsv::String;
    output_dir::String = "efdr_output",
    global_qval_threshold::Float64 = 0.01,
    analyze_proteins::Bool = true,
    plot_formats::Vector{Symbol} = [:pdf, :png]
)
    # Step 1: Load data
    # Step 2: Filter and process
    # Step 3: Calculate q-values
    # Step 4: Calculate entrapment FDR
    # Step 5: Generate outputs
    
    return analysis_results
end
```

## Phase 6: Visualization (Day 5)

### 6.1 Plotting Functions
Create `src/plotting/visualization.jl`:

```julia
"""
Create entrapment FDR comparison plot
"""
function plot_entrapment_comparison(
    df::DataFrame;
    title::String,
    output_path::String
)
    # Use notebook's color scheme: RGB(0.39215686, 0.58431373, 0.92941176)
    # X-axis: Traditional FDR
    # Y-axis: Entrapment FDR
    # Diagonal reference line
end
```

## Phase 7: Testing and Documentation (Day 6)

### 7.1 Test Suite
Create comprehensive tests:
- Unit tests for each function
- Integration test with sample data
- Comparison with notebook outputs

### 7.2 Documentation
- Update README with usage examples
- Create migration guide from notebook
- Document all public functions

## Implementation Timeline
- **Day 1**: Core module structure and data loaders
- **Day 2**: Pairing system implementation
- **Day 3**: EFDR methods
- **Day 4**: Analysis pipeline
- **Day 5**: Visualization and outputs
- **Day 6**: Testing and documentation

## Risk Mitigation
1. **Data format variations**: Implement flexible column detection
2. **Performance concerns**: Profile and optimize O(n²) operations
3. **Backward compatibility**: Maintain notebook function signatures

## Success Criteria
1. Module processes identical input files as notebook
2. Results match notebook output within numerical precision
3. Performance is equal or better than notebook
4. Clear migration path documented

## Future Enhancements
1. Support for additional file formats (mzML, mzIdentML)
2. Parallel processing for multiple files
3. Interactive visualization options
4. Integration with existing proteomics pipelines