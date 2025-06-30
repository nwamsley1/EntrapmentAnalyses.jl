# Pioneer.jl EntrapmentAnalysis Module - Technical Documentation

## Overview
The Pioneer.jl EntrapmentAnalysis module is a production-ready Julia package for empirical false discovery rate (eFDR) analysis using entrapment sequences in proteomics data. This is a modularized, object-oriented implementation of the entrapment analysis concepts demonstrated in the Jupyter notebook.

## Architecture Overview

### Module Structure
```
EntrapmentAnalysis.jl         # Main module definition and exports
├── src/
│   ├── api.jl               # High-level API function (run_efdr_analysis)
│   ├── core/               
│   │   ├── efdr_methods.jl  # EFDR calculation implementations
│   │   ├── entrapment_pairing.jl  # Pairing logic
│   │   └── scoring.jl       # Score manipulation utilities
│   ├── analysis/
│   │   ├── efdr_analysis.jl # Comparison and evaluation functions
│   │   └── calibration.jl   # Calibration error calculations
│   └── plotting/
│       └── efdr_plots.jl    # Visualization functions
```

## Key Differences from Notebook Implementation

### 1. **Data Input Format**
- **Notebook**: Reads Parquet files containing PSM results and TSV spectral library
- **Module**: Expects Arrow format files for both precursor results and library
- **Rationale**: Arrow provides better performance and type safety for columnar data

### 2. **Object-Oriented Design**
- **Notebook**: Functional approach with standalone functions
- **Module**: Introduces abstract types and structs for EFDR methods
  ```julia
  abstract type EFDRMethod end
  struct CombinedEFDR{T<:Real, I<:Integer} <: EFDRMethod
  struct PairedEFDR{T<:Real, I<:Integer} <: EFDRMethod
  ```
- **Benefits**: Type safety, multiple dispatch, extensibility

### 3. **Pairing Algorithm**
- **Notebook**: Uses dictionaries to map sequences to pair IDs
  ```julia
  pair_dict[(mod_seq = mod_seq, z = z)] = pair_id
  ```
- **Module**: More sophisticated pairing based on:
  - base_pep_id (base peptide identifier)
  - prec_charge
  - is_decoy
  - mod_key (modification pattern)
- **Key Difference**: Module handles multiple entrapment groups per base peptide

### 4. **Modification Handling**
- **Notebook**: Works directly with peptide sequences
- **Module**: Introduces `getModKey()` function to extract and normalize modifications
  ```julia
  getModKey("(5,M,Unimod:4)(5,M,Unimod:35)") # => "Unimod:4;Unimod:35"
  ```

### 5. **EFDR Calculation Methods**

#### Combined EFDR
Both implementations use the formula: `(Nε*(1 + 1/r)) / (Nε + Nτ)`
- Notebook calls it "empirical FDR"
- Module implements as `CombinedEFDR` struct with `calculate_efdr()` method

#### Paired EFDR
Formula: `(Nε + Nεsτ + 2*Nετs) / (Nε + Nτ)`
- **Notebook**: Inline implementation in `getPairedEmpiricalFdr()`
- **Module**: Structured as `PairedEFDR` with cleaner variable names
- Both use O(n²) algorithm but module is more readable

### 6. **Data Processing Pipeline**

**Notebook Pipeline**:
1. Load Parquet files
2. Filter decoys
3. Calculate local/global q-values
4. Group by file/channel
5. Calculate entrapment FDR per group
6. Generate plots

**Module Pipeline**:
1. Load Arrow files
2. Filter non-targets (if target column exists)
3. Add mod_key column
4. Assign entrapment pairs
5. Add pair IDs to results
6. Add original target scores
7. Calculate EFDR for multiple score/q-value pairs
8. Generate comprehensive analysis report

### 7. **Output Generation**
- **Notebook**: 
  - TSV files with results
  - PDF plots
  - Manual plot generation
- **Module**: 
  - Arrow format results
  - Multiple plot formats (PNG, PDF)
  - Automated markdown report with embedded visualizations
  - Calibration analysis

### 8. **Error Handling and Validation**
- **Notebook**: Minimal error checking
- **Module**: 
  - Column existence validation
  - Vector length verification in constructors
  - Warning messages for missing columns
  - Graceful handling of missing data

### 9. **Flexibility and Configuration**
- **Notebook**: Hard-coded for 9-plex TMT analysis
- **Module**: 
  - Configurable EFDR methods
  - Multiple score/q-value pairs
  - Adjustable r_lib parameter
  - Customizable output formats

### 10. **Performance Optimizations**
- **Module** uses:
  - Type-stable structs with parametric types
  - Pre-allocated arrays
  - Efficient sorting with `sortperm`
  - Arrow format for faster I/O

## API Comparison

### Notebook Function
```julia
getPairedEmpiricalFdr(df, score_col, channel_id_col, mod_seq_col, 
                      charge_col, is_original_dict, pair_id_dict)
```

### Module API
```julia
run_efdr_analysis(prec_results_path, library_precursors_path;
                  output_dir="efdr_out",
                  method_types=[CombinedEFDR, PairedEFDR],
                  score_qval_pairs=[(:global_prob, :global_qval)],
                  r_lib=1.0,
                  plot_formats=[:png, :pdf],
                  verbose=true)
```

## Key Improvements in Module

1. **Automated Report Generation**: Creates markdown reports with analysis summary
2. **Calibration Analysis**: Calculates mean absolute calibration errors
3. **Multiple Method Support**: Can run multiple EFDR methods in parallel
4. **Better Visualization**: Automated plot generation with multiple formats
5. **Extensibility**: Easy to add new EFDR methods by extending abstract type
6. **Testing Infrastructure**: Includes test suite (not present in notebook)
7. **Documentation**: Comprehensive docstrings and examples

## Usage Pattern Differences

### Notebook (Interactive)
- Load data manually
- Run analysis step by step
- Manually create visualizations
- Copy/paste functions as needed

### Module (Production)
- Single function call for complete analysis
- Automated output generation
- Reusable across projects
- Version controlled and tested

## Technical Notes for Implementation

1. **Entrapment Group Handling**: Module supports multiple entrapment groups (IDs > 0), while notebook assumes binary (0 or 1)

2. **Sorting Strategy**: Both use q-value primary, score secondary sorting, but module implements via `sortperm` for efficiency

3. **Missing Data**: Module handles missing values more gracefully with Union types

4. **Column Naming**: Module uses systematic naming for EFDR columns: `{score_col}_{method}_efdr`

5. **Decoy Filtering**: Module checks for 'target' column instead of 'decoy' column (inverse logic)

## Migration Guide

To migrate from notebook to module:
1. Convert data files from Parquet to Arrow format
2. Ensure library has required columns: base_pep_id, prec_charge, is_decoy, structural_mods
3. Replace manual function calls with `run_efdr_analysis()`
4. Update column references (e.g., 'decoy' → 'target')
5. Adapt to new output structure (markdown report + plots)