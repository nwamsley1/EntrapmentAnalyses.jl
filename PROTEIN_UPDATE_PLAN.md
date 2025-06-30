# Protein-Level Analysis Update Plan

## Overview
Update `run_protein_efdr_analysis` to match all the improvements made to the precursor-level analysis.

## Step-by-Step Implementation Plan

### Step 1: Update Q-value Calculation
- [x] Apply `calculate_qvalues_per_file!` to PSM data before protein rollup
- [x] Ensure protein-level q-values are calculated per file after rollup
- [x] Remove any global q-value calculations

### Step 2: Add Local Q-value Filtering
- [x] Add `local_qval_threshold` parameter (default 0.05)
- [x] Filter PSMs by local q-value before protein rollup
- [x] Add `global_qval_threshold` parameter for consistency

### Step 3: Replace Dictionary-Based Approach with Vector-Based
- [x] Remove `create_pairing_dictionaries` function call
- [x] Use `compute_pairing_vectors` instead for efficiency
- [x] Remove old `calculate_protein_efdr!` and `build_protein_pairing_vectors`

### Step 4: Implement Per-File EFDR Calculation
- [x] Group protein data by file (if multiple files)
- [x] Calculate both combined and paired EFDR per file
- [x] Store results using parent DataFrame indices (not SubDataFrame)

### Step 5: Add Combined EFDR Support
- [x] Add combined EFDR calculation to protein analysis
- [x] Initialize both EFDR columns (combined and paired)
- [x] Calculate both methods in the same loop

### Step 6: Update Visualization
- [x] Add `xlim` parameter based on min(local_qval, global_qval)
- [x] Pass xlim to plot functions
- [x] Plots already have dynamic y-axis scaling

### Step 7: Add Report Generation
- [x] Call `generate_analysis_report` with protein-level data
- [x] Ensure proper column names for protein-level analysis
- [x] Generate comprehensive markdown report

### Step 8: Code Structure Changes

#### Current Structure:
```julia
function run_protein_efdr_analysis(...)
    # Load data
    # Create pairing dictionaries (OLD)
    # Prepare protein data
    # Calculate protein EFDR (single method)
    # Save results
    # Create single plot
end
```

#### New Structure:
```julia
function run_protein_efdr_analysis(...)
    # Load data
    # Calculate per-file q-values
    # Filter by local q-value threshold
    # Compute pairing vectors (NEW)
    # Prepare protein data with per-file grouping
    # For each file:
        # Calculate combined EFDR
        # Calculate paired EFDR
    # Save results
    # Generate comprehensive report with all plots
end
```

### Step 9: Function Signature Updates
- [ ] Add `local_qval_threshold` parameter
- [ ] Add `global_qval_threshold` parameter
- [ ] Ensure backward compatibility

### Step 10: Testing
- [ ] Test with single file
- [ ] Test with multiple files
- [ ] Verify both EFDR methods produce non-zero results
- [ ] Check visualization outputs
- [ ] Validate markdown report generation

## Key Considerations

1. **Protein Rollup Timing**: Ensure q-value filtering happens at PSM level before protein rollup
2. **Pairing Consistency**: Protein sequences must map back to library pairs correctly
3. **Performance**: Protein-level analysis typically has fewer entries, so O(n²) is less of a concern
4. **Column Names**: Protein analysis uses different column names (e.g., `Protein_Qvalue` vs `local_qvalue`)

## Dependencies
- Ensure `protein_analysis.jl` functions are compatible with per-file processing
- Update any helper functions that assume global processing
- Maintain backward compatibility where possible