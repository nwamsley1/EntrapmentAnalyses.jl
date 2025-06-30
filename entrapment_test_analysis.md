# Entrapment Test Analysis - Technical Documentation

## Overview
This Jupyter notebook implements an entrapment-based false discovery rate (FDR) analysis for mass spectrometry proteomics data from a 9-plex TMT experiment. The analysis uses a paired entrapment strategy where original peptide sequences are paired with shuffled "entrapment" sequences to empirically estimate FDR.

## Key Concepts for Understanding This Code

### Entrapment Strategy
- **Original sequences**: Real peptide sequences from the proteome
- **Entrapment sequences**: Shuffled versions of original sequences that shouldn't exist in nature
- **Paired analysis**: Each entrapment sequence is paired with its corresponding original sequence
- The ratio of entrapment hits to total hits provides an empirical FDR estimate

### Data Structure
- The analysis processes multiple Parquet files containing peptide-spectrum matches (PSMs)
- Each PSM has:
  - `stripped_seq`: Peptide sequence
  - `z`: Charge state
  - `PredVal`: Prediction score (higher is better)
  - `decoy`: Boolean indicating if it's a decoy hit
  - `channel`: TMT channel (plex identifier)
  - `file_name`: Run identifier
  - `protein`: Protein group assignment

## Core Functions

### 1. `monotonize!(values)`
Converts q-values to FDR by ensuring monotonicity. Traverses the array in reverse to ensure that FDR values never decrease as you go down the ranked list.

### 2. `getPairedEmpiricalFdr()`
The heart of the entrapment analysis. This function:
- Takes a dataframe sorted by score (descending)
- For each score threshold, calculates:
  - N_e: Number of entrapment hits
  - N_t: Number of target (non-entrapment) hits
  - N_est: Entrapment hits scoring above threshold but below their paired original
  - N_ets: Entrapment hits scoring above their paired original
- Returns empirical FDR using formula: `(N_e + N_est + 2*N_ets)/(N_t + N_e)`

### 3. `initEntrapPairsDict()`
Creates two dictionaries from the spectral library:
- `pair_dict`: Maps (sequence, charge) → pair_id
- `is_original_dict`: Maps (sequence, charge) → boolean (true if original, false if entrapment)

## Analysis Workflow

### Phase 1: Library Processing
1. Loads the spectral library containing both original and entrapment sequences
2. Identifies target sequences (EntrapmentGroupId = 0) vs entrapment sequences (EntrapmentGroupId = 1)
3. Creates pairing dictionaries for the paired analysis

### Phase 2: Precursor-Level Analysis
1. **Data Loading**: Loads multiple Parquet files and combines them
2. **Decoy-based FDR**: 
   - Calculates local q-values using target-decoy approach
   - Calculates global q-values (best PSM per precursor)
   - Filters to 1% global FDR
3. **Entrapment FDR**:
   - Removes decoys from the dataset
   - Groups by file/run
   - Calculates paired empirical FDR for each run
   - Generates visualization comparing decoy-based FDR to entrapment FDR

### Phase 3: Protein Group Analysis
1. **Protein Rollup**: Selects best-scoring PSM per protein group, per channel, per run
2. **Protein FDR**: Calculates protein-level q-values using target-decoy approach
3. **Protein Entrapment FDR**: 
   - Applies paired entrapment analysis at protein level
   - Generates similar visualization for protein groups

## Key Outputs

### Files Generated:
1. `9plex_precursor_entrapment.tsv`: Precursor-level results with entrapment FDR
2. `9plex_pg_entrapment.tsv`: Protein group-level results with entrapment FDR
3. `9plex_entrapment_paired_precursors_jmod.pdf`: Scatter plot comparing FDR methods for precursors
4. `9plex_entrapment_paired_pgs_jmod.pdf`: Scatter plot comparing FDR methods for protein groups

### Visualization Details:
- X-axis: Traditional target-decoy FDR
- Y-axis: Entrapment-based FDR
- Diagonal line: Perfect agreement between methods
- Blue lines: Individual runs with 75% transparency
- Points above diagonal indicate entrapment FDR > decoy FDR (more conservative)

## Important Implementation Notes

1. **Missing Values**: The code handles missing values in FDR calculations (see monotonize! overload)
2. **Tie Breaking**: When scores are tied, targets are ranked above decoys
3. **Multiple Runs**: Each run is processed separately for FDR calculation
4. **Channel Handling**: If no channel information exists, dummy channel 0 is added
5. **Color Scheme**: Uses consistent blue color (RGB(0.39215686, 0.58431373, 0.92941176)) for all runs

## Dependencies
- CSV, DataFrames, Arrow, Plots, Test, Dictionaries
- Parquet2 (for reading Parquet files)
- Custom TMT data in specific directory structure

## Performance Considerations
- Loads all data into memory before processing
- Uses dictionaries for O(1) lookup of pair relationships
- Processes each run independently for parallelization potential