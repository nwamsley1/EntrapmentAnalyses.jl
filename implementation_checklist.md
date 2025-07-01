# Implementation Checklist: Plex-Specific Pairing

## Phase 1: Add Type Definitions
- [x] Add struct definitions to pairing.jl
  - [x] Define PeptideKey struct
  - [x] Define PlexPairKey struct
  - [x] Define ScorePair struct
  - [x] Add hash and equality methods for structs
- [ ] Test that types work correctly

## Phase 2: Implement Helper Functions
- [x] Add init_entrapment_pairs_dict function
  - [ ] Test with sample library data
- [x] Add add_plex_complement_scores! function
  - [x] Implement first pass to populate is_original and pair_id
  - [x] Implement per-file processing loop
  - [x] Implement score dictionary building
  - [x] Implement complement score assignment
  - [ ] Test with sample data

## Phase 3: Update compute_pairing_vectors
- [x] Backup current compute_pairing_vectors function
- [x] Replace with new implementation
  - [x] Call init_entrapment_pairs_dict
  - [x] Call add_plex_complement_scores!
  - [x] Extract vectors for return value
  - [x] Include complement_scores in return tuple
- [ ] Test compatibility with existing code

## Phase 4: Update API Usage
- [x] Update run_efdr_analysis in api.jl
  - [x] Locate paired EFDR calculation section
  - [x] Replace get_complement_scores calls with pre-computed values
  - [x] Use pairing_info.complement_scores instead
- [x] Update run_protein_efdr_analysis similarly
- [ ] Test that EFDR calculations still work

## Phase 5: Testing & Validation
- [ ] Create test data with multiple plexes
- [ ] Verify complement scores are plex-specific
- [ ] Compare results with notebook output
- [ ] Run existing tests to ensure compatibility

## Phase 6: Cleanup
- [ ] Mark get_complement_scores as deprecated
- [ ] Update documentation
- [ ] Commit all changes

## Notes
- Make incremental commits after each major step
- Test frequently with small datasets
- Ask for help if any step is unclear