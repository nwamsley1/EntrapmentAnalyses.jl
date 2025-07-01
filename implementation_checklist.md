# Implementation Checklist: Plex-Specific Pairing

## Phase 1: Add Type Definitions
- [ ] Add struct definitions to pairing.jl
  - [ ] Define PeptideKey struct
  - [ ] Define PlexPairKey struct
  - [ ] Define ScorePair struct
  - [ ] Add hash and equality methods for structs
- [ ] Test that types work correctly

## Phase 2: Implement Helper Functions
- [ ] Add init_entrapment_pairs_dict function
  - [ ] Test with sample library data
- [ ] Add add_plex_complement_scores! function
  - [ ] Implement first pass to populate is_original and pair_id
  - [ ] Implement per-file processing loop
  - [ ] Implement score dictionary building
  - [ ] Implement complement score assignment
  - [ ] Test with sample data

## Phase 3: Update compute_pairing_vectors
- [ ] Backup current compute_pairing_vectors function
- [ ] Replace with new implementation
  - [ ] Call init_entrapment_pairs_dict
  - [ ] Call add_plex_complement_scores!
  - [ ] Extract vectors for return value
  - [ ] Include complement_scores in return tuple
- [ ] Test compatibility with existing code

## Phase 4: Update API Usage
- [ ] Update run_efdr_analysis in api.jl
  - [ ] Locate paired EFDR calculation section
  - [ ] Replace get_complement_scores calls with pre-computed values
  - [ ] Use pairing_info.complement_scores instead
- [ ] Update run_protein_efdr_analysis similarly
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