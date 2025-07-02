# Implementation Checklist: Refactor compute_pairing_vectors

## Phase 1: Update pairing.jl
- [x] Update compute_pairing_vectors docstring with column descriptions
- [x] Rename function to compute_pairing_vectors!
- [x] Add entrap_label column creation with descriptive message
- [x] Add complement_indices column (deprecated)
- [x] Update return statement
- [ ] Test that columns are added correctly

## Phase 2: Update run_efdr_analysis in api.jl
- [ ] Add detailed comment at Step 7 explaining columns
- [ ] Update function docstring to document added columns
- [ ] Change compute_pairing_vectors call to compute_pairing_vectors!
- [ ] Remove pairing_info variable
- [ ] Remove entrap_labels computation (use column)
- [ ] Update multi-file loop to use columns
- [ ] Update single-file case to use columns

## Phase 3: Update run_protein_efdr_analysis in api.jl
- [ ] Add detailed comment where pairing is computed
- [ ] Update function docstring
- [ ] Change compute_pairing_vectors call to compute_pairing_vectors!
- [ ] Update multi-file loop to use columns
- [ ] Update single-file case to use columns

## Phase 4: Cleanup and Testing
- [ ] Remove all references to pairing_info
- [ ] Verify no broken references
- [ ] Test with sample data
- [ ] Update CLAUDE.md documentation

## Commits
- [ ] Commit Phase 1 changes
- [ ] Commit Phase 2 changes
- [ ] Commit Phase 3 changes
- [ ] Commit final cleanup