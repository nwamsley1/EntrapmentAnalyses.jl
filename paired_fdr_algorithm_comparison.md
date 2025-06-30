# Paired FDR Algorithm Comparison: Java vs Julia Implementation

## Overview

Both implementations calculate paired empirical FDR, but there are some key differences in their approaches and formulas.

## Key Algorithm Differences

### 1. **Terminology and Variable Names**

| Java | Julia | Description |
|------|-------|-------------|
| `n_entrapment_targets` | `Nε` | Number of entrapment hits |
| `n_targets` | `Nτ` | Number of target (original) hits |
| `n_p_t_s` | `Nετs` | Entrapments scoring above threshold where paired target also scores above threshold |
| `n_p_s_t` | `Nεsτ` | Entrapments scoring above threshold where paired target scores below threshold |

### 2. **Core Formula Differences**

#### Java Implementation (Line 211):
```java
fdp_1b = (n_entrapment_targets + 2.0 * n_p_t_s + n_p_s_t) / (n_entrapment_targets + n_targets)
```

#### Julia Implementation:
```julia
efdr = (Nε + Nεsτ + 2*Nετs) / (Nτ + Nε)
```

**These formulas are mathematically equivalent!** The difference is just in variable naming:
- Java's `n_entrapment_targets` = Julia's `Nε`
- Java's `n_p_t_s` = Julia's `Nετs` 
- Java's `n_p_s_t` = Julia's `Nεsτ`

### 3. **Counting Logic Differences**

#### Java Approach:
The Java code processes matches from **low to high confidence** (line 99: `for(int k=i;k>=0;k--)`), building two sets:
- `cur_hits`: Sequences seen so far (from worst to current score)
- `all_hits`: All sequences in the dataset

For each entrapment:
- If paired target is in `cur_hits`: increment `n_p_t_s` (both > threshold)
- If paired target is NOT in `all_hits`: increment `n_p_s_t` (entrapment > threshold > target)

#### Julia Approach:
The Julia code processes matches in **sorted order by q-value**, and for each threshold `s`:
- For each entrapment with score ≥ s:
  - If entrapment wins (e > t) AND target score ≥ s: increment `Nετs`
  - If entrapment score ≥ s AND target score < s: increment `Nεsτ`

### 4. **Key Algorithmic Difference**

The main difference is in how they determine the relationship between scores:

**Java**: Uses set membership to determine if paired sequences are above/below threshold
- Builds sets incrementally from worst to best scores
- Uses "not in all_hits" to mean the paired target scored below threshold

**Julia**: Directly compares scores and uses pre-computed relationships
- Pre-computes whether entrapment wins/loses against its pair
- Explicitly checks if scores are above/below threshold `s`

### 5. **Handling of Missing Pairs**

**Java**: 
- Only counts if paired sequence exists in the dataset
- Uses `!all_hits.contains(paired_target)` to identify pairs below threshold

**Julia**:
- Explicitly handles missing pairs with `complement_scores[i] = -1`
- Pre-computes complement indices, making missing pairs explicit

### 6. **Performance Characteristics**

**Java**: 
- O(n²) time complexity
- Uses HashSets for O(1) lookup
- Processes in reverse score order

**Julia**:
- O(n²) time complexity  
- Pre-computes relationships for efficiency
- Uses vector operations where possible
- Includes progress monitoring

## Conclusion

**The algorithms are functionally equivalent** in terms of the paired FDR formula they implement. The differences are in:

1. **Implementation approach**: Java uses set membership, Julia uses direct score comparison
2. **Processing order**: Java goes from worst to best, Julia follows q-value order
3. **Missing pair handling**: Julia is more explicit about missing pairs
4. **Performance optimizations**: Julia pre-computes relationships, Java uses hash sets

Both correctly implement the paired empirical FDR method as described in the literature, just with different coding strategies.