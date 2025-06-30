# CRITICAL FINDINGS: Algorithm Differences and Implementation Bug

## The Three Implementations

### 1. **Notebook Implementation** (Reference)
```julia
if (entrapment_score >= s) & (s > original_score)
    N_est += 1
elseif (entrapment_score > original_score) & (original_score >= s)
    N_ets += 1
end
```
- `N_est`: Entrapment ≥ threshold > original
- `N_ets`: Entrapment > original AND original ≥ threshold (entrapment must WIN)

### 2. **Java Implementation**
- `n_p_s_t`: Entrapment above threshold, paired target below (position-based)
- `n_p_t_s`: Both above threshold (position-based)
- **Does NOT check who wins!**

### 3. **Our Julia Module** (Has a bug!)
```julia
if score_rels[idx] == :entrap_wins
    if o_score >= s
        Nετs += 1  # Correct: entrap wins and both above
    end
elseif score_rels[idx] != :unpaired && scores[idx] >= s
    if o_score < s
        Nεsτ += 1  # BUG: Only counts if entrap DOESN'T win!
    end
end
```

## The Key Differences

### Notebook vs Java:
The notebook's paired FDR method is **more sophisticated**:
- It only counts entrapments that actually WIN in the `N_ets` category
- Java counts ANY pair where both are above threshold

This leads to different formulas effectively:
- **Notebook**: (N_e + N_est + 2*N_ets) / (N_t + N_e) where N_ets requires winning
- **Java**: (N_e + N_est + 2*N_ets) / (N_t + N_e) where N_ets is just position-based

### Our Implementation Bug:
We incorrectly implemented the notebook's logic. The fix should be:

```julia
if is_original[idx]
    Nτ += 1
else
    Nε += 1
    
    # Check N_est condition (entrap >= s > original)
    if scores[idx] >= s && complement_scores[idx] < s
        Nεsτ += 1
    end
    
    # Check N_ets condition (entrap > original AND original >= s)
    if score_rels[idx] == :entrap_wins && complement_scores[idx] >= s
        Nετs += 1
    end
end
```

## Impact

1. **Java vs Notebook**: Java will generally give HIGHER FDR estimates because it counts more pairs in the penalty categories (doesn't require entrapment to win)

2. **Our Bug**: We're UNDER-counting Nεsτ, which would lead to LOWER FDR estimates than intended

3. **Which is "correct"?**: The notebook implementation appears to follow the Noble lab paper more precisely, requiring entrapments to actually outperform their paired targets to be counted as false positives.

## Action Items

1. **Fix our implementation** to match the notebook exactly
2. **Document** that our implementation follows the notebook (Noble lab) version, NOT the Java version
3. **Add tests** to verify the specific counting logic