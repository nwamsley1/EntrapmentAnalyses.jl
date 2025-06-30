# Deep Analysis: Java vs Julia Paired FDR Implementation

## Critical Realization

The Java and Julia implementations are **NOT** doing the same thing!

## Java Implementation Analysis

The Java `FDPCalc` class:
1. Takes a single position `i` as input
2. Calculates FDR **only at that specific threshold** (position i)
3. Has **O(n) complexity for a single threshold**
4. Would need to be called n times (once per position) for full analysis, giving O(n²) total

### Java Algorithm Flow (for threshold at position i):

```java
// First loop: Count everything up to position i
for(int k=0; k<=i; k++) {
    // Count n_targets and n_entrapment_targets
    // Build all_hits set (all sequences up to position i)
}

// Second loop: Process from position i backwards to 0
for(int k=i; k>=0; k--) {
    // Build cur_hits incrementally (from worst to best)
    // For each entrapment:
    //   - If paired target in cur_hits: both > threshold (n_p_t_s++)
    //   - If paired target NOT in all_hits: only entrapment > threshold (n_p_s_t++)
}
```

**Key insight**: The Java code uses set membership to determine score relationships:
- `cur_hits`: Contains sequences from position i (worst) up to current k (better)
- `all_hits`: Contains ALL sequences from 0 to i
- If paired target is NOT in `all_hits`, it must be at position > i (worse score)

## Julia Implementation Analysis

Our Julia implementation:
1. Takes ALL positions at once
2. Calculates EFDR for **every threshold in a single pass**
3. Has **O(n²) complexity inherently** due to nested loops

### Julia Algorithm Flow:

```julia
for i in 1:n  # For each threshold position
    for j in 1:i  # For each position up to threshold
        # Count based on direct score comparisons
    end
end
```

## The Fundamental Difference

### Java: "Where is the paired sequence?"
- Uses **position-based** logic via set membership
- Determines if paired sequence is above/below threshold by checking which set it's in
- Efficiently handles one threshold at a time

### Julia: "What are the scores?"
- Uses **score-based** logic with direct comparisons
- Explicitly compares scores against threshold
- Computes all thresholds in one pass

## Are They Equivalent?

**YES**, but only if:

1. **Java** is called for every position i from 0 to n-1
2. Both are working with the same sorted order
3. The score at position i represents the same threshold in both

The Java approach is actually quite clever:
- By building sets incrementally, it avoids explicit score comparisons
- The "NOT in all_hits" check elegantly identifies sequences below threshold
- It's designed for incremental/streaming calculation

## Correcting My Previous Analysis

I was wrong to say they have the same complexity for a single run. The Java implementation is O(n) per threshold, while our Julia is O(n²) for all thresholds.

However, the **formulas are still equivalent**:
- Java: `(n_entrapment_targets + 2*n_p_t_s + n_p_s_t) / (n_entrapment_targets + n_targets)`
- Julia: `(Nε + Nεsτ + 2*Nετs) / (Nτ + Nε)`

The difference is in how they determine which category each entrapment falls into.