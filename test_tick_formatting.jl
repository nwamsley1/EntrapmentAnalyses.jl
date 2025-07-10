using EntrapmentAnalyses
using DataFrames

# Create test data with various ranges
test_data_small = DataFrame(
    local_qvalue = Float32.(range(0, 0.05, length=100)),
    combined_entrapment_fdr = Float32.(range(0, 0.045, length=100) .+ 0.001 * randn(100)),
    paired_entrapment_fdr = Float32.(range(0, 0.048, length=100) .+ 0.001 * randn(100)),
    decoy = fill(false, 100),
    entrapment_group = rand(Bool, 100)
)

test_data_large = DataFrame(
    local_qvalue = Float32.(range(0, 0.5, length=100)),
    combined_entrapment_fdr = Float32.(range(0, 0.45, length=100) .+ 0.01 * randn(100)),
    paired_entrapment_fdr = Float32.(range(0, 0.48, length=100) .+ 0.01 * randn(100)),
    decoy = fill(false, 100),
    entrapment_group = rand(Bool, 100)
)

# Test with small range data
println("Testing with small range data (0-0.05)...")
plot_efdr_comparison_both_methods(test_data_small, output_path="test_small_range.pdf")

# Test with larger range data
println("\nTesting with larger range data (0-0.5)...")
plot_efdr_comparison_both_methods(test_data_large, output_path="test_large_range.pdf")

println("\nPlots saved. Check test_small_range.pdf and test_large_range.pdf")