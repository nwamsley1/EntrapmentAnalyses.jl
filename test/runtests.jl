using Test
using EntrapmentAnalyses

@testset "EntrapmentAnalyses.jl" begin
    # Unit tests
    include("unit/test_pairing.jl")
    include("unit/test_efdr_methods.jl") 
    include("unit/test_qvalue.jl")
    include("unit/test_paired_efdr_validation.jl")
    
    # Integration tests
    include("integration/test_full_pipeline.jl")
end
