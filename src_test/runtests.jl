using Base.Test
ROOT = pwd()
include(joinpath("..", "src_SOShierarchy", "SOShierarchy.jl"))

@testset "Global real tests" begin
    include("sos_example1.jl")
    include("sos_example2.jl")
    include("sos_example3.1_dense_WB2.jl")
    include("sos_example3.2_dense_WB2.jl")
    include("sos_example4.1_dense_WB5.jl")
    include("sos_example4.2_dense_WB5.jl")
    include("sos_example5_sparse_WB5_rankrel.jl")
    # include("sos_example6_matpower_rankrel.jl")
    # include("sos_example7_matpower_ordre2.jl")
end