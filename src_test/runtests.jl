using Base.Test
ROOT = pwd()
include(joinpath("..", "src_SOShierarchy", "SOShierarchy.jl"))

@testset "Global real tests" begin
    include("sos_example1.jl")
    include("sos_example2.jl")
    include("sos_example3_WB2_rankrel.jl")
    include("sos_example4_WB2_order2.jl")
    include("sos_example5_WB5_rankrel.jl")
end
