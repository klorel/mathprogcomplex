using Base.Test
ROOT = pwd()
include(joinpath("..", "src_SOShierarchy", "SOShierarchy.jl"))

@testset "Global real tests" begin
    include("sos_example1.jl")
    include("sos_example2.jl")
end
