using Mosek, DataStructures

include("src_SOShierarchy/SDP_types.jl")
include("src_SOShierarchy/run_mosek.jl")

instance = read_SDPInstance("body.sdp", "types.sdp", "rhs.sdp")

println("BLOCKS size:        $(size(instance.BLOCKS))")
println("VAR_TYPES size:     $(size(instance.VAR_TYPES))")
println("CONSTRAINTS size:   $(size(instance.CONSTRAINTS))")

sdp = SDP_Problem()

set_constraints!(sdp, instance)
set_blocks!(sdp, instance)
set_matrices!(sdp, instance)


for (cstr, block) in sdp.name_to_block
    println("  - $cstr -> $block")
end

for (name, ctr) in sdp.name_to_ctr
    println("  * $name \t $ctr")
end

for (name, ctr) in sdp.matrices
    println("  x $name \t $ctr")
end

@show sdp.linear
@show sdp.cst_ctr



primal=Dict{Tuple{String,String,String}, Float64}()
dual=Dict{String, Float64}()

solve_mosek(sdp::SDP_Problem, primal::Dict{Tuple{String,String,String}, Float64}, dual::Dict{String, Float64})

println("Primal:")
for ((k1, k2, k3), v) in primal
    @printf("%10s, %10s, %10s -> %f", k1, k2, k3, v)
end

println("Dual:")
for (k, v) in dual
    @printf("%10s -> %f", k, v)
end