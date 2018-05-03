using Mosek, DataStructures

include("src_SOShierarchy/SOShierarchy.jl")
include("src_SOShierarchy/SDP_types.jl")
include("src_SOShierarchy/run_mosek.jl")

instance = read_SDPInstance(".")

println("VAR_TYPES size:     $(size(instance.VAR_TYPES))")
println("BLOCKS size:        $(size(instance.BLOCKS))")
println("LINEAR size:        $(size(instance.LINEAR))")
println("CONST size:         $(size(instance.CONST))")

sdp = SDP_Problem()

set_constraints!(sdp, instance)
set_blocks!(sdp, instance)
set_matrices!(sdp, instance)
set_linear!(sdp, instance)
set_const!(sdp, instance)


for (cstr, block) in sdp.name_to_sdpblock
    println("  - $cstr -> $block")
end

for (name, ctr) in sdp.name_to_ctr
    println("  * $name \t $ctr")
end

for (name, ctr) in sdp.matrices
    println("  s $name \t $ctr")
end

for (name, ctr) in sdp.linear
    println("  l $name \t $ctr")
end

for (name, ctr) in sdp.cst_ctr
    println("  c $name \t $ctr")
end




primal=Dict{Tuple{String,String,String}, Float64}()
dual=Dict{String, Float64}()

solve_mosek(sdp::SDP_Problem, primal::Dict{Tuple{String,String,String}, Float64}, dual::Dict{String, Float64})

println()
println("Primal:")
for ((k1, k2, k3), v) in primal
    @printf("%20s, %20s, %20s -> %f\n", k1, k2, k3, v)
end

println("Dual:")
for (k, v) in dual
    @printf("%20s -> %f\n", k, v)
end