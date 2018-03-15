ROOT = pwd()
include(joinpath(ROOT, "src_PowSysMod", "PowSysMod_body.jl"))
include("momentmatrix_dense.jl")

#######################################
# Example
v1 = Variable("v1", Complex)
v2 = Variable("v2", Complex)
vars = Set([v1,v2])

println("⚇ vars: $vars")

d, k = 2, 1
p = abs2(v1) - 0.9^2

mm = MomentMatrix(vars, d-k)
println("\n⚇ Moment matrix: \n$mm")

println("⚇ Localizing matrix wrt p = $p:\n$(mm*p)")

mmb = convertMMtobase(mm*p, d, k)

mmb.basis[Exponent(Dict(v1=>Degree(1,1)))]
