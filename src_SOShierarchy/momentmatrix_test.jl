ROOT = pwd()
include(joinpath(ROOT, "src_PowSysMod", "PowSysMod_body.jl"))
include("momentmatrix_dense.jl")

#######################################
# Example
v1 = Variable("v1", Complex)
v2 = Variable("v2", Complex)
vars = Set([v1,v2])

println("⚇ Variables: $vars")

d, k = 2, 1
p = abs2(v1) - 0.9^2

mm = MomentMatrix(vars, d-k)
println("\n⚇ Moment matrix of order d-k=$d-$k: \n$mm")

println("⚇ Localizing matrix wrt p = $p:\n$(mm*p)\n")

println("⚇ Converting localizing matrix to the moment basis :")
mmb = convertMMtobase(mm*p, d, k)

for (moment, matrix) in mmb.basis
    println("$moment\t ⟶  $matrix")
end
