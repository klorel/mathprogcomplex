using JuMP, SCS

ROOT = pwd()
include(joinpath(ROOT, "src_PowSysMod", "PowSysMod_body.jl"))
include(joinpath(ROOT, "src_SOShierarchy", "func_definitions.jl"))

########################################
# Construction du problème type
z1 = Variable("z1", Complex)
z2 = Variable("z2", Complex)
problemraw = Problem()
add_variable!(problemraw, z1); add_variable!(problemraw, z2)
set_objective!(problemraw, 3-abs2(z1)-0.5im*z1*conj(z2)^2+0.5im*z2^2*conj(z1))
add_constraint!(problemraw, "eq1", (abs2(z1) - 0.25*z1^2 - 0.25*conj(z1)^2) == 1)
add_constraint!(problemraw, "eq2", (abs2(z1) + abs2(z2)) == 3)
add_constraint!(problemraw, "eq3", (im*z2 - im*conj(z2)) == 0)
add_constraint!(problemraw, "ineq", (z2 + conj(z2)) >> 0)

print(problemraw)

########################################
# Normalizing pb and setting relaxation order by constraint
problem = normalize_problem(problemraw);

print(problem)
# ▶ variables: z
# ▶ objective: (0.5)*z + (0.5)*conj(z)
# ▶ constraints:
#    ineq_hi: 0 < 4 + (-1.0)*conj(z) * z
# moment_cstr: 0 < 1.0
relax_ctx = set_relaxation(problem, issparse = false, ismultiordered = false, d = 3)

########################################
# Construction du sparsity pattern
sparsity_pattern = max_cliques = 0
########################################
# Relaxation degree par clique and variables par constrainte
varsbycstr, cliquevarsbycstr, orderbyclique = 0, 0, 0
########################################
# Calcul des matrices B_i et pose du probleme
momentmatrices = compute_momentmat(problem, max_cliques, cliquevarsbycstr, orderbyclique, relax_ctx)

B_i = compute_Bibycstr(problem, momentmatrices, max_cliques, cliquevarsbycstr, orderbyclique, relax_ctx)

SDP_SOS = build_SDP_SOS(problem, max_cliques, B_i, cliquevarsbycstr, orderbyclique, relax_ctx);


m, Zi, yα_re, yα_im = make_JuMPproblem(SDP_SOS, SCSSolver(eps=1e-8, verbose=false))

println("-----> SDP_SOS problem size: ", Base.summarysize(m)/1024, " ko")
println("-----> JuMP problem size: ", Base.summarysize(m)/1024, " ko")

########################################
# Calcul d'une solution par un solveur
println("\n-----> Starting solve")
solve(m)

println("\n\nObjective value: ", getobjectivevalue(m), "\n")
for (cstrname, mmb) in B_i
    println("$cstrname \t= ", getvalue(Zi[cstrname]), "\n")
end

println("\n\n----->Lagrange multipliers : yα =")
yα = - getdual(yα_re) - im*getdual(yα_im)
print_cmat(transpose(yα))
