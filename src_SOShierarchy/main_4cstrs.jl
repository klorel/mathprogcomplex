using JuMP, SCS

ROOT = pwd()
include(joinpath(ROOT, "src_PowSysMod", "PowSysMod_body.jl"))
include(joinpath("func_definitions.jl"))

# function main()

########################################
# Construction du problème type
# OPFpbs = load_OPFproblems(MatpowerInput, joinpath("data_Matpower", "matpower", "WB2.m"))
# problemraw = build_globalpb!(OPFpbs)
#
# println("WB2 problem built")

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
problem = normalize_problem(problemraw)

print(problem)
relax_ctx = set_relaxation(problem, issparse = false, ismultiordered = false, d = 2)

########################################
# Construction du sparsity pattern
sparsity_pattern = compute_sparsitypattern(problem, relax_ctx)

# Extension chordale et détection des cliques maximales
compute_chordalextension!(sparsity_pattern)
max_cliques = compute_maxcliques(sparsity_pattern)

########################################
# Relaxation degree par clique and variables par constrainte
varsbycstr = compute_varsbycstr(problem)
cliquevarsbycstr = compute_varsbycstr(sparsity_pattern, max_cliques, varsbycstr)
orderbyclique = compute_cliqueorders(sparsity_pattern, varsbycstr, max_cliques, relax_ctx)

########################################
# Calcul des matrices B_i et pose du probleme
momentmatrices = compute_momentmat(problem, max_cliques, cliquevarsbycstr, orderbyclique, relax_ctx)

for (name, mat) in momentmatrices
    println("$name   →\n$mat")
end

B_i = compute_Bibycstr(problem, momentmatrices, max_cliques, cliquevarsbycstr, orderbyclique, relax_ctx)

for (expo, mat) in B_i["ineq_lo"].basis
    println("$expo  ->\n$mat")
end

# SDP_SOS = build_SDP_SOS(problem, max_cliques, B_i, cliquevarsbycstr, orderbyclique, relax_ctx)

m = Model(solver = SCSSolver())

## Variables
Zi = Dict{String, Array{JuMP.Variable, 2}}()
for (cstrname, mmb) in B_i
    n = size(collect(values(mmb.basis))[1], 1)
    Zi[cstrname*"_Re"] = @variable(m, [1:n, 1:n], SDP)
    Zi[cstrname*"_Im"] = @variable(m, [1:n, 1:n], SDP)
    println("Added two SDP var $cstrname, size $n")
end

## Setting objective
expo = Exponent()

f0 = haskey(problem.objective, expo) ? problem.objective[expo] : 0
objectif = f0
for (i, mmbi) in B_i
    # println("-- $i -- obj = $objectif")
    if haskey(mmbi.basis, expo)
        # println("$(mmbi.basis[expo])")
        objectif += trace(Zi[i*"_Re"] * real(mmbi.basis[expo]) - Zi[i*"_Im"] * imag(mmbi.basis[expo]))
        # println("⫐ $(trace(Zi[i*"_Re"] * real(mmbi.basis[expo]) - Zi[i*"_Im"] * imag(mmbi.basis[expo])))")
        # println("⫐ $(trace(Zi[i*"_Im"] * real(mmbi.basis[expo]) + Zi[i*"_Re"] * imag(mmbi.basis[expo])))")
    end
end
# println("final obj: $objectif")
@objective(m, Max, objectif)

## Setting constraints
d = relax_ctx.di["moment_cstr"]
vars = Set([Variable(varname, vartype) for (varname, vartype) in problem.variables])
realexpos = compute_exponents(vars, d)
conjexpos = compute_exponents(vars, d, compute_conj=true)

for re in realexpos, im in conjexpos
    expo = product(re, im)

    f_α = haskey(problem.objective, expo) ? problem.objective[expo] : 0
    pss_re, pss_im = 0, 0
    for (i, mmbi) in B_i
        if haskey(mmbi.basis, expo)
            Biα = mmbi.basis[expo]
            pss_re += trace(Zi[i*"_Re"] * real(Biα) - Zi[i*"_Im"] * imag(Biα))
            pss_im += trace(Zi[i*"_Re"] * imag(Biα) + Zi[i*"_Im"] * real(Biα))
        end
    end

    if !(real(f_α) == 0 && pss_re == 0)
        @constraint(m, real(f_α) == pss_re)
    else
        println("-> Expo $expo provided no constraint")
    end

    if !(imag(f_α) == 0 && pss_im == 0)
        @constraint(m, imag(f_α) == pss_im)
    else
        println("-> Expo $expo provided no constraint")
    end
end


print(m)


########################################
# Calcul d'une solution par un solveur
m = make_JuMPproblem(SDP_SOS, SCSSolver())

solve(m)

# end

# main()
