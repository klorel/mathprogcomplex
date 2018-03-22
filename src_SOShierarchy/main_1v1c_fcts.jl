workspace()

using JuMP, SCS

ROOT = pwd()
include(joinpath(ROOT, "src_PowSysMod", "PowSysMod_body.jl"))
include(joinpath(ROOT, "src_SOShierarchy", "func_definitions.jl"))

########################################
# Construction du problème type
z = Variable("z", Complex)
problemraw = Problem()
add_variable!(problemraw, z)
set_objective!(problemraw, imag(z))
add_constraint!(problemraw, "ineq", abs2(z) << 4)

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
relax_ctx = set_relaxation(problem, issparse = false, ismultiordered = false, d = 2)

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

m = Model(solver = SCSSolver())

## Variables
Zi = Dict{String, Array{JuMP.Variable, 2}}()
for (cstrname, n) in SDP_SOS.variables
    var = @variable(m, [1:2*n,1:2*n], SDP, basename=cstrname)
    Zi[cstrname] = var
    cstr1 = @constraint(m, [i=1:n,j=1:n], var[i, j] == var[i+n, j+n])
    for i in 1:n, j in 1:n
        println("$i, $j  --> $(cstr1[i, j])")
    end
    cstr2 = @constraint(m, [i=1:n,j=1:i], var[n+i, j] == - var[n+j, i])
    for i in 1:n, j in 1:i
        println("$i, $j --> $(cstr2[i, j])")
    end
end

function formconstraint(SDPform, Zi)
    l_re = real(SDPform.scal)
    l_im = imag(SDPform.scal)
    for (i, Biα) in SDPform.mats
        l_re -= vecdot(Λreal(Zi[i]), real(Biα)) - vecdot(Λimag(Zi[i]), imag(Biα))
        l_im -= vecdot(Λreal(Zi[i]), imag(Biα)) + vecdot(Λimag(Zi[i]), real(Biα))
    end
    return l_re, l_im
end

## Setting objective
expo = Exponent()

objectif_re, _ = formconstraint(SDP_SOS.objective, Zi)
@objective(m, Max, objectif_re)

## Setting constraints
d = relax_ctx.di["moment_cstr"]
vars = Set([Variable(varname, vartype) for (varname, vartype) in problem.variables])
realexpos = compute_exponents(vars, d)
conjexpos = compute_exponents(vars, d, compute_conj=true)

#Constraint storage
@constraintref myCons_re[1:length(realexpos), 1:length(realexpos)]
@constraintref myCons_im[1:length(conjexpos), 1:length(conjexpos)]

expo2int = Dict{Exponent, Int}()
int2expo = Dict{Int, Exponent}()
i = 1
for expo in sort(collect(realexpos))
    expo2int[expo] = i
    int2expo[i] = expo
    i += 1
end

myCons_re[1,1] = @constraint(m, 1==1)
myCons_im[1,1] = @constraint(m, 1==1)

#Building constraints
for re in realexpos, im in conjexpos
    expo = product(re, im)

    ## NOTE: (alpha, beta) and (beta, alpha) should provide the same constraint, to be checked.
    if expo != Exponent()
        println("--> $(expo2int[re]), $(expo2int[conj(im)])")
        l_re, l_im = formconstraint(SDP_SOS.constraints[expo], Zi)

        if l_re != 0
            myCons_re[expo2int[re], expo2int[conj(im)]] = @constraint(m, l_re == 0)
            println(myCons_re[expo2int[re], expo2int[conj(im)]])
        else
            println("-> Expo $expo provided no constraint")
        end

        if l_im != 0
            myCons_im[expo2int[re], expo2int[conj(im)]] = @constraint(m, l_im == 0)
            println(myCons_im[expo2int[re], expo2int[conj(im)]])
        else
            println("-> Expo $expo provided no constraint")
        end
    end
end

print(m)


########################################
# Calcul d'une solution par un solveur
println("-----> Starting solve")
solve(m)

println("Objective value: ", getobjectivevalue(m))
for (cstrname, mmb) in B_i
    println("$cstrname = ", getvalue(Zi[cstrname]))
end

println(getdual(myCons_re[1:3, 1:3]))
println(-getdual(myCons_im[1:3, 1:3]))