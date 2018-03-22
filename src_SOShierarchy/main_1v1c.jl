workspace()

using JuMP, SCS

ROOT = pwd()
include(joinpath(ROOT, "src_PowSysMod", "PowSysMod_body.jl"))
include(joinpath("src_SOShierarchy", "func_definitions.jl"))

########################################
# Construction du problème type
z = Variable("z", Complex)
problemraw = Problem()
add_variable!(problemraw, z)
set_objective!(problemraw, (z+conj(z))/2)
add_constraint!(problemraw, "ineq", abs2(z) << 4)

print(problemraw)

########################################
# Normalizing pb and setting relaxation order by constraint
problem = normalize_problem(problemraw);

print(problem)
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

m = Model(solver = SCSSolver())
## Variables
Zi = Dict{String, Array{JuMP.Variable, 2}}()
# for (cstrname, mmb) in B_i
#     n = size(collect(values(mmb.basis))[1], 1)
#     Zi[cstrname*"_Re"] = @variable(m, [1:n, 1:n], SDP)
#     Zi[cstrname*"_Im"] = @variable(m, [1:n, 1:n], SDP)
#     println("Added two SDP var $cstrname, size $n")
# end

function Λreal(A::AbstractMatrix)
    n, m = size(A)
    return A[1:Int(n/2), 1:Int(m/2)]
end

function Λimag(A::AbstractMatrix)
    n, m = size(A)
    return A[(Int(n/2+1)):n, 1:Int(m/2)]
end

Zi["ineq_hi"] = @variable(m, ineq_lo[1:2*2, 1:2*2], SDP)
@constraint(m, ineq_sdp1[i=1:2,j=1:2], ineq_lo[i, j] == ineq_lo[i+2, j+2])
@constraint(m, ineq_sdp2[i=1:2,j=1:2], ineq_lo[i+2, j] == -ineq_lo[j, i+2])
Zi["moment_cstr"] = @variable(m, moment_cstr[1:3*2, 1:3*2], SDP)
@constraint(m, moment_cstr1[i=1:3,j=1:2], moment_cstr[i, j] == moment_cstr[i+3, j+3])
@constraint(m, moment_cstr2[i=1:3,j=1:3], moment_cstr[i+3, j] == -moment_cstr[j, i+3])

## Setting objective
expo = Exponent()

f0 = haskey(problem.objective, expo) ? problem.objective[expo] : 0
objectif = f0
for (i, mmbi) in B_i
    if haskey(mmbi.basis, expo)

        Biα = mmbi.basis[expo]
        objectif -= vecdot(Λreal(Zi[i]), real(Biα)) + vecdot(Λimag(Zi[i]), imag(Biα))
    end
end
@objective(m, Max, objectif)

## Setting constraints
d = relax_ctx.di["moment_cstr"]
vars = Set([Variable(varname, vartype) for (varname, vartype) in problem.variables])
realexpos = compute_exponents(vars, d)
conjexpos = compute_exponents(vars, d, compute_conj=true)

@constraintref myCons_re[1:3, 1:3]
@constraintref myCons_im[1:3, 1:3]

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

for re in realexpos, im in conjexpos
    expo = product(re, im)

    ## NOTE: alpha, beta and beta, alpha should provide the same constraint, to be checked.
    if expo != Exponent()
        f_α = haskey(problem.objective, expo) ? problem.objective[expo] : 0
        pss_re, pss_im = 0, 0
        for (i, mmbi) in B_i
            if haskey(mmbi.basis, expo)
                # println("$expo --> $i")
                Biα = mmbi.basis[expo]
                pss_re += vecdot(Λreal(Zi[i]), real(Biα)) - vecdot(Λimag(Zi[i]), imag(Biα))
                pss_im += vecdot(Λreal(Zi[i]), imag(Biα)) + vecdot(Λimag(Zi[i]), real(Biα))
            end
        end

        println("--> $(expo2int[re]), $(expo2int[conj(im)])")
        # println("$(expo2int[re]), $(expo2int[conj(im)]) --> $(myCons_re[expo2int[re], expo2int[conj(im)]])")
        if !(real(f_α) == 0 && pss_re == 0)
            cstr = @constraint(m, real(f_α) == pss_re)
            myCons_re[expo2int[re], expo2int[conj(im)]] = cstr
            println(myCons_re[expo2int[re], expo2int[conj(im)]])
        else
            println("-> Expo $expo provided no constraint")
        end

        if !(imag(f_α) == 0 && pss_im == 0)
            cstr = @constraint(m, imag(f_α) == pss_im)
            myCons_im[expo2int[re], expo2int[conj(im)]] = cstr
            println(myCons_im[expo2int[re], expo2int[conj(im)]])
        else
            println("-> Expo $expo provided no constraint")
        end
    end
end

for re in realexpos, im in conjexpos
    println("$(expo2int[re]), $(expo2int[conj(im)]) --> $(myCons_re[expo2int[re], expo2int[conj(im)]])")
end

print(m)


########################################
# Calcul d'une solution par un solveur
solve(m)

println("Objective value: ", getobjectivevalue(m))
println("ineq_lo = ", getvalue(ineq_lo))
println("moment_cstr = ", getvalue(moment_cstr))

println(getdual(myCons_re[1:3, 1:3]))

myCons_re
myCons_im

println(-getdual(myCons_re[1:3, 1:3]))
println(-getdual(myCons_im[1:3, 1:3]))