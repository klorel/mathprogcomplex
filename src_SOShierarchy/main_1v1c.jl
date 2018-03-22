using JuMP, SCS

ROOT = pwd()
include(joinpath(ROOT, "src_PowSysMod", "PowSysMod_body.jl"))
include(joinpath("src_SOShierarchy", "func_definitions.jl"))

# function main()

########################################
# Construction du problème type

z = Variable("z", Complex)
problemraw = Problem()
add_variable!(problemraw, z)
set_objective!(problemraw, -z-conj(z))
add_constraint!(problemraw, "ineq", abs2(z) << 1)

print(problemraw)

########################################
# Normalizing pb and setting relaxation order by constraint
problem = normalize_problem(problemraw);

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

function Λ(A::AbstractMatrix)
    n = size(A,1)
    B = zeros(2n,2n)
    B[1:n, 1:n] = real(A)
    B[1:n, (n+1):(2n)] = -imag(A)
    B[(n+1):2n, 1:n] = imag(A)
    B[(n+1):2n, (n+1):2n] = real(A)
    return B
end

function Λreal(A::AbstractMatrix)
    n, m = size(A)
    return A[1:Int(n/2), 1:Int(m/2)]
end

function Λimag(A::AbstractMatrix)
    n, m = size(A)
    return A[(Int(n/2+1)):n, 1:Int(m/2)]
end

Zi["ineq_lo"] = @variable(m, ineq_lo[1:2*2, 1:2*2], SDP)
@constraint(m, ineq_sdp1[i=1:2,j=1:2], ineq_lo[i, j] == ineq_lo[i+2, j+2])
@constraint(m, ineq_sdp2[i=1:2,j=1:2], ineq_lo[i+2, j] == -ineq_lo[j, i+2])
Zi["moment_cstr"] = @variable(m, moment_cstr[1:3*2, 1:3*2], SDP)
@constraint(m, moment_cstr1[i=1:3,j=1:2], moment_cstr[i, j] == moment_cstr[i+3, j+3])
@constraint(m, moment_cstr2[i=1:3,j=1:3], moment_cstr[i+3, j] == -moment_cstr[j, i+3])

# macro def_hermitianmatrix(m, n, matname)
#     return quote
#         local var = @variable($m, $matname[1:2*$n, 1:2*$n], SDP)
#         @constraint($m, $matname[i=1:$n, j=1:$n], $matname[i, j] == $matname[i+$n, j+$n])
#         # @constraint(m, $(matname)_sdp2[i=1:$n, j=1:$n], $matname[i+$n, j] == -$matname[i, j+$n])
#         return var
#     end
# end

# @macroexpand @def_hermitianmatrix(m, 2, toto)

# macro time(ex)
#     return quote
#         local t0 = time()
#         local val = $ex
#         local t1 = time()
#         println("elapsed time: ", t1-t0, " seconds")
#         val
#     end
# end

## Setting objective
expo = Exponent()

f0 = haskey(problem.objective, expo) ? problem.objective[expo] : 0
objectif = f0
for (i, mmbi) in B_i
    # println("-- $i -- obj = $objectif")
    if haskey(mmbi.basis, expo)
        # println("$(mmbi.basis[expo])")
        # objectif -= vecdot(Zi[i], Λ(mmbi.basis[expo]))

        Biα = mmbi.basis[expo]
        objectif -= vecdot(Λreal(Zi[i]), real(Biα)) + vecdot(Λimag(Zi[i]), imag(Biα))
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

    ## NOTE: alpha, beta and beta, alpha should provide the same constraint, to be checked.
    if expo != Exponent()
        f_α = haskey(problem.objective, expo) ? problem.objective[expo] : 0
        pss_re, pss_im = 0, 0
        for (i, mmbi) in B_i
            if haskey(mmbi.basis, expo)
                println("$expo --> $i")
                Biα = mmbi.basis[expo]
                pss_re += vecdot(Λreal(Zi[i]), real(Biα)) - vecdot(Λimag(Zi[i]), imag(Biα))
                pss_im += vecdot(Λreal(Zi[i]), imag(Biα)) + vecdot(Λimag(Zi[i]), real(Biα))
            end
        end

        # println("$expo -- $f_α -- $pss")
        # @constraint(m, f_α == pss)

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
end


print(m)


########################################
# Calcul d'une solution par un solveur
# m = make_JuMPproblem(SDP_SOS, SCSSolver())

solve(m)

# end

# main()
