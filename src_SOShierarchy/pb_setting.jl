"""
    norm_pb = normalize_problem(problem)

    Transform the problem so it contains only inequalities s.a. g_i(z) ≥ 0, add the moment constraint g_0(z) = 1 ≥ 0.
"""
# NOTE: Care for integer variables will have to be taken at some point.
function normalize_problem(problem)
    println("\n=== normalize_problem(problem)")
    println("Transform the problem so it contains only ineqs s.a. g_i(z) ≥ 0, add the moment constraint g_0(z)=1≥0.")

    normpb = Problem()
    normpb.variables = copy(problem.variables)
    normpb.objective = copy(problem.objective)

    for (cstrname, cstr) in get_constraints(problem)
        if cstr.lb != -Inf-im*Inf
            add_constraint!(normpb, cstrname*"_lo", 0 << (cstr.p - cstr.lb))
        end
        if cstr.ub != Inf+im*Inf
            add_constraint!(normpb, cstrname*"_hi", 0 << (cstr.ub - cstr.p))
        end
    end

    add_constraint!(normpb, "moment_cstr", 0 << Exponent())

    exposet = Set()
    nb_expotot = 0
    degbycstr = Int64[]
    for (cstrname, cstr) in get_constraints(normpb)
        push!(degbycstr, max(cstr.p.degree.explvar, cstr.p.degree.conjvar))
        for (expo, λ) in cstr.p
            push!(exposet, expo)
            nb_expotot += 1
        end
    end
    println("-> Problem with:")
    println("-> Nb of variables:                        $(length(normpb.variables))")
    println("-> Nb of constraints:                      $(length(normpb.constraints))")
    println("-> Nb of monomials:                        $nb_expotot ($(length(exposet)) different)")
    println("-> Max total degree by constraint:         $(mean(degbycstr)) / $(std(degbycstr)) (mean/std)")
    return normpb
end


mutable struct RelaxationContext
    issparse
    ismultiordered
    di
    ki
end

"""
    relax_ctx = set_relaxation(problem; issparse = false, ismultiordered = false, d = 1)

    Build a `relax_ctx` object containing relaxation choices and problem features : order by constraint, relaxation order by constraint.
"""
function set_relaxation(problem; issparse = false, ismultiordered = false, d = 1)
    println("\n=== set_relaxation(problem, issparse = $issparse, ismultiordered = $ismultiordered, d = $d)")
    println("Build a `relax_ctx` object containing relaxation choices and problem features : polynomial degree by constraint, relaxation order by constraint...")

    ## TODO: Problem should be normalized
    nb_densecstrs = 0
    maxdeg_densecstr = Float32[]
    ki = Dict{String, Int}()
    di = Dict{String, Int}()

    for (cstrname, cstr) in problem.constraints
        ki[cstrname] = max(cstr.p.degree.explvar, cstr.p.degree.conjvar)
        di[cstrname] = d

        if di[cstrname] > ki[cstrname]
            nb_densecstrs += 1
            push!(maxdeg_densecstr, ki[cstrname])
        end
    end

    println("-> Number of constraints s.t. di > ki:     $nb_densecstrs / $(length(problem.constraints))")
    println("-> Max total degree on such constraints:   $(mean(maxdeg_densecstr)) / $(std(maxdeg_densecstr)) (mean/std)")
    println("All variables appearing in such constraints will be linked in the sparsity pattern, which will largely densify it.")
    return RelaxationContext(issparse, ismultiordered, di, ki)
end
