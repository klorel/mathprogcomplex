"""
    relax_ctx = set_relaxation(pb::Problem; ismultiordered=false, issparse=false, leveragesymmetries=true, hierarchykind=:Complex, renamevars=false, di=SortedDict{String, Int}(), d=-1)

    Build a `relax_ctx` object containing relaxation choices and problem features : order by constraint, relaxation order by constraint...
"""
function set_relaxation(pb::Problem; ismultiordered=false, 
                                     issparse=false, 
                                     leveragesymmetries=true,
                                     hierarchykind=:Complex,
                                     renamevars=false,
                                     di=SortedDict{String, Int}(),
                                     d=-1)
    println("\n=== set_relaxation(pb; ismultiordered=$ismultiordered, issparse=$issparse, leveragesymmetries=$leveragesymmetries, hierarchykind=$hierarchykind, renamevars=$renamevars, di=SortedDict of length $(length(di)), d=$d)")

    # Compute each constraint degree
    ki = SortedDict{String, Int}()
    for (cstrname, cstr) in pb.constraints
        ki[cstrname] = max(cstr.p.degree.explvar, cstr.p.degree.conjvar)
    end

    # Check that either d or di was provided as input
    ((di == SortedDict{String, Int}()) ⊻ (d==-1)) || error("RelaxationContext(): Either di or d should be provided as input, not both.")

    if d!=-1
        for (cstr, ki_) in ki
            (ki_ <= d) || warn("RelaxationContext(): Provided d ($d) is lower than constraint $cstr order ($ki_). \nUsing value $ki_, hierarchy may be multiordered.")
            di[cstr] = max(ki_, d)
        end
    else
        (SortedSet(keys(di)) == SortedSet(keys(ki))) || error("RelaxationContext(): Provided di doesn't match the set of constraint names.")
        for (cstr, ki_) in ki
            di_ = di[cstr]
            (ki_ <= di_) || warn("RelaxationContext(): Provided di ($di_) is lower than constraint $cstr order ($ki_). \nUsing value $ki_.")
        end
    end

    # Check that all variables have a type fitting the hierarchy kind
    for (varname, vartype) in pb.variables
        (vartype<:Int) && error("set_relaxation() : variable $varname,$vartype is integer, unfit for SOS relaxation.\nConsider relaxing it and adding a complementarity constraint.")
        (hierarchykind==:Complex) && !(vartype<:Complex) && error("set_relaxation() : variable $varname,$vartype should be complex for complex hierarchy.")
        (hierarchykind==:Real) && !(vartype<:Real) && error("set_relaxation() : variable $varname,$vartype should be real for real hierarchy.")
    end

    # log intel
    nb_densecstrs = 0
    maxdeg_densecstr = Float32[]
    for (cstr, ki_) in ki
        if di[cstr] > ki_
            nb_densecstrs += 1
            push!(maxdeg_densecstr, ki_)
        end
    end

    @printf("-> Number of constraints s.t. di > ki:     %i / %i\n", nb_densecstrs, length(pb.constraints))
    @printf("-> Max total degree on such constraints:   %.1f / %.1f (mean/std)\n", mean(maxdeg_densecstr), std(maxdeg_densecstr))
    @printf("All variables appearing in such constraints will be linked in the sparsity pattern, which will largely densify it.\n")
    return RelaxationContext(ismultiordered, issparse, leveragesymmetries, hierarchykind, renamevars, di, ki)
end


"""
    norm_pb = normalize_problem(problem)

    Transform the problem so it contains only ineqs or eqs s.a. g_i(z) ≥ 0 or g_j(z) = 0, add the moment constraint g_0(z) = 1 ≥ 0.
"""
# NOTE: Care for integer variables will have to be taken at some point.
function normalize_problem(problem)
    println("\n=== normalize_problem(problem)")
    println("Transform the problem so it contains only ineqs or eqs s.a. g_i(z) ≥ 0 or g_j(z) = 0, add the moment constraint g_0(z)=1≥0.")

    normpb = Problem()
    normpb.variables = deepcopy(problem.variables)
    normpb.objective = deepcopy(problem.objective)

    for (cstrname, cstr) in get_constraints(problem)
        if cstr.lb != cstr.ub
            if cstr.lb != -Inf-im*Inf
                add_constraint!(normpb, cstrname*"_lo", 0 << (cstr.p - cstr.lb))
            end
            if cstr.ub != Inf+im*Inf
                add_constraint!(normpb, cstrname*"_hi", 0 << (cstr.ub - cstr.p))
            end
        else
            add_constraint!(normpb, cstrname, (cstr.p - cstr.ub) == 0)
        end
    end

    add_constraint!(normpb, "moment_cstr", 0 << Exponent())

    exposet = SortedSet()
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



function print(io::IO, relctx::RelaxationContext)
    print(io, "RelaxationContext:\n")
    print(io, "ismultiordered         : $(relctx.ismultiordered)\n")
    print(io, "issparse               : $(relctx.issparse)\n")
    print(io, "leveragesymmetries     : $(relctx.leveragesymmetries)\n")
    print(io, "hierarchykind          : $(relctx.hierarchykind)\n")
    print(io, "renamevars             : $(relctx.renamevars)\n")
    print(io, "di                     : $(relctx.di)\n")
    print(io, "ki                     : $(relctx.ki)")
end