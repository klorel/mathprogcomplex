"""
    relax_ctx = set_relaxation(pb::Problem; ismultiordered=false, issparse=false, symmetries=Set(), hierarchykind=:Complex, renamevars=false, di=Dict{String, Int}(), d=-1)

    Build a `relax_ctx` object containing relaxation choices and problem features : order by constraint, relaxation order by constraint...
"""
function set_relaxation(pb::Problem; ismultiordered::Bool=false,
                                     issparse::Bool=false,
                                     symmetries::Array{Type, 1}=Type[],
                                     hierarchykind::Symbol=:Complex,
                                     renamevars::Bool=false,
                                     di::Dict{String, Int}=Dict{String, Int}(),
                                     d::Int=-1)
    println("\n=== set_relaxation(pb; ismultiordered=$ismultiordered, issparse=$issparse, symmetries=$symmetries, hierarchykind=$hierarchykind, renamevars=$renamevars, di=Dict of length $(length(di)), d=$d)")

    # Check that all variables have a type fitting the hierarchy kind
    for (varname, vartype) in pb.variables
        (vartype<:Int) && error("set_relaxation() : variable $varname,$vartype is integer, unfit for SOS relaxation.\nConsider relaxing it and adding a complementarity constraint.")
        (hierarchykind==:Complex) && !(vartype<:Complex) && error("set_relaxation() : variable $varname,$vartype should be complex for complex hierarchy.")
        (hierarchykind==:Real) && !(vartype<:Real) && error("set_relaxation() : variable $varname,$vartype should be real for real hierarchy.")
    end

    # Compute each constraint degree
    ki = Dict{String, Int}()
    ki[get_momentcstrname()] = 0
    for (cstrname, cstr) in pb.constraints
        cstrtype = get_cstrtype(cstr)
        if cstrtype == :ineqdouble
            cstrname_lo, cstrname_up = get_cstrname(cstrname, cstrtype)
            ki[cstrname_lo] = max(cstr.p.degree.explvar, cstr.p.degree.conjvar)
            ki[cstrname_up] = max(cstr.p.degree.explvar, cstr.p.degree.conjvar)
        else
            ki[get_cstrname(cstrname, cstrtype)] = max(cstr.p.degree.explvar, cstr.p.degree.conjvar)
        end
    end

    # Store each SDP multiplier type
    cstrtypes = Dict{String, Symbol}()
    for (cstrname, cstr) in pb.constraints
        cstrtype = get_cstrtype(cstr)
        if cstrtype == :ineqdouble
            cstrname_lo, cstrname_up = get_cstrname(cstrname, cstrtype)
            cstrtypes[cstrname_lo] = (hierarchykind==:Complex ? :SDPC : :SDP)
            cstrtypes[cstrname_up] = (hierarchykind==:Complex ? :SDPC : :SDP)
        elseif cstrtype == :eq
            cstrtypes[get_cstrname(cstrname, cstrtype)] = (hierarchykind==:Complex ? :SymC : :Sym)
        else
            cstrtypes[get_cstrname(cstrname, cstrtype)] = (hierarchykind==:Complex ? :SDPC : :SDP)
        end
    end
    cstrtypes[get_momentcstrname()] = (hierarchykind==:Complex ? :SDPC : :SDP)

    # Relaxation order management
    di_relax = Dict{String, Int}()
    !((di == Dict{String, Int}()) && (d==-1)) || error("RelaxationContext(): Either di or d should be provided as input.")

    for (cstrname, cstr) in pb.constraints
        cur_order = haskey(di, cstrname) ? di[cstrname] : d

        # Check provided di is suitable wrt constraint degree, add
        cstrtype = get_cstrtype(cstr)
        if cstrtype == :ineqdouble
            cstrname_lo, cstrname_up = get_cstrname(cstrname, cstrtype)
            cur_ki = ki[cstrname_lo]
            (0 ≤ cur_order-ceil(cur_ki/2)) || warn("RelaxationContext(): Provided order ($cur_order) is lower than constraint $cstrname order ($cur_ki). \nUsing value ceil($cur_ki/2).")
            # (cur_ki <= cur_order) || warn("RelaxationContext(): Provided order ($cur_order) is lower than constraint $cstrname order ($cur_ki). \nUsing value $cur_ki, hierarchy may be multiordered.")
            di_relax[cstrname_lo] = cur_order
            di_relax[cstrname_up] = cur_order
            # di_relax[cstrname_lo] = max(cur_order, ceil(cur_ki/2))
            # di_relax[cstrname_up] = max(cur_order, ceil(cur_ki/2))
        else # :eq, :ineqlo, :ineqhi
            cur_ki = ki[get_cstrname(cstrname, cstrtype)]
            (0 ≤ cur_order-ceil(cur_ki/2)) || warn("RelaxationContext(): Provided order ($cur_order) is lower than constraint $cstrname order ($cur_ki). \nUsing value ceil($cur_ki/2).")
            # (cur_ki <= cur_order) || warn("RelaxationContext(): Provided order ($cur_order) is lower than constraint $cstrname order ($cur_ki). \nUsing value $cur_ki, hierarchy may be multiordered.")
            di_relax[get_cstrname(cstrname, cstrtype)] = cur_order
            # di_relax[get_cstrname(cstrname, cstrtype)] = max(cur_order, ceil(cur_ki/2))
        end
    end

    # Moment constraint relaxation order
    if haskey(di, get_momentcstrname())
        di_relax[get_momentcstrname()] = di[get_momentcstrname()]
    elseif d!=-1
        di_relax[get_momentcstrname()] = d
    else
        di_relax[get_momentcstrname()] = maximum(values(di_relax))
    end
    # Objective polynomial must be representable by moment matrix
    obj_degree = max(pb.objective.degree.explvar, pb.objective.degree.conjvar)
    # if obj_degree > di_relax[get_momentcstrname()]
    #     warn("RelaxationContext(): Moment matrix order $(di_relax[get_momentcstrname()]) is lower than objective degree ($obj_degree). \nUsing value $obj_degree, hierarchy may be multiordered.")
    #     di_relax[get_momentcstrname()] = ceil(obj_degree/2)
    # end

    relax_ctx = RelaxationContext(ismultiordered, issparse, Set{DataType}(), hierarchykind, renamevars, di_relax, ki, cstrtypes)

    # Check whether the problem has the suggested symmetries
    pbsymmetries = Set{DataType}()
    isa(symmetries, Array) || error("set_relaxation(): symmetries should be an Array of types.")
    for symtype in symmetries
        if has_symmetry(relax_ctx, pb, symtype)
            push!(pbsymmetries, symtype)
        end
    end
    relax_ctx.symmetries = pbsymmetries

    # log intel
    nb_densecstrs = 0
    maxdeg_densecstr = Float32[]
    for (cstr, ki_) in ki
        if di_relax[cstr] > ki_
            nb_densecstrs += 1
            push!(maxdeg_densecstr, ki_)
        end
    end


    exposet = Set()
    nb_expotot = 0
    degbycstr = Int64[]
    for (cstrname, cstr) in get_constraints(pb)
        push!(degbycstr, max(cstr.p.degree.explvar, cstr.p.degree.conjvar))
        for (expo, λ) in cstr.p
            push!(exposet, expo)
            nb_expotot += 1
        end
    end
    println("=> Problem with:")
    println("-> Nb of variables:                        $(length(pb.variables))")
    println("-> Nb of constraints:                      $(length(pb.constraints))")
    println("-> Nb of monomials:                        $nb_expotot ($(length(exposet)) different)")
    @printf("-> Max total degree by constraint:         %.1f / %.1f (mean/std)\n", mean(degbycstr), std(degbycstr))

    println("=> Relaxation characteristics:")
    @printf("-> Number of constraints s.t. di > ki:     %i / %i\n", nb_densecstrs, length(pb.constraints))
    @printf("-> Max total degree on such constraints:   %.1f / %.1f (mean/std)\n", mean(maxdeg_densecstr), std(maxdeg_densecstr))
    @printf("All variables appearing in such constraints will be linked in the sparsity pattern, which will largely densify it.\n")
    return relax_ctx
end



function print(io::IO, relctx::RelaxationContext)
    print(io, "RelaxationContext:\n")
    print(io, "ismultiordered         : $(relctx.ismultiordered)\n")
    print(io, "issparse               : $(relctx.issparse)\n")
    print(io, "symmetries             : $(relctx.symmetries)\n")
    print(io, "hierarchykind          : $(relctx.hierarchykind)\n")
    print(io, "renamevars             : $(relctx.renamevars)\n")
    for (cstrname, di) in relctx.di
    print(io, "di                     : $cstrname  \t=> $di\n")
    end
    for (cstrname, ki) in relctx.ki
    print(io, "ki                     : $cstrname  \t=> $ki\n")
    end
    for (cstrname, cstrtype) in relctx.cstrtypes
    print(io, "bar var types          : $cstrname  \t=> $(string(cstrtype))\n")
    end
end