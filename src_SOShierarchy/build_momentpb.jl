"""
    mm = MomentMatrix(vars::SortedSet{Variable}, d, symmetries)

    Build the moment matrix corresponding to the moment of degree up to `d` of the `vars` polynomial algebra.
    Only monomials featuring all `symmetries` appear in the moment matrix.
"""
function MomentMatrix(relax_ctx, vars::SortedSet{Variable}, d::Int, symmetries::SortedSet{DataType}, matrixkind::Symbol)
    mm = SortedDict{Tuple{Exponent, Exponent}, AbstractPolynomial}()
    realexpos = compute_exponents(vars, d)
    conjexpos = compute_exponents(vars, d, compute_conj=true)
    for cexp in conjexpos
        for rexp in realexpos
            expo = cexp*rexp
            issym = true
            for sym in symmetries
                issym = issym && has_symmetry(relax_ctx, expo, sym)
            end
            if !issym
                continue
            end
            if cexp ≥ rexp
                mm[(cexp, rexp)] = cexp*rexp
            end
        end
    end
    return MomentMatrix(mm, SortedSet(vars), d, matrixkind)
end

function copy(mm::MomentMatrix)
    return MomentMatrix(copy(mm.mm), mm.vars, mm.order, mm.matrixkind)
end

function print(io::IO, mm::MomentMatrix)
    for (key, val) in mm.mm
        println(io, "($(key[1]), $(key[2])) ⟶  $val")
    end
    print(io, " $(mm.matrixkind)")
end


##########################
# Moment matrix algebra
function product!(p::T, mm::MomentMatrix) where T<:AbstractPolynomial
    for (key, val) in mm.mm
        mm.mm[key] = p * mm.mm[key]
    end
end
function product!(λ::Number, mm::MomentMatrix)
    for (key, val) in mm.mm
        mm.mm[key] = λ*val
    end
end

function *(p::T, mm::MomentMatrix) where T<:AbstractPolynomial
    mm_copy = copy(mm)
    product!(p, mm_copy)
    return mm_copy
end
*(mm::MomentMatrix, p::T) where T<:AbstractPolynomial = p*mm

function *(λ::Number, mm::MomentMatrix)
    mm_copy = copy(mm)
    product!(λ, mm_copy)
    return mm_copy
end
*(mm::MomentMatrix, λ::Number) = λ*mm

function evaluate(mm::MomentMatrix, pt::Point)
    mm_eval = SortedDict{Tuple{Exponent, Exponent}, AbstractPolynomial}()
    for (key, p) in mm.mm
        res = evaluate(p, pt)
        if res == Polynomial()
            delete!(mm_eval, key)
        else
            mm_eval[key] = res
        end
    end
    return MomentMatrix(mm_eval, setdiff(mm.vars, SortedSet(keys(pt))), mm.order, mm.matrixkind)
end


"""
    momentrelaxation = MomentRelaxationPb(relax_ctx, problem, moment_param::SortedDict{String, Tuple{SortedSet{String}, Int}}, max_cliques::SortedDict{String, SortedSet{Variable}})

    Compute the `momentrelaxation` of `problem` corresponding to the clique decomposition `max_cliques` and parameters `moment_param`.
"""
function MomentRelaxationPb(relax_ctx, problem, momentmat_param::SortedDict{String, Int},
                                                localizingmat_param::SortedDict{String, Tuple{SortedSet{String}, Int}},
                                                max_cliques::SortedDict{String, SortedSet{Variable}};
                                                verbose=false)
    println("\n=== MomentRelaxationPb(relax_ctx, problem, moment_param::SortedDict{String, Tuple{SortedSet{String}, Int}}, max_cliques::SortedDict{String, SortedSet{Variable}})")
    println("Compute the moment and localizing matrices associated with the problem constraints and clique decomposition and return a MomentRelaxationPb object.")

    verbose && println("Computing overlapping variables...")
    ## Locate clique overlapping variables
    vars_overlap = SortedDict{Variable, SortedSet{String}}()

    # Collect all variables
    variables = SortedSet{Variable}()
    for (clique, clvars) in max_cliques
        union!(variables, clvars)
    end

    # Collect cliques by variable
    for var in variables
        for (clique, clique_vars) in max_cliques
            if var in clique_vars
                haskey(vars_overlap, var) || (vars_overlap[var] = SortedSet{String}())

                insert!(vars_overlap[var], clique)
            end
        end
    end

    # Delete variables appearing in one clique only
    for (var, cliques) in vars_overlap
        length(cliques) > 1 || delete!(vars_overlap, var)
    end

    verbose && for (var, cliques) in vars_overlap println("$var, $(collect(cliques))") end


    ## Build moment matrix
    momentmatrices = SortedDict{Tuple{String, String}, MomentMatrix}()
    for (cliquename, vars) in max_cliques

        # Compute set of variables in cliquename
        vars_inclique = SortedSet{Variable}()
        for var in vars
            var_inclique = var
            if haskey(vars_overlap, var)
                @assert haskey(vars_overlap[var], cliquename)
                var_inclique = get_varinclique(var, cliquename)
            end
            insert!(vars_inclique, var_inclique)
        end
        verbose && println("Building moment matrix: clique $cliquename, $(length(vars)) variables")
        # warn("Moment matrix vars : $vars_inclique")
        dcl = momentmat_param[cliquename]
        momentmatrices[(get_momentcstrname(), cliquename)] = MomentMatrix(relax_ctx, vars_inclique, dcl, relax_ctx.symmetries, :SDP)
    end

    ## Build localizing matrices
    for (cstrname, cstr) in problem.constraints

        cstrtype = get_cstrtype(cstr)
        if cstrtype == :ineqdouble
            cstrname_lo, cstrname_up = get_cstrname(cstrname, cstrtype)

            # Deal with lower inequality
            clique_keys, order = localizingmat_param[cstrname_lo]
            vars, cliquename = collect_cliquesvars(clique_keys, max_cliques, vars_overlap)

            verbose && print("Building localizing matrix: constraint $cstrname_lo - $cliquename, $(length(vars)) variables : building moment matrix... ")
            mmt = MomentMatrix(relax_ctx, vars, order, relax_ctx.symmetries, relax_ctx.cstrtypes[cstrname_lo])
            verbose && println("computing product with polynomial...")
            momentmatrices[(cstrname_lo, cliquename)] = mmt * (cstr.p - cstr.lb)

            # Deal with upper inequality, no recomputing of variables or moment matrix if possible
            clique_keys_up, order_up = localizingmat_param[cstrname_up]
            if collect(clique_keys) != collect(clique_keys_up)
                warn("clique keys different from lower and upper side of double constraint")
                vars, cliquename = collect_cliquesvars(clique_keys_up, max_cliques, vars_overlap)

                verbose && print("Building localizing matrix: constraint $cstrname_up - $cliquename, $(length(vars)) variables : building moment matrix... ")
                mmt = MomentMatrix(relax_ctx, vars, order_up, relax_ctx.symmetries, relax_ctx.cstrtypes[cstrname_up])
            elseif order_up != order
                warn("order different from lower and upper side of double constraint")

                verbose && print("Building localizing matrix: constraint $cstrname_up - $cliquename, $(length(vars)) variables : building moment matrix... ")
                mmt = MomentMatrix(relax_ctx, vars, order_up, relax_ctx.symmetries, relax_ctx.cstrtypes[cstrname_up])
            end

            verbose && collect(clique_keys) == collect(clique_keys_up) && order_up == order && print("Building localizing matrix: constraint $cstrname_up - $cliquename, $(length(vars)) variables : ")
            verbose && println("computing product with polynomial...")
            momentmatrices[(cstrname_up, cliquename)] = mmt * (cstr.ub - cstr.p)

        else
            # either cstrtype == :ineqlo, :ineqhi, :eq
            clique_keys, order = localizingmat_param[get_cstrname(cstrname, cstrtype)]
            vars, cliquename = collect_cliquesvars(clique_keys, max_cliques, vars_overlap)

            verbose && print("Building localizing matrix: constraint $(get_cstrname(cstrname, cstrtype)) - $cliquename, $(length(vars)) variables : building moment matrix... ")
            mmt = MomentMatrix(relax_ctx, vars, order, relax_ctx.symmetries, relax_ctx.cstrtypes[get_cstrname(cstrname, cstrtype)])
            verbose && println("computing product with polynomial...")
            warn("Variables are $vars")
            momentmatrices[(get_cstrname(cstrname, cstrtype), cliquename)] = mmt * get_normalizedpoly(cstr, cstrtype)
        end
    end

    ## Treat objective polynomial
    mmtobj = Polynomial()
    for (expo, coeff) in problem.objective
        mmtexpo = Exponent()
        for (var, degree) in expo
            varmmt = var
            if haskey(vars_overlap, var)
                varmmt = get_varinclique(var, first(collect(vars_overlap[var]))) # NOTE: proper choice to be made here ?
            end
            add_expod!(mmtexpo, Exponent(varmmt=>degree))
        end
        add!(mmtobj, mmtexpo * coeff)
    end

    return MomentRelaxationPb(mmtobj, momentmatrices, vars_overlap)
end


function print(io::IO, momentrelax::MomentRelaxationPb)
    println(io, "Moment Relaxation Problem:")
    println(io, "▶ Objective: ", momentrelax.objective)
    println(io, "▶ Constraints:")
    for ((cstrname, blocname), mmtmat) in momentrelax.constraints
        println(io, " → $cstrname, $blocname")
        println(io, mmtmat)
    end
    println(io, "▶ Variables clique overlap:")
    if length(momentrelax.vars_overlap) > 0
        for (var, cliquenames) in momentrelax.vars_overlap
            print(io, " → $var : ")
            for clique in cliquenames print(io, "$clique, ") end
            @printf(io, "\b\b \n")
        end
    else
        println(io, "  None")
    end
end