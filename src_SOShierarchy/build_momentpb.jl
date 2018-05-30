"""
    mm = MomentMatrix(vars::SortedSet{Variable}, d, symmetries)

    Build the moment matrix corresponding to the moment of degree up to `d` of the `vars` polynomial algebra.
    Only monomials featuring all `symmetries` appear in the moment matrix.
"""
function MomentMatrix(relax_ctx, vars::SortedSet{Variable}, d::Int, symmetries::SortedSet{DataType},
                                                                    matrixkind::Symbol;
                                                                    default_clique::String="",
                                                                    var_to_cliques::SortedDict{Variable, SortedSet{String}}=SortedDict{Variable, SortedSet{String}}())
    mm = SortedDict{Tuple{Exponent, Exponent}, SortedDict{Moment, Number}}()

    ## Computing exponents for available variables
    realexpos = compute_exponents(vars, d)
    conjexpos = compute_exponents(vars, d, compute_conj=true)
    for cexp in conjexpos
        for rexp in realexpos
            expo = product(cexp, rexp)

            ## Checking if current exponent has required symmetries
            issym = true
            for sym in symmetries
                issym = issym && has_symmetry(relax_ctx, expo, sym)
            end

            ## Storing only lower triangular matrix
            if issym && cexp ≥ rexp
                @assert default_clique!="" || var_to_cliques!=SortedDict{Variable, SortedSet{String}}()

                # Get exponent clique
                expo_clique = default_clique
                if var_to_cliques!=SortedDict{Variable, SortedSet{String}}()
                    expo_clique = get_exponentclique(expo, var_to_cliques)
                end

                mm[(cexp, rexp)] = SortedDict{Moment, Number}(Moment(expo, expo_clique)=>1)
            end
        end
    end
    return MomentMatrix(mm, SortedSet(vars), d, matrixkind)
end

function copy(mm::MomentMatrix)
    return MomentMatrix(copy(mm.mm), mm.vars, mm.order, mm.matrixkind)
end

function print(io::IO, mm::MomentMatrix)
    keylen = maximum(x->length("($(x[1]), $(x[2])) ⟶  "), keys(mm.mm))

    maxorder = 0

    for (key, momentpoly) in mm.mm
        print_string(io, "($(key[1]), $(key[2])) ⟶  ", keylen)

        (moment_first, val_first) = first(momentpoly)
        println(io, "$(moment_first.clique) -- $(moment_first.conj_part) × $(moment_first.expl_part) × $val_first")

        for (moment, val) in momentpoly
            maxorder = max(maxorder, moment.expl_part.degree.explvar)

            (moment_first, val_first) == (moment, val) && continue
            println(io, " "^(keylen+1), "$(moment.clique) -- $(moment.conj_part) × $(moment.expl_part) × $val")
        end
    end
    info(io, "maxorder is $maxorder")
    print(io, " $(mm.matrixkind)")
end

"""
    cliquename = get_exponentclique(expo, var_to_cliques)

    Determine which clique expo fits in, that is which cliques contain all variables of expo.
    Error if no such clique are found.
"""
function get_exponentclique(expo, var_to_cliques)
    cliques = SortedSet{String}()

    ## If expo is one, return default clique
    expo == Exponent() && return "clique_un"

    union!(cliques, var_to_cliques[first(expo)[1]])
    for (var, deg) in expo
        cliques = intersect(cliques, var_to_cliques[var])
    end

    length(cliques) == 0 && error("get_exponentclique(): $expo is split amongst several cliques.\nMaximal cliques provided are not suitable for this relaxation.")

    clique = first(cliques)
    # length(cliques) > 1 && warn("get_exponentclique(): $expo appears in $(length(cliques)) cliques : $cliques.\nChoosing first one $clique") ## TODO: better logging system
    return clique
end

##########################
## Moment matrix algebra
##########################

## AbstractPolynomial types
function product!(mm::MomentMatrix, p::T, var_to_cliques) where T<:Union{AbstractPolynomial, Number}
    for (key, momentpoly) in mm.mm
        mm.mm[key] = product(momentpoly, p, var_to_cliques)
    end
end

function product(momentpoly::SortedDict{Moment, Number}, p::T, var_to_cliques) where T<:Union{AbstractPolynomial, Number}
    return product(momentpoly, convert(Polynomial, p), var_to_cliques)
end

function product(momentpoly::SortedDict{Moment, Number}, p::Polynomial, var_to_cliques)
    resmpoly = SortedDict{Moment, Number}()

    for (expo, val1) in p
        for (moment, val2) in momentpoly
            resmoment = product(moment, expo, var_to_cliques)

            haskey(resmpoly, resmoment) || (resmpoly[resmoment] = 0)
            resmpoly[resmoment] += val1*val2
            (resmpoly[resmoment] == 0) && delete!(resmpoly, resmoment)
        end
    end

    return resmpoly
end

function product(mmt, expo, var_to_cliques)
    respoly = mmt.conj_part * mmt.expl_part * expo
    resexpo = first(respoly)[1]

    clique = get_exponentclique(resexpo, var_to_cliques)
    return Moment(resexpo, clique)
end

# function evaluate(mm::MomentMatrix, pt::Point)
#     mm_eval = SortedDict{Tuple{Exponent, Exponent}, AbstractPolynomial}()
#     for (key, p) in mm.mm
#         res = evaluate(p, pt)
#         if res == Polynomial()
#             delete!(mm_eval, key)
#         else
#             mm_eval[key] = res
#         end
#     end
#     return MomentMatrix(mm_eval, setdiff(mm.vars, SortedSet(keys(pt))), mm.order, mm.matrixkind)
# end


"""
    momentrelaxation = MomentRelaxationPb(relax_ctx, problem, moment_param::SortedDict{String, Tuple{SortedSet{String}, Int}}, max_cliques::SortedDict{String, SortedSet{Variable}})

    Compute the `momentrelaxation` of `problem` corresponding to the clique decomposition `max_cliques` and parameters `moment_param`.
"""
function MomentRelaxationPb(relax_ctx, problem, momentmat_param::SortedDict{String, Int},
                                                localizingmat_param::SortedDict{String, Tuple{SortedSet{String}, Int}},
                                                max_cliques::SortedDict{String, SortedSet{Variable}})
# function MomentRelaxationPb(relax_ctx, problem, momentmat_param::SortedDict{String, Int}, localizingmat_param::SortedDict{String, Tuple{SortedSet{String}, Int}}, max_cliques::SortedDict{String, SortedSet{Variable}})
    println("\n=== MomentRelaxationPb(relax_ctx, problem, moment_param::SortedDict{String, Tuple{SortedSet{String}, Int}}, max_cliques::SortedDict{String, SortedSet{Variable}})")
    println("Compute the moment and localizing matrices associated with the problem constraints and clique decomposition and return a MomentRelaxationPb object.")

    var_to_cliques = SortedDict{Variable, SortedSet{String}}()
    for (clique, vars) in max_cliques
        for var in vars
            haskey(var_to_cliques, var) || (var_to_cliques[var] = SortedSet{String}())
            push!(var_to_cliques[var], clique)
        end
    end

    ## Building linear-in-moments objective
    objective = SortedDict{Moment, Number}()
    for (expo, val) in problem.objective
        clique = get_exponentclique(expo, var_to_cliques)
        objective[Moment(expo, clique)] = val
    end


    ## Building linear matrix inequalities
    momentmatrices = SortedDict{Tuple{String, String}, MomentMatrix}()

    ## Build moment matrix
    for (cliquename, vars) in max_cliques
        dcl = momentmat_param[cliquename]
        momentmatrices[(get_momentcstrname(), cliquename)] = MomentMatrix(relax_ctx, vars, dcl, relax_ctx.symmetries,
                                                                                                relax_ctx.cstrtypes[get_momentcstrname()],
                                                                                                default_clique = cliquename)
    end

    ## Build localizing matrices
    for (cstrname, cstr) in problem.constraints

        cstrtype = get_cstrtype(cstr)
        if cstrtype == :ineqdouble
            cstrname_lo, cstrname_up = get_cstrname(cstrname, cstrtype)

            # Deal with lower inequality
            clique_keys, order = localizingmat_param[cstrname_lo]
            vars, cliquename = collect_cliquesvars(clique_keys, max_cliques)
            length(clique_keys) == 1 || error("MomentRelaxationPb(): constraint $cstrname spans several cliques ($clique_keys).\nNot supported yet.")

            mmt = MomentMatrix(relax_ctx, vars, order, relax_ctx.symmetries,
                                                       relax_ctx.cstrtypes[cstrname_lo],
                                                       var_to_cliques = var_to_cliques)
            # print_with_color(:green, "$cstrname, :Lo\n") ##NOTE: find better logging system.
            product!(mmt, cstr.p - cstr.lb, var_to_cliques)
            momentmatrices[(cstrname_lo, cliquename)] = mmt

            # Deal with upper inequality
            clique_keys, order = localizingmat_param[cstrname_up]
            vars, cliquename = collect_cliquesvars(clique_keys, max_cliques)
            length(clique_keys) == 1 || error("MomentRelaxationPb(): constraint $cstrname spans several cliques ($clique_keys).\nNot supported yet.")

            mmt = MomentMatrix(relax_ctx, vars, order, relax_ctx.symmetries,
                                                       relax_ctx.cstrtypes[cstrname_up],
                                                       var_to_cliques = var_to_cliques)
            # print_with_color(:green, "$cstrname, :Up\n") ##NOTE: find better logging system.
            product!(mmt, cstr.ub - cstr.p, var_to_cliques)
            momentmatrices[(cstrname_up, cliquename)] = mmt

            # # Deal with upper inequality, no recomputing of variables or moment matrix if possible
            # clique_keys_up, order_up = localizingmat_param[cstrname_up]
            # if collect(clique_keys) != collect(clique_keys_up)
            #     warn("clique keys different from lower and upper side of double constraint")
            #     length(clique_keys_up) == 1 || error("MomentRelaxationPb(): constraint $cstrname spans several cliques ($clique_keys).\nNot supported yet.")
            #     vars, cliquename = collect_cliquesvars(clique_keys_up, max_cliques)

            #     mmt = MomentMatrix(relax_ctx, vars, order_up, relax_ctx.symmetries,
            #                                                   relax_ctx.cstrtypes[cstrname_hi],
            #                                                   var_to_cliques = var_to_cliques)
            # elseif order_up != order
            #     warn("order different from lower and upper side of double constraint")
            #     mmt = MomentMatrix(relax_ctx, vars, order_up, relax_ctx.symmetries,
            #                                                   relax_ctx.cstrtypes[cstrname_hi],
            #                                                   var_to_cliques = var_to_cliques)
            # end


        else
            # either cstrtype == :ineqlo, :ineqhi, :eq
            clique_keys, order = localizingmat_param[get_cstrname(cstrname, cstrtype)]
            vars, cliquename = collect_cliquesvars(clique_keys, max_cliques)
            # length(clique_keys) == 1 || error("MomentRelaxationPb(): constraint $cstrname spans several cliques ($clique_keys).\nNot supported yet.")

            mmt = MomentMatrix(relax_ctx, vars, order, relax_ctx.symmetries,
                                                       relax_ctx.cstrtypes[get_cstrname(cstrname, cstrtype)],
                                                       var_to_cliques = var_to_cliques)


            # print_with_color(:green, "$cstrname, :Up\n") ##NOTE: find better logging system.
            product!(mmt, get_normalizedpoly(cstr, cstrtype), var_to_cliques)

            momentmatrices[(get_cstrname(cstrname, cstrtype), cliquename)] = mmt
        end
    end

    ## Locate clique overlapping moments
    expo_to_cliques = SortedDict{Exponent, SortedSet{String}}()

    # Collect Exponents per clique (moment matrix)
    for ((ctrobj, clique), mmtmat) in momentmatrices
        if ctrobj == get_momentcstrname()

            for (key, momentpoly) in mmtmat.mm
                for (moment, coeff) in momentpoly
                    expo = product(moment.conj_part, moment.expl_part)
                    haskey(expo_to_cliques, expo) || (expo_to_cliques[expo] = SortedSet{String}())
                    push!(expo_to_cliques[expo], moment.clique)
                end
            end
        end
    end

    nb_expos = length(expo_to_cliques)
    for (expo, cliques) in expo_to_cliques
        length(cliques) > 1 || delete!(expo_to_cliques, expo)
    end

    nb_overlap_expos = length(expo_to_cliques)
    if nb_overlap_expos > 0
        info("Nb exponents coupled: $nb_overlap_expos (over $nb_expos)")
    end

    return MomentRelaxationPb(objective, momentmatrices, expo_to_cliques)
end


function print(io::IO, momentrelax::MomentRelaxationPb)
    println(io, "Moment Relaxation Problem:")
    println(io, "▶ Objective: ")
    momentlen = maximum(x->length(string(x)), keys(momentrelax.objective))
    for (moment, coeff) in momentrelax.objective
        print_string(io, string(moment), momentlen)
        println(io, " $coeff")
    end

    println(io, "▶ Constraints:")
    for ((cstrname, blocname), mmtmat) in momentrelax.constraints
        println(io, " → $cstrname, $blocname")
        println(io, mmtmat)
    end

    println(io, "▶ Moments clique overlap:")
    if length(momentrelax.moments_overlap) > 0
        mmtlength = maximum(x->length(string(x)), keys(momentrelax.moments_overlap))
        for (moment, cliquenames) in momentrelax.moments_overlap
            print(io, " → ")
            print_string(io, string(moment), mmtlength)
            for clique in cliquenames print(io, "$clique, ") end
            @printf(io, "\b\b \n")
        end
    else
        print(io, "  None")
    end
end