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
function MomentRelaxationPb(relax_ctx, problem, momentmat_param::SortedDict{String, Int}, localizingmat_param::SortedDict{String, Tuple{SortedSet{String}, Int}}, max_cliques::SortedDict{String, SortedSet{Variable}})
    println("\n=== MomentRelaxationPb(relax_ctx, problem, moment_param::SortedDict{String, Tuple{SortedSet{String}, Int}}, max_cliques::SortedDict{String, SortedSet{Variable}})")
    println("Compute the moment and localizing matrices associated with the problem constraints and clique decomposition and return a MomentRelaxationPb object.")

    momentmatrices = SortedDict{Tuple{String, String}, MomentMatrix}()

    ## Build moment matrix
    for (cliquename, vars) in max_cliques
        dcl = momentmat_param[cliquename]
        momentmatrices[(get_momentcstrname(), cliquename)] = MomentMatrix(relax_ctx, vars, dcl, relax_ctx.symmetries, relax_ctx.cstrtypes[get_momentcstrname()])
    end

    ## Build localizing matrices
    for (cstrname, cstr) in problem.constraints

        cstrtype = get_cstrtype(cstr)
        if cstrtype == :ineqdouble
            cstrname_lo, cstrname_up = get_cstrname(cstrname, cstrtype)

            # Deal with lower inequality
            clique_keys, order = localizingmat_param[cstrname_lo]
            vars, cliquename = collect_cliquesvars(clique_keys, max_cliques)

            mmt = MomentMatrix(relax_ctx, vars, order, relax_ctx.symmetries, relax_ctx.cstrtypes[cstrname_lo])
            momentmatrices[(cstrname_lo, cliquename)] = mmt * (cstr.p - cstr.lb)

            # Deal with upper inequality, no recomputing of variables or moment matrix if possible
            clique_keys_up, order_up = localizingmat_param[cstrname_up]
            if collect(clique_keys) != collect(clique_keys_up)
                warn("clique keys different from lower and upper side of double constraint")
                vars, cliquename = collect_cliquesvars(clique_keys_up, max_cliques)

                mmt = MomentMatrix(relax_ctx, vars, order_up, relax_ctx.symmetries, relax_ctx.cstrtypes[cstrname_hi])
            elseif order_up != order
                warn("order different from lower and upper side of double constraint")
                mmt = MomentMatrix(relax_ctx, vars, order_up, relax_ctx.symmetries, relax_ctx.cstrtypes[cstrname_hi])
            end

            momentmatrices[(cstrname_up, cliquename)] = mmt * (cstr.ub - cstr.p)

        else
            # either cstrtype == :ineqlo, :ineqhi, :eq
            clique_keys, order = localizingmat_param[get_cstrname(cstrname, cstrtype)]
            vars, cliquename = collect_cliquesvars(clique_keys, max_cliques)

            mmt = MomentMatrix(relax_ctx, vars, order, relax_ctx.symmetries, relax_ctx.cstrtypes[get_cstrname(cstrname, cstrtype)])
            momentmatrices[(get_cstrname(cstrname, cstrtype), cliquename)] = mmt * get_normalizedpoly(cstr, cstrtype)
        end
    end

    ## Locate clique overlapping moments
    moments_overlap = SortedDict{Exponent, SortedSet{String}}()

    # Collect Exponents per clique (moment matrix)
    clique_to_moments = SortedDict{String, SortedSet{Exponent}}()
    for ((ctrobj, clique), mmtmat) in momentmatrices
        if ctrobj == get_momentcstrname()
            cur_expos = SortedSet{Exponent}()

            for (key, moment) in mmtmat.mm
                @assert typeof(first(moment)[1]) == Exponent
                push!(cur_expos, first(moment)[1])
            end
            clique_to_moments[clique] = cur_expos
        end
    end

    cliques = SortedSet{String}(keys(max_cliques))
    for clique_i in cliques, clique_j in cliques
        if clique_i < clique_j
            expos_i = clique_to_moments[clique_i]
            expos_j = clique_to_moments[clique_j]

            for expo in intersect(expos_i, expos_j)
                haskey(moments_overlap, expo) || (moments_overlap[expo] = SortedSet{String}())
                union!(moments_overlap[expo], [clique_i, clique_j])
            end
        end
    end

    return MomentRelaxationPb(problem.objective, momentmatrices, moments_overlap)
end


function print(io::IO, momentrelax::MomentRelaxationPb)
    println(io, "Moment Relaxation Problem:")
    println(io, "▶ Objective: ", momentrelax.objective)
    println(io, "▶ Constraints:")
    for ((cstrname, blocname), mmtmat) in momentrelax.constraints
        println(io, " → $cstrname, $blocname")
        println(io, mmtmat)
    end
    println(io, "▶ Moments clique overlap:")
    mmtlength = maximum(x->length(string(x)), keys(momentrelax.moments_overlap))
    if length(momentrelax.moments_overlap) > 0
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