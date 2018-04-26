"""
    mm = MomentMatrix(vars::SortedSet{Variable}, d, symmetries)

    Build the moment matrix corresponding to the moment of degree up to `d` of the `vars` polynomial algebra.
    Only monomials featuring all `symmetries` appear in the moment matrix.
"""
function MomentMatrix(relax_ctx, vars::SortedSet{Variable}, d::Int, symmetries::SortedSet{DataType})
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
    return MomentMatrix(mm, SortedSet(vars), d)
end

function copy(mm::MomentMatrix)
    return MomentMatrix(copy(mm.mm), mm.vars, mm.order)
end

function print(io::IO, mm::MomentMatrix)
    for (key, val) in mm.mm
        println(io, "($(key[1]), $(key[2])) ⟶  $val")
    end
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
    return MomentMatrix(mm_eval, setdiff(mm.vars, SortedSet(keys(pt))), mm.order)
end


"""
    momentrelaxation = MomentRelaxationPb(relax_ctx, problem, moment_param::SortedDict{String, Tuple{SortedSet{String}, Int}}, max_cliques::SortedDict{String, SortedSet{Variable}})

    Compute the `momentrelaxation` of `problem` corresponding to the clique decomposition `max_cliques` and parameters `moment_param`.
"""
function MomentRelaxationPb(relax_ctx, problem, moment_param::SortedDict{String, Tuple{SortedSet{String}, Int}}, max_cliques::SortedDict{String, SortedSet{Variable}})
    println("\n=== MomentRelaxationPb(relax_ctx, problem, moment_param::SortedDict{String, Tuple{SortedSet{String}, Int}}, max_cliques::SortedDict{String, SortedSet{Variable}})")
    println("Compute the moment and localizing matrices associated with the problem constraints and clique decomposition and return a MomentRelaxationPb object.")

    momentmatrices = SortedDict{Tuple{String, String}, MomentMatrix}()

    ## Build moment matrix
    # NOTE: sparsity work tbd here : several moment matrices ?
    clique_keys, order = moment_param[get_momentcstrname()]
    vars, cliquename = collect_cliquesvars(clique_keys, max_cliques)
    
    momentmatrices[(get_momentcstrname(), cliquename)] = MomentMatrix(relax_ctx, vars, order, relax_ctx.symmetries)

    ## Build localizing matrices
    for (cstrname, cstr) in problem.constraints

        cstrtype = get_cstrtype(cstr)
        if cstrtype == :ineqdouble
            cstrname_lo, cstrname_up = get_cstrname(cstrname, cstrtype)
            
            # Deal with lower inequality
            clique_keys, order = moment_param[cstrname_lo]
            vars, cliquename = collect_cliquesvars(clique_keys, max_cliques)
            
            mmt = MomentMatrix(relax_ctx, vars, order, relax_ctx.symmetries)
            momentmatrices[(cstrname_lo, cliquename)] = mmt * (cstr.p - cstr.lb)

            # Deal with upper inequality, no recomputing of variables or moment matrix if possible
            clique_keys_up, order_up = moment_param[cstrname_up]
            if collect(clique_keys) != collect(clique_keys_up)
                warn("clique keys different from lower and upper side of double constraint")
                vars, cliquename = collect_cliquesvars(clique_keys_up, max_cliques)

                mmt = MomentMatrix(relax_ctx, vars, order_up, relax_ctx.symmetries)
            elseif order_up != order
                warn("order different from lower and upper side of double constraint")
                mmt = MomentMatrix(relax_ctx, vars, order_up, relax_ctx.symmetries)
            end
            
            momentmatrices[(cstrname_up, cliquename)] = mmt * (cstr.ub - cstr.p)

        else
            # either cstrtype == :ineqlo, :ineqhi, :eq
            clique_keys, order = moment_param[get_cstrname(cstrname, cstrtype)]
            vars, cliquename = collect_cliquesvars(clique_keys, max_cliques)
            
            mmt = MomentMatrix(relax_ctx, vars, order, relax_ctx.symmetries)
            momentmatrices[(get_cstrname(cstrname, cstrtype), cliquename)] = mmt * get_normalizedpoly(cstr, cstrtype)
        end
    end


    return MomentRelaxationPb(problem.objective, momentmatrices)
end


function print(io::IO, momentrelax::MomentRelaxationPb)
    println(io, "Moment Relaxation Problem:")
    println(io, "▶ Objective: ", momentrelax.objective)
    println(io, "▶ Constraints:")
    for ((cstrname, blocname), mmtmat) in momentrelax.constraints
        println(io, "--> $cstrname, $blocname")
        println(io, mmtmat)
    end
end