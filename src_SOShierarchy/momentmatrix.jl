"""
    mm = MomentMatrix(vars::OrderedSet{Variable}, d, symmetries)

    Build the moment matrix corresponding to the moment of degree up to `d` of the `vars` polynomial algebra. 
    Only monomials featuring all `symmetries` appear in the moment matrix.
"""
function MomentMatrix(relax_ctx, vars::OrderedSet{Variable}, d::Int, symmetries::OrderedSet{DataType})
    mm = SortedDict{Tuple{Exponent, Exponent}, AbstractPolynomial}()
    realexpos = compute_exponents(vars, d)
    conjexpos = compute_exponents(vars, d, compute_conj=true)
    for cexp in conjexpos
        for rexp in realexpos
            expo = cexp*rexp
            hassyms = true
            for sym in symmetries
                hassyms = hassyms && has_symmetry(relax_ctx, expo, sym)
            end
            if hassyms
                mm[(cexp, rexp)] = cexp*rexp
            end
        end
    end
    return MomentMatrix(mm, copy(vars), d)
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
    return MomentMatrix(mm_eval, setdiff(mm.vars, OrderedSet(keys(pt))), mm.order)
end


"""
    momentmatrices = compute_momentmat(problem, max_cliques, cliquevarsbycstr, orderbyclique, relax_ctx)

    Compute the moment and localizing matrices associated with the problem constraints and clique decomposition.
"""

function MomentRelaxationPb(relax_ctx, problem, moment_param::SortedDict{String, Tuple{OrderedSet{String}, Int}}, max_cliques::SortedDict{String, OrderedSet{Variable}})
    println("\n=== MomentRelaxationPb(relax_ctx, problem, moment_param::SortedDict{String, Tuple{OrderedSet{String}, Int}}, max_cliques::SortedDict{String, OrderedSet{Variable}})")
    println("Compute the moment and localizing matrices associated with the problem constraints and clique decomposition and return a MomentRelaxationPb object.")

    momentmatrices = SortedDict{Tuple{String, String}, MomentMatrix}()

    for (cstrname, (clique_keys, order)) in moment_param
        # Collect variables involved in constraint
        vars = OrderedSet{Variable}()
        blocname = ""
        for clique_key in clique_keys
            union!(vars, max_cliques[clique_key])
            blocname = blocname*clique_key*"_"
        end
        momentmatrices[(cstrname, blocname[1:end-1])] = MomentMatrix(relax_ctx, vars, order, relax_ctx.symmetries) * problem.constraints[cstrname].p
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