"""
    mm = MomentMatrix(vars::SortedSet{Variable}, d)

    Build the moment matrix corresponding to the moment of degree up to `d` of the `vars` polynomial algebra.
"""
function MomentMatrix(vars::SortedSet{Variable}, d::Int)
    mm = SortedDict{Tuple{Exponent, Exponent}, AbstractPolynomial}()
    realexpos = compute_exponents(vars, d)
    conjexpos = compute_exponents(vars, d, compute_conj=true)
    for cexp in conjexpos
        for rexp in realexpos
            mm[(cexp, rexp)] = cexp*rexp
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
    return MomentMatrix(mm_eval, setdiff(mm.vars, SortedSet(keys(pt))), mm.order)
end


"""
    momentmatrices = compute_momentmat(problem, max_cliques, cliquevarsbycstr, orderbyclique, relax_ctx)

    Compute the moment and localizing matrices associated with the problem constraints and clique decomposition.
"""
function compute_momentmat(relax_ctx, problem, moment_param::SortedDict{String, Tuple{SortedSet{String}, Int}}, max_cliques::SortedDict{String, SortedSet{Variable}})
    println("\n=== compute_momentmat(relax_ctx, problem, moment_param, max_cliques)")
    println("Compute the moment and localizing matrices associated with the problem constraints and clique decomposition.")

    # NOTE: Things will have to be slightly extended to support the several SDP sparse moment constraint (cstr key will not suffise)
    momentmatrices = SortedDict{Tuple{String, String}, MomentMatrix}()

    for (cstrname, (clique_keys, order)) in moment_param
        # Collect variables involved in constraint
        vars = SortedSet{Variable}()
        blocname = ""
        for clique_key in clique_keys
            union!(vars, max_cliques[clique_key])
            blocname = blocname*clique_key*"_"
        end
        momentmatrices[(cstrname, blocname[1:end-1])] = MomentMatrix(vars, order) * problem.constraints[cstrname].p
        println("****************************************************************************************")
        println(MomentMatrix(vars, order))
        println(problem.constraints[cstrname].p)
        println(momentmatrices[(cstrname, blocname[1:end-1])])
        println("****************************************************************************************")
    end

    println("xxxxx which indicators ? xxxxx")

    return momentmatrices
end


function MomentRelaxationPb(relax_ctx, problem, moment_param::SortedDict{String, Tuple{SortedSet{String}, Int}}, max_cliques::SortedDict{String, SortedSet{Variable}})
    momentmatrices = compute_momentmat(relax_ctx, problem, moment_param, max_cliques)

    if relax_ctx.leveragesymmetries && has_phasesymmetry(relax_ctx, problem)
        println("MomentRelaxationPb(): Problem is phase-shift invariant. Leveraging to reduce the number of moments.")
        enforce_phaseinvariance!(relax_ctx, momentmatrices)
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