"""
    mm = MomentMatrix(vars::Set{Variable}, d)

    Build the moment matrix corresponding to the moment of degree up to `d` of the `vars` polynomial algebra.
"""
function MomentMatrix(vars::Set{Variable}, d::Int)
    mm = Dict{Tuple{Exponent, Exponent}, AbstractPolynomial}()
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
    mm_eval = Dict{Tuple{Exponent, Exponent}, AbstractPolynomial}()
    for (key, p) in mm.mm
        mm_eval[key] = evaluate(p, pt)
    end
    return MomentMatrix(mm_eval, setdiff(mm.vars, keys(pt)), mm.order)
end


"""
    momentmatrices = compute_momentmat(problem, max_cliques, cliquevarsbycstr, orderbyclique, relax_ctx)

    Compute the moment and localizing matrices associated with the problem constraints and clique decomposition.
"""
function compute_momentmat(relax_ctx, problem, moment_param::Dict{String, Tuple{Set{String}, Int}}, max_cliques::Dict{String, Set{Variable}})
    println("\n=== compute_momentmat(relax_ctx, problem, moment_param, max_cliques)")
    println("Compute the moment and localizing matrices associated with the problem constraints and vlique decomposition.")

    # NOTE: Things will have to be slightly extended to support the several SDP sparse moment constraint (cstr key will not suffise)
    momentmatrices = Dict{String, MomentMatrix}()

    for (cstrname, (clique_keys, order)) in moment_param
        # Collect variables involved in constraint
        vars = Set{Variable}()
        for clique_key in clique_keys
            union!(vars, max_cliques[clique_key])
        end
        momentmatrices[cstrname] = MomentMatrix(vars, order) * problem.constraints[cstrname].p
    end

    println("xxxxx which indicators ? xxxxx")

    return momentmatrices
end
