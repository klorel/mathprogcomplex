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
    return MomentMatrix(mm, SortedSet(vars), d, matrixkind)::MomentMatrix
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
    return 0
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

    return resmpoly::SortedDict{Moment, Number}
end

function product(moment::Moment, expo::Exponent, var_to_cliques)
    resexpo = product(moment.conj_part, moment.expl_part)
    product!(resexpo, expo)

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
