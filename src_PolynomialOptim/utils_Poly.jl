"""
    exponents = get_exponents(variables, dmax::Int)

Compute the set of all exponents in `variables` variables, of degree up to
`dmax`.
"""
function compute_exponents(variables::SortedSet{Variable}, dmax::Int; compute_conj=false)
    cur_order = Set{Exponent}([Exponent()])
    result = copy(cur_order)
    prev_order = Set{Exponent}()
    for i=1:dmax
        prev_order = copy(cur_order)
        cur_order = Set{Exponent}()
        for var in variables
            if compute_conj
                union!(cur_order, Set([product(conj(var), elt) for elt in prev_order]))
            else
                union!(cur_order, Set([product(Exponent(var), elt) for elt in prev_order]))
            end
        end
        union!(result, cur_order)
    end
    return result
end


"""
    ishomo = is_homogeneous(p, kind)

    Check wether `p`, of kind `:Real` or `:Complex` is homogeneous or not.
"""
function is_homogeneous(p::Polynomial, kind::Symbol)
    ishomo = true
    for (expo, λ) in p
        ishomo = ishomo && is_homogeneous(expo, kind)
    end
    return ishomo
end

"""
    ishomo = is_homogeneous(expo, kind)

    Check wether `expo`, of kind `:Real` or `:Complex` is homogeneous or not.
"""
function is_homogeneous(expo::Exponent, kind::Symbol)
    explsum, conjsum = get_sumdegs(expo)
    if kind == :Real
        (conjsum != 0) && error("is_homogeneous(): Exponent $expo has conjugated variables for a real hierarchy.")
        return explsum % 2 == 0
    elseif kind == :Complex
        return explsum == conjsum
    else
        error("is_homogeneous(expo, kind): kind should be either :Real or :Complex ($kind here).")
    end
end

"""
    explsum, conjsum = get_sumdegs(expo)

    Compute `|α|`, `|β|` the sum of the real variables exponents and conjugated variables exponents.
"""
function get_sumdegs(expo::Exponent) 
    explsum = conjsum = 0
    for (var, deg) in expo
        explsum += deg.explvar
        conjsum += deg.conjvar
    end
    return explsum, conjsum
end

"""
    isdoublesided(cstr::Constraint)

    Test whether the constraint has two finite bounds.
"""
function isdoublesided(cstr::Constraint)
    return (cstr.lb != -Inf-im*Inf) && (cstr.ub != Inf+im*Inf)
end