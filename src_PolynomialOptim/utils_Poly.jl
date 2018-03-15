"""
    exponents = get_exponents(variables, dmax::Int)

Compute the set of all exponents in `variables` variables, of degree up to
`dmax`.
"""
function compute_exponents(variables, dmax::Int; compute_conj=false)
    cur_order = Set{Exponent}([Exponent()])
    result = copy(cur_order)
    prev_order = Set{Exponent}()
    for i=1:dmax
        prev_order = copy(cur_order)
        cur_order = Set{Exponent}()
        for var in variables
            if compute_conj
                union!(cur_order, Set([product(Exponent(Dict(var=>Degree(0,1))), elt) for elt in prev_order]))
            else
                union!(cur_order, Set([product(Exponent(var), elt) for elt in prev_order]))
            end
        end
        union!(result, cur_order)
    end
    return result
end
