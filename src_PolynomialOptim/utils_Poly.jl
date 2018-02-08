"""
    exponents = get_exponents(variables, dmax::Int)

Compute the set of all exponents in `variables` variables, of degree up to
`dmax`.
"""
function compute_exponents(variables, dmax::Int)
    cur_order = Set{Exponent}([Exponent()])
    result = copy(cur_order)
    prev_order = Set{Exponent}()
    for i=1:d
        prev_order = copy(cur_order)
        cur_order = Set{Exponent}()
        for var in input_vars
            union!(cur_order, Set([product(Exponent(var), elt) for elt in prev_order]))
        end
        union!(result, cur_order)
    end
    return result
end
