include("src_SOShierarchy/SOShierarchy.jl")

x1 = Variable("x1", Real)
x2 = Variable("x2", Real)

import Base.Ordering
import Base.lt
import DataStructures.eq

sset = SortedSet{Variable}([x1, x2])

@show compute_exponents(sset, 2)


immutable ExpoOrdering <: Ordering
end

lt(::Base.Order.ForwardOrdering, a::Exponent, b::Exponent) = isless(a, b)
lt(::ExpoOrdering, a, b) = isless(a, b)
# eq(::CaseInsensitive, a, b) = isequal(lowercase(a), lowercase(b))

ex1 = Exponent(x1)
ex2 = Exponent(x2)
ex1x2 = Exponent(SortedDict{Variable, Degree}(x1=>Degree(1,0),x2=>Degree(1,0)))
ex1x1 = Exponent(SortedDict{Variable, Degree}(x1=>Degree(2,0)))
ex2x2 = Exponent(SortedDict{Variable, Degree}(x2=>Degree(2,0)))

@show SortedSet{Exponent}([Exponent(), ex1x1, ex2x2, ex1x2, ex1, ex2])
@show SortedSet{Exponent}([Exponent(), ex1x1, ex2x2, ex1x2, ex1, ex2], ExpoOrdering())

println("-----------")

dmax = 2
variables = sset
compute_conj=false

cur_order = SortedSet{Exponent}([Exponent()])
result = deepcopy(cur_order)
prev_order = SortedSet{Exponent}()
for i=1:dmax
    prev_order = deepcopy(cur_order)
    cur_order = SortedSet{Exponent}()
    for var in variables
        if compute_conj
            union!(cur_order, SortedSet{Exponent}([product(conj(var), elt) for elt in prev_order]))
        else
            union!(cur_order, SortedSet{Exponent}([product(Exponent(var), elt) for elt in prev_order]))
        end
    end
    # union!(result, cur_order)
    u = union(result, cur_order)
    result = deepcopy(u)
    println(" it $i")
    for e in result print(" $e") end
    println()
    for e in cur_order print(" $e") end
    println("n$(typeof(cur_order)) \n")
    
end
result

println("result is $result")
@show result