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

ex1x2 = Exponent(SortedDict{Variable, Degree}(x1=>Degree(1,0),x2=>Degree(1,0)))
@show ex1x2.expo
@show ex1x2.degree
exbis = product(Exponent(x1), Exponent(x2))
@show exbis.expo
@show exbis.degree

unsortarray = [Exponent(), ex1x1, ex2x2, ex1x2, ex1, ex2]
x = SortedSet{Exponent}(unsortarray)
@show x
@show typeof(x)
@show SortedSet{Exponent}(unsortarray, ExpoOrdering())

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
    @show typeof(result)
    for e in result print(" $e") end
    println()
    @show typeof(cur_order)
    for e in cur_order print(" $e") end
    println("n$(typeof(cur_order)) \n")
end
result

println("result is $result")
@show typeof(result)

@show SortedSet{Exponent}(collect(result))

println("unsortarray =")
println(unsortarray)

println("\nSortedSet(shuffle(unsortarray))")
println(typeof(unsortarray))
sharr = shuffle(unsortarray)
for e in sharr print(" $e") end; println()
for e in SortedSet(sharr, ExpoOrdering()) print(" $e") end; println(); println()
sharr = shuffle(unsortarray)
for e in sharr print(" $e") end; println()
for e in SortedSet(sharr, ExpoOrdering()) print(" $e") end; println(); println()
sharr = shuffle(unsortarray)
for e in sharr print(" $e") end; println()
for e in SortedSet(sharr, ExpoOrdering()) print(" $e") end; println(); println()
sharr = shuffle(unsortarray)
for e in sharr print(" $e") end; println()
for e in SortedSet(sharr, ExpoOrdering()) print(" $e") end; println(); println()
sharr = shuffle(unsortarray)
for e in sharr print(" $e") end; println()
for e in SortedSet(sharr, ExpoOrdering()) print(" $e") end; println(); println()


println("collect(result)")
y = collect(result)
println(y, " ", typeof(y))

println("\nSortedSet(y, ExpoOrdering())")
for e in SortedSet(y, ExpoOrdering()) print(" $e") end

sharr = shuffle(y)
println(typeof(sharr))
for e in sharr print(" $e") end; println()
for e in SortedSet(sharr, ExpoOrdering()) print(" $e") end; println(); println()
sharr = shuffle(y)
for e in sharr print(" $e") end; println()
for e in SortedSet(sharr, ExpoOrdering()) print(" $e") end; println(); println()
sharr = shuffle(y)
for e in sharr print(" $e") end; println()
for e in SortedSet(sharr, ExpoOrdering()) print(" $e") end; println(); println()
sharr = shuffle(y)
for e in sharr print(" $e") end; println()
for e in SortedSet(sharr, ExpoOrdering()) print(" $e") end; println(); println()
sharr = shuffle(y)
for e in sharr print(" $e") end; println()
for e in SortedSet(sharr, ExpoOrdering()) print(" $e") end; println(); println()

for exp1 in SortedSet(y), exp2 in SortedSet(y)
    println("isless($exp1, $exp2) = $(isless(exp1,exp2))")
end
