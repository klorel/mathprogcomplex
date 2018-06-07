ROOT=pwd()
include(joinpath("src_SOShierarchy", "SOShierarchy.jl"))

using BenchmarkTools

@benchmark (dict = Dict([i=>rand() for i=1:1000000]))

dict = Dict([i=>rand() for i=1:1000000])



function add_simple(d1, indmax)
    for i=1:indmax
        !haskey(d1, i) && (d1[i] = 0.0)
        d1[i] += 5.0
    end
end

function add_simple_inbnds(d1, indmax)
    @propagate_inbounds for i=1:indmax
        !haskey(d1, i) && (d1[i] = 0.0)
        d1[i] += 5.0
    end
end

function add_custom(d1, indmax)
    for i=1:indmax
        addindex!(d1, 5.0, i)
    end
end

## Validation
# println("Validation")
# println("d1:")
d1 = deepcopy(dict);
add_simple(d1, 2*1000000);
# compare(d1, dict)

# println("\nd2:")
d2 = deepcopy(dict);
add_custom(d2, 2*1000000);
# compare(d2, dict)

@show d1 == d2

println("\nTimings:")
println("Simple")
d1 = deepcopy(dict);
@benchmark add_simple(d1, 1*1000000)

println("\ncustom")
d2 = deepcopy(dict);
@benchmark add_custom(d2, 1*1000000)
