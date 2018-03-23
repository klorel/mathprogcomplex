include("pb_setting.jl")

# """
#     di = set_uniformdi(problem, d)
#
#     Return a `Dict{String, Int}` corresponding to the relaxation degree of each constraint, uniformly initialized to `d` over constraints.
# """
# function set_uniformdi(problem::Problem, d::Int)
#     println("\n=== set_uniformdi(problem::Problem, d::Int)")
#     println("Return a `Dict{String, Int}` corresponding to the relaxation degree of each constraint, uniformly initialized to `d` over constraints.")
#     return Dict{String, Int}([(cstrname, d) for cstrname in keys(problem.constraints)])
# end
#
#
# """
#     ki = get_ki(pb::Problem)
#
#     Compute the degree of each constraint of `pb`.
# """
# function get_ki(problem::Problem)
#     println("\n=== get_ki(problem::Problem)")
#     println("Compute the degree of each constraint of `pb`.")
#
#     ki = Dict{String, Int}()
#     for (cstrname, cstr) in problem.constraints
#         deg = cstr.p.degree
#         ki[cstrname] = max(deg.explvar, deg.conjvar)
#     end
#
#     println("-> Nb of variables:            $(length(problem.variables))")
#     println("-> Nb of poly constraints:     $(length(problem.constraints))")
#     return ki
# end
#
# """
#     check_di_ki!(di, ki)
#
#     Enforce the condition di - ki ≥ 0 on each constraint (if needed).
# """
# function check_di_ki!(di, ki)
#     println("\n=== check_di_ki!(di, ki)")
#     println("Enforce the condition di - ki ≥ 0 on each constraint (if needed).")
#
#     excessorders = 0
#     for (cstrname, d_cstr) in di
#         if (d_cstr - ki[cstrname] < 0)
#             d_cstr = ki[cstrname]
#         end
#         (d_cstr > ki[cstrname]) && excessorders += 1
#     end
#
#     println("-> number of di > ki:        xx / xx constraints")
# end



"""
    SparsityPattern

    Type for storing and working on sparsitty patterns.
"""
type SparsityPattern end

"""
    sparsity_pattern = compute_sparsitypattern(problem, di, ki)

    Compute the sparsity_pattern corresponding to the given partial orders.
"""
function compute_sparsitypattern(problem::Problem, relax_ctx)
    println("\n=== compute_sparsitypattern(problem::Problem, relax_ctx)")
    println("Compute the sparsity_pattern corresponding to the given partial orders.")
    println("-> Nb vertex / edges:      xx / xx")
    println("-> Order by vertex:        xx / xx (mean/std)")
    sparsity_pattern = SparsityPattern()
    return sparsity_pattern::SparsityPattern
end

"""
    compute_chordalextension!(sparsity_pattern)

    Compute the chordal extension on the provided sparsity_pattern.
"""
function compute_chordalextension!(sparsity_pattern::SparsityPattern)
    println("\n=== compute_chordalextension!(sparsity_pattern::SparsityPattern)")
    println("Compute the chordal extension on the provided sparsity_pattern.")
    println("-> Nb added edges:                 xx")
    println("-> Nb vertex / edges:              xx / xx")
    println("-> Order by vertex:                xx / xx (mean/std)")
end

"""
    maxcliques = compute_maxcliques(sparsity_pattern)

    Compute a `Array{Set{Variable}}` describing the maximum cliques on the provided sparsity_pattern.
"""
function compute_maxcliques(sparsity_pattern::SparsityPattern)
    println("\n=== compute_maxcliques(sparsity_pattern::SparsityPattern)")
    println("Compute a `Array{Set{Variable}}` describing the maximum cliques on the provided sparsity_pattern.")
    println("-> Nb cliques / nb vertex:         xx / xx")
    println("-> Nb vertex by clique:            xx / xx (mean/std)")
    maxcliques = Array{Set{Variable}}()
    return maxcliques::Array{Set{Variable}}
end


"""
    varsbycstr = compute_varsbycstr(problem)

    Compute a `Dict{String, Set{Variable}}` providing the set of variables involved in each constraint.
"""
function compute_varsbycstr(problem::Problem)
    println("\n=== compute_varsbycstr(problem::Problem)")
    println("Compute a `Dict{String, Set{Variable}}` providing the set of variables involved in each constraint.")
    println("-> Nb poly constraints:                xx")
    println("-> Nb variables by poly constraint:    xx / xx (mean/std)")
    varsbycstr = Dict{String, Set{Variable}}()
    return varsbycstr::Dict{String, Set{Variable}}
end

"""
    cliquevarsbycstr = compute_varsbycstr(sparsity_pattern, max_cliques, varsbycstr)

    Compute a `Dict{String, Set{Variable}}` providing the set of variables involved in the SDP localizing matrix corresponding to each constraint.
"""
function compute_varsbycstr(sparsity_pattern, max_cliques, varsbycstr)
    println("\n=== compute_varsbycstr(sparsity_pattern, max_cliques, varsbycstr)")
    println("Compute a `Dict{String, Set{Variable}}` providing the set of variables involved in the SDP localizing matrix corresponding to each constraint.")
    println("-> Nb SDP constraints:                 xx")
    println("-> Nb variables by SDP constraint:     xx / xx (mean/std)")
    cliquevarsbycstr = Dict{String, Set{Variable}}()
    return cliquevarsbycstr::Dict{String, Set{Variable}}
end

"""
    orderbyclique = compute_cliqueorders(sparsity_pattern, di, varsbycstr, max_cliques)

    Compute a `Array{Int}` providing the relaxation order corresponding to each clique.
"""
function compute_cliqueorders(sparsity_pattern, varsbycstr, max_cliques, relax_ctx)
    println("\n=== compute_cliqueorders(sparsity_pattern, varsbycstr, max_cliques, relax_ctx)")
    println("Compute a `Array{Int}` providing the relaxation order corresponding to each clique.")
    println("-> Nb cliques:                         xx")
    println("-> Relaxation order by clique:         xx / xx (mean/std)")
    orderbyclique = Dict{Int, Int}()
    return orderbyclique::Dict{Int, Int}
end


include("momentmatrix_dense.jl")
include("compute_Bi.jl")


include("build_SDP_SOS.jl")

include("export_JuMP.jl")


function print_cmat(mat::AbstractArray, round = 1e-3)
    for i=1:size(mat, 1)
        for j=1:size(mat, 2)
            re, im = real(mat[i, j]), imag(mat[i, j])
            @printf("% 5.4f", re)
            @printf(" ")
            @printf("%+5.4fim", im)
            @printf("   ")
        end
        @printf("\n")
    end
end