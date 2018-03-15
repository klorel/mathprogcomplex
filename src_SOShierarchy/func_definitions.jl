"""
    k_i = get_k_i(pb::Problem)

    Compute the degree of each constraint of `pb`.
"""
function get_ki(problem::Problem)
    println("\n=== get_ki(problem::Problem)")
    println("Compute the degree of each constraint of `pb`.")
    println("-> Nb of variables:            xx")
    println("-> Nb of poly constraints:     xx")
    return Dict{String, Int}()
end

"""
    check_di_ki!(d_i, k_i)

    Enforce the condition d_i - k_i ≥ 0 on each constraint (if needed).
"""
function check_di_ki!(d_i, k_i)
    println("\n=== check_di_ki!(d_i, k_i)")
    println("Enforce the condition d_i - k_i ≥ 0 on each constraint (if needed).")
    println("-> number of d_i > k_i:        xx / xx constraints")
end



"""
    SparsityPattern

    Type for storing and working on sparsitty patterns.
"""
type SparsityPattern end

"""
    sparsity_pattern = compute_sparsitypattern(problem, d_i, k_i)

    Compute the sparsity_pattern corresponding to the given partial orders.
"""
function compute_sparsitypattern(problem::Problem, d_i, k_i)
    println("\n=== compute_sparsitypattern(problem::Problem, d_i, k_i)")
    println("Compute the sparsity_pattern corresponding to the given partial orders.")
    println("-> Nb vertex / edges:              xx / xx")
    println("-> Order by vertex mean / std:     xx / xx")
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
    println("-> Nb added edges: xx")
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
    println("-> Nb variables by clique:         xx / xx (mean/std)")
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
    println("-> Nb variables by poly constraints:   xx / xx (mean/std)")
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
    println("-> Nb variables by SDP constraints:    xx / xx (mean/std)")
    cliquevarsbycstr = Dict{String, Set{Variable}}()
    return cliquevarsbycstr::Dict{String, Set{Variable}}
end

"""
    orderbyclique = compute_cliqueorders(sparsity_pattern, d_i, varsbycstr, max_cliques)

    Compute a `Array{Int}` providing the relaxation order corresponding to each clique.
"""
function compute_cliqueorders(sparsity_pattern, d_i, varsbycstr, max_cliques)
    println("\n=== compute_cliqueorders(sparsity_pattern, d_i, varsbycstr, max_cliques)")
    println("Compute a `Array{Int}` providing the relaxation order corresponding to each clique.")
    println("-> Nb cliques:                         xx")
    println("-> Relaxation order by clique:         xx / xx (mean/std)")
    orderbyclique = Dict{Int, Int}()
    return orderbyclique::Dict{Int, Int}
end




"""
    B_i_dict = compute_Bibycstr(problem, max_cliques, cliquevarsbycstr, orderbyclique)

    Compute the decomposition of localizing matrix corresponding to each constraint on the moment variable basis, yielding several matrices B_i,α,β for each constraint i.
"""
function compute_Bibycstr(problem, max_cliques, cliquevarsbycstr, orderbyclique)
    println("\n=== compute_Bibycstr(problem, max_cliques, cliquevarsbycstr, orderbyclique)")
    println("Compute the decomposition of localizing matrix corresponding to each constraint on the moment variable basis, yielding several matrices B_i,α,β for each constraint i.")
    println("-> Nb matrices:                        xx")
    println("-> Nb matrices by constraint:          xx / xx (mean/std)")
    println("-> Sparsity degree of the matrices:    xx / xx (mean/std)")
    return
end


"""
    SDP_SOS = build_SDP_SOS(problem, d_i, max_cliques, B_i, cliquevarsbycstr, orderbyclique)

    Build the primal SDP corresponding to the dual SOS hierarchy for the provided problem.
"""
function build_SDP_SOS(problem, d_i, max_cliques, B_i, cliquevarsbycstr, orderbyclique)
    println("\n=== build_SDP_SOS(problem, d_i, max_cliques, B_i, cliquevarsbycstr, orderbyclique)")
    println("Build the primal SDP corresponding to the dual SOS hierarchy for the provided problem.")
    println("-> Nb of SDP variables:                                xx")
    println("-> Size of the SDP variables:                          xx/xx (mean/std)")
    println("-> Nb of eq. SDP constraints:                          xx")
    println("-> Nb of coupling y_α,β variables:                     xx")
    println("-> Nb of constraint involved by coupling variables:    xx/xx (mean/std)")
    return
end


"""
    m = make_JuMPproblem(SDP_SOS)

    Convert the SDP_SOS problem into a JuMP problem
"""
function make_JuMPproblem(SDP_SOS)
    println("\n=== make_JuMPproblem(SDP_SOS)")
    println("Convert the SDP_SOS problem into a JuMP problem")
    return
end
