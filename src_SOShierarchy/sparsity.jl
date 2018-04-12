"""
    moments_param = build_sparsity(relax_ctx, problem)

    Build the sparsitty pattern and variables decomposition for laying out the moment or SOS hierarchy
"""
function build_sparsity(relax_ctx, problem, max_cliques::SortedDict{String, OrderedSet{Variable}})

    if relax_ctx.issparse == false
        (length(max_cliques) == 1) || error("build_sparsity(): Relaxation is not sparse, one clique is expected (not $(length(max_cliques)))")
        
        moments_param = SortedDict{String, Tuple{OrderedSet{String}, Int}}()
        for (cstr, di) in relax_ctx.di
            ki = relax_ctx.ki[cstr]
            moments_param[cstr] = (OrderedSet(["clique1"]), di-ki)
        end
        return moments_param

    else
        error("build_sparsity(): Sparse hierarchy not handled yet.")
        
        # TODO
    end
end



function get_maxcliques(relax_ctx, problem)
    if !relax_ctx.issparse
        vars = OrderedSet{Variable}([Variable(name, kind) for (name, kind) in problem.variables])
        return SortedDict{String, OrderedSet{Variable}}("clique1"=>vars)
    else
        error("Sparse relaxation is not supported yet")
    end
end

#################################################################################
## Old stuff


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

    Compute a `Array{OrderedSet{Variable}}` describing the maximum cliques on the provided sparsity_pattern.
"""
function compute_maxcliques(sparsity_pattern::SparsityPattern)
    println("\n=== compute_maxcliques(sparsity_pattern::SparsityPattern)")
    println("Compute a `Array{OrderedSet{Variable}}` describing the maximum cliques on the provided sparsity_pattern.")
    println("-> Nb cliques / nb vertex:         xx / xx")
    println("-> Nb vertex by clique:            xx / xx (mean/std)")
    maxcliques = Array{OrderedSet{Variable}}()
    return maxcliques::Array{OrderedSet{Variable}}
end


"""
    varsbycstr = compute_varsbycstr(problem)

    Compute a `SortedDict{String, OrderedSet{Variable}}` providing the set of variables involved in each constraint.
"""
function compute_varsbycstr(problem::Problem)
    println("\n=== compute_varsbycstr(problem::Problem)")
    println("Compute a `SortedDict{String, OrderedSet{Variable}}` providing the set of variables involved in each constraint.")
    println("-> Nb poly constraints:                xx")
    println("-> Nb variables by poly constraint:    xx / xx (mean/std)")
    varsbycstr = SortedDict{String, OrderedSet{Variable}}()
    return varsbycstr::SortedDict{String, OrderedSet{Variable}}
end

"""
    cliquevarsbycstr = compute_varsbycstr(sparsity_pattern, max_cliques, varsbycstr)

    Compute a `SortedDict{String, OrderedSet{Variable}}` providing the set of variables involved in the SDP localizing matrix corresponding to each constraint.
"""
function compute_varsbycstr(sparsity_pattern, max_cliques, varsbycstr)
    println("\n=== compute_varsbycstr(sparsity_pattern, max_cliques, varsbycstr)")
    println("Compute a `SortedDict{String, OrderedSet{Variable}}` providing the set of variables involved in the SDP localizing matrix corresponding to each constraint.")
    println("-> Nb SDP constraints:                 xx")
    println("-> Nb variables by SDP constraint:     xx / xx (mean/std)")
    cliquevarsbycstr = SortedDict{String, OrderedSet{Variable}}()
    return cliquevarsbycstr::SortedDict{String, OrderedSet{Variable}}
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
    orderbyclique = SortedDict{Int, Int}()
    return orderbyclique::SortedDict{Int, Int}
end