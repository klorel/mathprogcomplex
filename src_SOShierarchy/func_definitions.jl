"""
    norm_pb = normalize_problem(problem)

    Transform the problem so it contains only inequalities s.a. g_i(z) ≥ 0, add the moment constraint g_0(z) = 1 ≥ 0.
"""
# NOTE: Care for integer variables will have to be taken at some point.
function normalize_problem(problem)
    println("\n=== normalize_problem(problem)")
    println("Transform the problem so it contains only ineqs s.a. g_i(z) ≥ 0, add the moment constraint g_0(z) = 1 ≥ 0.")

    normpb = Problem()
    normpb.variables = copy(problem.variables)
    normpb.objective = copy(problem.objective)

    for (cstrname, cstr) in get_constraints(problem)
        if cstr.lb != -Inf-im*Inf
            add_constraint!(normpb, cstrname*"_lo", 0 << (cstr.p - cstr.lb))
        end
        if cstr.ub != Inf+im*Inf
            add_constraint!(normpb, cstrname*"_hi", 0 << (cstr.ub - cstr.p))
        end
    end


    exposet = Set()
    nb_expotot = 0
    degbycstr = Int64[]
    for (cstrname, cstr) in get_constraints(normpb)
        push!(degbycstr, max(cstr.p.degree.explvar, cstr.p.degree.conjvar))
        for (expo, λ) in cstr.p
            push!(exposet, expo)
            nb_expotot += 1
        end
    end
    println("-> Problem with:")
    println("-> Nb of variables:                        $(length(normpb.variables))")
    println("-> Nb of constraints:                      $(length(normpb.constraints))")
    println("-> Nb of monomials:                        $nb_expotot ($(length(exposet)) different)")
    println("-> Max total degree by constraint:         $(mean(degbycstr)) / $(std(degbycstr)) (mean/std)")
    return normpb
end

"""
    relax_ctx = set_relaxation(problem; issparse = false, ismultiordered = false, d = 1)

    Build a `relax_ctx` object containing relaxation choices and problem features : order by constraint, relaxation order by constraint.
"""
function set_relaxation(problem; issparse = false, ismultiordered = false, d = 1)
    println("\n=== set_relaxation(problem, issparse = false, ismultiordered = false, d = 1)")
    println("Build a `relax_ctx` object containing relaxation choices and problem features : polynomial degree by constraint, relaxation order by constraint...")

    println("-> Number of constraints s.t. di > ki:     xx")
    println("-> Max total degree on such constraints:   xx / xx (mean/std)")
    println("All variables appearing in such constraints will be linked in the sparsity pattern, which will largely densify it.")
    return
end

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
    SDP_SOS = build_SDP_SOS(problem, di, max_cliques, B_i, cliquevarsbycstr, orderbyclique)

    Build the primal SDP corresponding to the dual SOS hierarchy for the provided problem.
"""
function build_SDP_SOS(problem, max_cliques, B_i, cliquevarsbycstr, orderbyclique, relax_ctx)
    println("\n=== build_SDP_SOS(problem, max_cliques, B_i, cliquevarsbycstr, orderbyclique, relax_ctx)")
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
