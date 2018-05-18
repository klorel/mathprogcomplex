"""
    momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem)
    Build the sparsitty pattern and variables decomposition for laying out the moment or SOS hierarchy
"""
function build_sparsity(relax_ctx, problem, max_cliques::SortedDict{String, SortedSet{Variable}})

    ((relax_ctx.issparse == false) && (length(max_cliques) > 1)) && error("build_sparsity(): Relaxation is not sparse, one clique is expected (not $(length(max_cliques)))")
    localizingmat_param = SortedDict{String, Tuple{SortedSet{String}, Int}}()
    for (ctrname, ctr) in problem.constraints
        ctrtype = get_cstrtype(ctr)
        ctrcliques = get_locctrcliques(ctr.p, max_cliques)

        if ctrtype == :ineqdouble
            ctrname_lo, ctrname_up = get_cstrname(ctrname, ctrtype)
            di_lo, ki_lo = relax_ctx.di[ctrname_lo], relax_ctx.ki[ctrname_lo]
            di_up, ki_up = relax_ctx.di[ctrname_up], relax_ctx.ki[ctrname_up]
            localizingmat_param[ctrname_lo] = (ctrcliques, di_lo-ki_lo)
            localizingmat_param[ctrname_up] = (ctrcliques, di_up-ki_up)
        else # :ineqlo, :ineqhi, :eq
            di, ki = relax_ctx.di[get_cstrname(ctrname, ctrtype)], relax_ctx.ki[get_cstrname(ctrname, ctrtype)]
            localizingmat_param[get_cstrname(ctrname, ctrtype)] = (ctrcliques, di-ki)
        end
    end


    momentmat_param = SortedDict{String, Int}()
    for (cliquename, cliquevars) in max_cliques
        cur_d = -1
        for (ctrname, (ctrcliques, _)) in localizingmat_param
            if length(ctrcliques) == 1 && cliquename == first(ctrcliques)
                cur_d = max(cur_d, relax_ctx.di[ctrname])
            end
        end
        if cur_d == -1 # Clique does not match a full constraint
            cur_d = minimum(values(relax_ctx.di))
        end
        momentmat_param[cliquename] = cur_d
    end

    return momentmat_param, localizingmat_param
end


"""
pvars = get_variables(p)

Collect all variables appearing in `p`.
"""
function get_variables(p::Polynomial)
    pvars = SortedSet{Variable}()
    for (expo, coeff) in p
        for (var, deg) in expo
            insert!(pvars, var)
        end
    end
    return pvars
end

"""
    locctrcliques = get_locctrcliques(ctrvars, max_cliques)

    Find a minimal set of cliques gathering all variables from `ctrvars`.
"""
function get_locctrcliques(p::Polynomial, max_cliques::SortedDict{String, SortedSet{Variable}})
    ctrvars = get_variables(p)

    # Build constraint variables to cliques dict
    var_to_cliques = SortedDict{Variable, SortedSet{String}}()
    for (clique, cliquevars) in max_cliques
        for var in intersect(cliquevars, ctrvars)
            haskey(var_to_cliques, var) || (var_to_cliques[var] = SortedSet{String}())
            insert!(var_to_cliques[var], clique)
        end
    end

    def_cliques = SortedSet{String}()
    unaffected_vars = SortedSet{Variable}(keys(var_to_cliques))

    keepon = true
    i = 0
    while keepon
        i += 1
        # Find which variables appear in one clique only, remove them from unaffected_vars
        for var in unaffected_vars
            cliques = var_to_cliques[var]
            if length(cliques) == 1
                insert!(def_cliques, first(cliques))
                delete!(unaffected_vars, var)
            end
        end

        # Remove variables involved in at least one clique previously selected
        for var in unaffected_vars
            inter = intersect(var_to_cliques[var], def_cliques)
            if !isempty(inter)
                var_to_cliques[var] = SortedSet{String}([first(inter)]) # NOTE: proper way to choose in inter here ?
                delete!(unaffected_vars, var)
            end
        end

        # Hopefully all variables are treated that way. Else repeat this process by choosing a clique. Again, which one ?
        if length(unaffected_vars) != 0
            warn("get_locctrcliques(): length(unaffected_vars) = $(length(unaffected_vars))")
            cliques_from_unaffvar = SortedDict{String, Int}()
            for var in unaffected_vars
                for clique in var_to_cliques[var]
                    haskey(cliques_from_unaffvar, clique) || (cliques_from_unaffvar[clique] = 0)
                    cliques_from_unaffvar[clique] += 1
                end
            end
            cur_clique, cur_pop = first(cliques_from_unaffvar)[1], first(cliques_from_unaffvar)[2]
            for (clique, pop) in cliques_from_unaffvar
                if pop > cur_pop
                    cur_clique = clique
                    cur_pop = pop
                end
            end
            insert!(def_cliques, cur_clique)
        else
            keepon = false
        end
    end

    return def_cliques
end

function get_maxcliques(relax_ctx, problem)
    if !relax_ctx.issparse
        vars = SortedSet{Variable}([Variable(name, kind) for (name, kind) in problem.variables])
        return SortedDict{String, SortedSet{Variable}}("clique1"=>vars)
    else
        error("Sparse relaxation is not supported yet")
    end
end

function get_WB5cliques(relax_ctx, problem)
    maxcliques = SortedDict{String, SortedSet{Variable}}()
    maxcliques["clique1"] = SortedSet{Variable}([
        Variable("BaseCase_1_VOLT_Im", Real),
        Variable("BaseCase_1_VOLT_Re", Real),
        Variable("BaseCase_2_VOLT_Im", Real),
        Variable("BaseCase_2_VOLT_Re", Real),
        Variable("BaseCase_3_VOLT_Im", Real),
        Variable("BaseCase_3_VOLT_Re", Real)])
    maxcliques["clique2"] = SortedSet{Variable}([
        Variable("BaseCase_2_VOLT_Im", Real),
        Variable("BaseCase_2_VOLT_Re", Real),
        Variable("BaseCase_3_VOLT_Im", Real),
        Variable("BaseCase_3_VOLT_Re", Real),
        Variable("BaseCase_4_VOLT_Im", Real),
        Variable("BaseCase_4_VOLT_Re", Real),
        Variable("BaseCase_5_VOLT_Im", Real),
        Variable("BaseCase_5_VOLT_Re", Real)])
    return maxcliques
end

"""
    vars, blocname = collect_cliquesvars(clique_keys, max_cliques)
    Collect variables of `cliques_keys` cliques, described in `max_cliques`
"""
function collect_cliquesvars(clique_keys, max_cliques)
    # Collect variables involved in constraint
    vars = SortedSet{Variable}()
    blocname = ""
    for clique_key in clique_keys
        union!(vars, max_cliques[clique_key])
        blocname = blocname*clique_key*"_"
    end
    return vars, blocname[1:end-1]
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
    Compute a `Array{SortedSet{Variable}}` describing the maximum cliques on the provided sparsity_pattern.
"""
function compute_maxcliques(sparsity_pattern::SparsityPattern)
    println("\n=== compute_maxcliques(sparsity_pattern::SparsityPattern)")
    println("Compute a `Array{SortedSet{Variable}}` describing the maximum cliques on the provided sparsity_pattern.")
    println("-> Nb cliques / nb vertex:         xx / xx")
    println("-> Nb vertex by clique:            xx / xx (mean/std)")
    maxcliques = Array{SortedSet{Variable}}()
    return maxcliques::Array{SortedSet{Variable}}
end


"""
    varsbycstr = compute_varsbycstr(problem)
    Compute a `SortedDict{String, SortedSet{Variable}}` providing the set of variables involved in each constraint.
"""
function compute_varsbycstr(problem::Problem)
    println("\n=== compute_varsbycstr(problem::Problem)")
    println("Compute a `SortedDict{String, SortedSet{Variable}}` providing the set of variables involved in each constraint.")
    println("-> Nb poly constraints:                xx")
    println("-> Nb variables by poly constraint:    xx / xx (mean/std)")
    varsbycstr = SortedDict{String, SortedSet{Variable}}()
    return varsbycstr::SortedDict{String, SortedSet{Variable}}
end

"""
    cliquevarsbycstr = compute_varsbycstr(sparsity_pattern, max_cliques, varsbycstr)
    Compute a `SortedDict{String, SortedSet{Variable}}` providing the set of variables involved in the SDP localizing matrix corresponding to each constraint.
"""
function compute_varsbycstr(sparsity_pattern, max_cliques, varsbycstr)
    println("\n=== compute_varsbycstr(sparsity_pattern, max_cliques, varsbycstr)")
    println("Compute a `SortedDict{String, SortedSet{Variable}}` providing the set of variables involved in the SDP localizing matrix corresponding to each constraint.")
    println("-> Nb SDP constraints:                 xx")
    println("-> Nb variables by SDP constraint:     xx / xx (mean/std)")
    cliquevarsbycstr = SortedDict{String, SortedSet{Variable}}()
    return cliquevarsbycstr::SortedDict{String, SortedSet{Variable}}
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
