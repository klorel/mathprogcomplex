get_cstrname_lower(cstrname::String) = cstrname*"_lo"
get_cstrname_upper(cstrname::String) = cstrname*"_hi"
get_cstrname_eq(cstrname::String) = cstrname*"_eq"

get_momentcstrname() = "moment_cstr"

function get_blockname(cstrname, cliquename, mmtrelax_pb)
    cstrcliques = Set(Iterators.filter(x->x[1] == cstrname, keys(mmtrelax_pb.constraints)))
    if length(cstrcliques) == 1
        return cstrname
    else
        cstrname == get_momentcstrname() || warn("get_blockname(): several cliques found for $cstrname : $cstrcliques")
        return cstrname*"_"*cliquename
    end
end

function get_cstrname(cstrname::String, cstrtype::Symbol)
    if cstrtype == :ineqdouble
        return get_cstrname_lower(cstrname), get_cstrname_upper(cstrname)
    elseif cstrtype == :eq
        return get_cstrname_eq(cstrname)
    elseif cstrtype == :ineqlo
        return cstrname
        # return get_cstrname_lower(cstrname)
    elseif cstrtype == :ineqhi
        return cstrname
        # return get_cstrname_upper(cstrname)
    else
        error("get_cstrname(): Unknown type $cstrtype.")
    end
end

###############################################################################
## Moment Relaxation problem

function get_normalizedpoly(cstr::Constraint, cstrtype::Symbol)
    if cstrtype == :eq
        return cstr.p - cstr.lb
    elseif cstrtype == :ineqlo
        return cstr.p - cstr.lb
    elseif cstrtype == :ineqhi
        return cstr.ub - cstr.p
    else
        error("get_normalizedpoly(): Unhandeld type $cstrtype.")
    end
end

"""
    ginorm = set_givars(gi, vars, vars_overlap)

    replacing variables occuring in several cliques by `vars` variables.
    NOTE: this *should* be done inplace for efficiency.
"""
function set_givars(gi::Polynomial, authorized_vars::SortedSet{Variable}, vars_overlap::SortedDict{Variable, SortedSet{String}}, ctr_cliques::SortedSet{String})
    gi_res = Polynomial()
    # info("poly input :")
    # println(gi)
    # info("Authorized vars :")
    # for var in authorized_vars print("$var  ") end
    # println()
    # info("vars_overlap:")
    # for (var, cliques) in vars_overlap
        # print("$var  ⇥ ")
        # for clique in cliques print("$clique ") end
        # println()
    # end

    for (expo, λ) in gi
        # warn(" - $expo")
        cur_expo = Exponent()
        for (var, degree) in expo
            # warn("   |-> $var")
            cur_var = var
            if !(var ∈ authorized_vars)
                @assert haskey(vars_overlap, var)
                inter = intersect(ctr_cliques, vars_overlap[var])
                # info("         cliques of var $var : $(vars_overlap[var])")
                # info("         cliques of var $var and ctr : $(inter)")
                clique = first(inter)
                cur_var = get_varinclique(var, clique)
            end
            cur_var in authorized_vars || println("((( $var  not $cur_var")
            @assert cur_var in authorized_vars
            add_expod!(cur_expo, Exponent(cur_var=>degree))
        end
        # warn("   |==> $cur_expo")
        add!(gi_res, cur_expo * λ)
    end
    # info("           : $gi_res")
    # info("          -> $(gi == gi_res)")
    return gi_res
end



function get_pbcstrname(cstrname::String)
    if ismatch(r".+(\_lo|\_hi|\_eq)", cstrname)
        return cstrname[1:end-3]
    else
        return cstrname
    end
end

function get_ccmultvar(var::Variable, clique1::String, clique2::String)
    return Variable("lagmult_cc_"*var.name*"_"*clique1*"_"*clique2, var.kind)
end

function get_varinclique(var::Variable, clique::String)
    return Variable(var.name*"_"*clique, var.kind)
end

function format_string(α::Exponent, β::Exponent)
    s = "$α,$β"
    return replace(s, " ", "_")
end

function format_string(s1::String, s2::String)
    s = "$s1,$s2"
    return replace(s, " ", "_")
end

function format_string(α::Exponent)
    s = "$α"
    return replace(s, " ", "_")
end

function change_eq_to_ineq!(problem::Problem)
    for (ctrname, ctr) in problem.constraints
        if get_cstrtype(ctr) == :eq
            rm_constraint!(problem, ctrname)
            add_constraint!(problem, get_cstrname_lower(ctrname), 0 << (ctr.p - ctr.lb))
            add_constraint!(problem, get_cstrname_upper(ctrname), 0 << (ctr.ub - ctr.p))
        end
    end
end

## Mosek
obj_key() = "1,1"