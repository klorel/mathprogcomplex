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

function get_pbcstrname(cstrname::String)
    if ismatch(r".+(\_lo|\_hi|\_eq)", cstrname)
        return cstrname[1:end-3]
    else
        return cstrname
    end
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