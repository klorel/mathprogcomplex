"""
    check_feasibility(elem::T, bus::String, elemid::String, elem_formulation::Symbol, scenario::String, point::Point) where T<:AbstractNodeLabel

Evaluate the constraints defined by `elemid` of type `elem` at `bus` for `scenario` (depending on `elem_formulation`) on `point` and return a dictionary of infeasible constraints SortedDict(ctrname => message with information)
Return nothing if there is not violated constraints.

"""
function check_feasibility(elem::T, bus::String, elemid::String, elem_formulation::Symbol, scenario::String, point::Point,epsilon::Float64) where T<:AbstractNodeLabel
    return
end


"""
    check_feasibility(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}, scenario::String, point::Point, epsilon::Float64) where T<:AbstractLinkLabel

Evaluate the constraints defined by `elemid` of type `elem` at `link` for `scenario` (depending on `elem_formulation`) on `point` and return a dictionary of infeasible constraints SortedDict(ctrname => message with information)
Return nothing if there is not violated constraints.

"""

function check_feasibility(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}, scenario::String, point::Point, epsilon::Float64) where T<:AbstractLinkLabel
    return
end
