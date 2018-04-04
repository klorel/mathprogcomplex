"""
    struct GOCLoad <: AbstractNodeLabel

Structure descendant of AbstractNodeLabel

# Fields
- `busname::String`
- `id::String`: necessary if several loads at a bus
- `load::Complex128`: Sd = Pd + im Qd

"""
struct GOCLoad <: AbstractNodeLabel
    busname::Int64
    id::String #label
    load::Complex128 #Sd = Pd+im*Qd
end

"""
    create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, scenario::String) where T <: GOCLoad

Create load variable at `bus` for `elemid` in `bus_vars` if `elem_formulation == : NewVar`. \n
Return nothing
"""
function create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, scenario::String) where T <: GOCLoad
    if elem_formulation == :NewVar
      bus_vars[elemid] = Variable(variable_name("Sload", bus, elemid, scenario), Complex)
    end
    return
end


"""
    [cstrname, polynom, lb, ub] = Snodal(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}) where T <: GOCLoad

Return the polynomial contribution in power of load `elemid` at `busid`. Will be used to construct power balance constraints in polynomial problem.\n
If `elem_formulation == :NbMinVar`, return value of Sload `[cstrname, Sload, 0, 0]`\n
If `elem_formulation == :NewVar`, return variable Sload `[cstrname, Polynomial(Sload), 0, 0]`\n
Return `["", Polynomial(), 0, 0]` otherwise.
"""
function Snodal(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}) where T <: GOCLoad
    cstrname = "LOAD"
    if elem_formulation == :NbMinVar
      Sload = element.load
      return [cstrname, Sload, 0, 0]
    elseif elem_formulation == :NewVar
      return [cstrname, Polynomial(bus_vars[elemid]), 0, 0]
    else
        return ["", Polynomial(), 0, 0]
    end
end


"""
    constraint(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T <: GOCLoad

Return a constraint to define load variable for `elemid` at `bus` if `elem_formulation == : NewVar`: `"LOAD"=> Sload==value` \n
Return empty dictionary of constraints otherwise.

"""
function constraint(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T <: GOCLoad
    if elem_formulation == :NewVar
      cstrname = get_LoadDef_cstrname()
      Sload = element.load
      return Dict{String, Constraint}(cstrname => Polynomial(bus_vars[elemid]) == Sload)
    end
    return Dict{String, Constraint}()
end
