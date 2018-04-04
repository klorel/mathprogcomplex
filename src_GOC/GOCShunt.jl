"""
    struct GOCShunt <: AbstractNodeLabel

Structure descendant of AbstractNodeLabel

# Fields
- `busname::String`
- `id::String`
- `shunt::Complex128`: shunt = Gs + im Bs

"""
struct GOCShunt <: AbstractNodeLabel
    busid::Int64
    id::String #label
    shunt::Complex128 # coeffS = Gs + im*Bs
end

"""
    create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, scenario::String) where T <: GOCShunt

Create shunt variable for shunt `elemid` at `bus` in `bus_vars` if `elem_formulation == :NewVar`. \n
Return nothing
"""
function create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, scenario::String) where T <: GOCShunt
    if elem_formulation == :NewVar
      bus_vars[elemid] = Variable(variable_name("Sshunt", bus, elemid, scenario), Complex)
    end
    return
end


"""
  [cstrname, polynom, lb, ub] = Snodal(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}) where T <: GOCShunt

Return the polynomial power contribution for shunt `elemid` at `busid`. Will be used to construct power balance constraints in polynomial problem.\n
If `elem_formulation == :NbMinVar`, return expression of Sshunt depending on the voltage `[cstrname, conj(shunt_value)|V|^2, 0, 0]`\n
If `elem_formulation == :NewVar`, return variable Sshunt `[cstrname, Polynomial(Shunt), 0, 0]`\n
Return `["", Polynomial(), 0, 0]` otherwise.
"""
function Snodal(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}) where T <: GOCShunt
    cstrname = "SHUNT"
    if elem_formulation == :NbMinVar
      p = conj(element.shunt) * abs2(bus_vars[volt_name()])
      return [cstrname, p, 0, 0]
    elseif elem_formulation == :NewVar
      return [cstrname,Polynomial(bus_vars[elemid]), 0, 0]
    end
    return ["", Polynomial(), 0, 0]

end


"""
    constraint(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T <: GOCShunt

Return a constraint to define shuntvariable for `elemid` at `bus` if `elem_formulation == : NewVar`: `"SHUNT"=> Sshunt-conj(shunt_value)|V|^2==0` \n
Return empty dictionary of constraints otherwise.

"""
function constraint(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T <: GOCShunt
    if elem_formulation == :NewVar
      cstrname = get_ShuntDef_cstrname()
      Sshunt = element.shunt* abs2(bus_vars[volt_name()])
      return Dict{String, Constraint}(cstrname => Polynomial(Sshunt - bus_vars[elemid]) == 0)
    end
   return Dict{String, Constraint}()
end
