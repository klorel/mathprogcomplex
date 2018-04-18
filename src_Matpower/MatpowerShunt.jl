"""
    struct MatpowerShunt <: AbstractNodeLabel

Structure descendant of AbstractNodeLabel

# Fields
- `Yshunt::Complex128`: Yshunt = Gs + im Bs
"""
struct MatpowerShunt <: AbstractNodeLabel
  Yshunt::Complex128 # coeffS = Gs + im*Bs
end


## 1. Variables creation
# function create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}) where T <: MyType
#   return
# end


"""
  [cstrname, polynom, lb, ub] = Snodal(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}) where T <: MatpowerShunt

Return the polynomial power contribution for shunt `elemid` at `busid`. Will be used to construct power balance constraints in polynomial problem.\n
Return expression of Sshunt depending on the voltage `["SHUNT", conj(shunt_value)|V|^2, 0, 0]`\n

"""
function Snodal(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}) where T <: MatpowerShunt
  p = conj(element.Yshunt) * abs2(bus_vars[volt_name()])
  return ["SHUNT", p, 0, 0]
end


## 3. Constraints creation
# function constraint(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T <: MyType
  # return SortedDict{String, Constraint}()
# end
