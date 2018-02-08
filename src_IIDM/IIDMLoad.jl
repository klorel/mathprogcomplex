## Element types
struct IIDMLoad <: AbstractNodeLabel
  id::String
  load::Complex128 # MW + im*MVar
end


## 1. Variables creation
function create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, scenario::String) where T <: IIDMLoad
  if elem_formulation == :NewVar
    bus_vars[elemid] = Variable(variable_name("Sload", bus, elemid, scenario), Complex)
  end
  return
end

# function create_vars!(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::Dict{String, Variable}) where T <: MyType
#   return
# end


## 2. Power balance
function Snodal(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}) where T <: IIDMLoad
  cstrname = "LOAD"
  if elem_formulation == :NbMinVar
    Sload = element.load
    return [cstrname, Polynomial(), -Sload, -Sload]
  elseif elem_formulation == :NewVar
    return [cstrname, Polynomial(bus_vars[elemid]), 0, 0]
  end
  return ["", Polynomial(), 0, 0]
end

# function Sorig(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::Dict{String, Variable}) where T<:MyType
#   return Polynomial()
# end

# function Sdest(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::Dict{String, Variable}) where T<:MyType
#   return Polynomial()
# end


## 3. Constraints creation
function constraint(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}) where T <: IIDMLoad
  if elem_formulation == :NewVar
    cstrname = "Cstr_$(busid)_$(elemid)"
    Sload = element.load
    p = Polynomial(bus_vars[elemid])
    return Dict{String, Constraint}(cstrname=> p == Sload)
  end
  return Dict{String, Constraint}()
end

# function constraint(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::Dict{String, Variable}) where T <: MyType
#   return [[], [], [], []]
# end
