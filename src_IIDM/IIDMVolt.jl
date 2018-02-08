## Element type
struct IIDMVolt <: AbstractNodeLabel
  busname::String
  busid::Int
  substationname::String
  nomV::Float64    # kV
  minV::Float64    # p.u.
  maxV::Float64    # p.u.
end

## 1. Variables creation
function create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, scenario::String) where T <: IIDMVolt
  bus_vars[elemid] = Variable(variable_name("VOLT", bus, "", scenario), Complex)
  return
end

# function create_vars!(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::Dict{String, Variable}) where T <: MyType
#   return
# end


## 2. Power balance
# function Snodal(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}) where T <: MyType
#   return ["", Polynomial(), 0, 0]
# end

# function Sorig(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::Dict{String, Variable}) where T<:MyType
#   return Polynomial()
# end

# function Sdest(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::Dict{String, Variable}) where T<:MyType
#   return Polynomial()
# end


## 3. Constraints creation
function constraint(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}) where T <: IIDMVolt
  cstrname = "VOLTM_$busid"
  Vmin = element.minV
  Vmax = element.maxV
  nomV = element.nomV
  p = abs2(bus_vars["Volt"])
  return Dict{String, Constraint}(cstrname=> Vmin^2 << p << Vmax^2)
end

# function constraint(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::Dict{String, Variable}) where T <: MyType
#   return [[], [], [], []]
# end


# ## 4. Measure accessor
# function add_measure(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}) where T <: IIDMVolt
#   return Dict("V_$(busid)" => abs2((1/element.nomV)*bus_vars["Volt"]) ==  element.V)
# end
#
# # function add_measure(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::Dict{String, Variable}) where T
# #   return Dict{String, Constraint}()
# # end
