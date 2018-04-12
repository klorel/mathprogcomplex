## Element type
struct IIDMLine_π <: AbstractLinkLabel
  resistance::Float64 # pu
  reactance::Float64 # pu
  susceptance_in::Complex # pu
  susceptance_out::Complex # pu
  basekV_orig::Float64 # kV
  basekV_dest::Float64 # kV
  permanent_limit::Float64
end


## 1. Variables creation
# function create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}) where T <: MyType
#   return
# end

function create_vars!(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}, scenario::String) where T <: IIDMLine_π
  orig, dest = link.orig, link.dest
  link_vars["Volt_orig"] = Variable(variable_name("VOLT", orig, "", scenario), Complex)
  link_vars["Volt_dest"] = Variable(variable_name("VOLT", dest, "", scenario), Complex)
  return
end


## 2. Power balance
# function Snodal(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}) where T <: MyType
#   return ["", Polynomial(), 0, 0]
# end

function Sorig(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T <: IIDMLine_π
  return link_vars["Volt_orig"] * conj(get_IorigMP(element, link_vars))
end

function Sdest(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T <: IIDMLine_π
  return link_vars["Volt_dest"] * conj(get_IdestMP(element, link_vars))
end


## 3. Constraints creation
# function constraint(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}) where T <: MyType
#   return [[], [], [], []]
# end

# function constraint(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T <: MyType
#   return [[], [], [], []]
# end


## 4. Measures
# function add_measure(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}) where T
#   return SortedDict{String, Constraint}()
# end

function add_measure(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T <: IIDMLine_π
  measure_equations = SortedDict{String, Constraint}()
  measure_equations["Sorig_$(link.orig)→$(link.dest)"] = Sorig(element, link, elemid, elem_formulation, link_vars) == element.S_orig
  measure_equations["Sdest_$(link.orig)→$(link.dest)"] = Sdest(element, link, elemid, elem_formulation, link_vars) == element.S_dest
  return measure_equations
end



## Utils functions
function get_IorigMP(element::IIDMLine_π, link_variables::SortedDict{String, Variable})
  rs = element.resistance
  xs = element.reactance
  bc_in = element.susceptance_in
  ys = 1/(rs+im*xs)
  return (ys + bc_in)*link_variables["Volt_orig"] - ys*link_variables["Volt_dest"]
end

function get_IdestMP(element::IIDMLine_π, link_variables::SortedDict{String, Variable})
  rs = element.resistance
  xs = element.reactance
  bc_out = element.susceptance_out
  ys = 1/(rs+im*xs)
  return -ys*link_variables["Volt_orig"] + (ys + bc_out)*link_variables["Volt_dest"]
end
