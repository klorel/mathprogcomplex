"""
    struct MatpowerVolt <: AbstractNodeLabel

Structure descendant of AbstractNodeLabel

# Fields
- `busname::String`
- `busid::Int64`
- `voltage_magnitude_min::Float64`: Vmin
- `voltage_magnitude_max::Float64`: Vmax
"""
struct MatpowerVolt <: AbstractNodeLabel
  busname::String
  busid::Int64
  voltage_magnitude_min::Float64 #Vmin
  voltage_magnitude_max::Float64 #Vmax
end


"""
    create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}, scenario::String) where T <: MatpowerVolt

Create voltage variable for `bus` in `bus_vars`.\n
Return nothing.
"""
function create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}, scenario::String) where T <: MatpowerVolt
  bus_vars[elemid] = Variable(variable_name("VOLT", bus, "", scenario), Complex)
  return
end


## 2. Power balance
# function Snodal(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}) where T <: MyType
#   return ["", Polynomial(), 0, 0]
# end


"""
    constraint(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T <: MatpowerVolt

Return the constraint of bounds on voltage magnitude at `busid`: "VOLTM" => Vmin^2 <= |V|^2 <= Vmax^2

"""
function constraint(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T <: MatpowerVolt
  cstrname = get_VoltM_cstrname()
  Vmin = element.voltage_magnitude_min
  Vmax = element.voltage_magnitude_max
  return SortedDict{String, Constraint}(cstrname => Vmin^2 << abs2(bus_vars[volt_name()]) <<Vmax^2)
end
