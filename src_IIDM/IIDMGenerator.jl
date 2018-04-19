## Element type
struct IIDMGenerator <: AbstractNodeLabel
  id::String
  minP::Float64 # MW
  maxP::Float64 # MW
  targetV::Float64 # p.u.
  volt_regulator_on::Bool
  S_bounds::Array{Tuple{Complex, Complex}} # MW+im*MVar
end

## 1. Variables creation
function create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}, scenario::String) where T <: IIDMGenerator
  bus_vars[elemid] = Variable(variable_name("Sgen", bus, elemid, scenario), Complex)
  return
end

# function create_vars!(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T <: MyType
#   return
# end

## 2. Power balance
function Snodal(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}) where T <: IIDMGenerator
  cstrname = "UNIT"
  return [cstrname, Polynomial(bus_vars[elemid]), 0, 0]
end

# function Sorig(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T<:MyType
#   return Polynomial()
# end

# function Sdest(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T<:MyType
#   return Polynomial()
# end

## 3. Constraints creation
function constraint(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}) where T <: IIDMGenerator
  bounds = sort(element.S_bounds, by=x->real(x[1]))
  length(bounds) == 2 || warn("$(length(bounds))x2 ≠ 4 points for power domain of $elemid, $busid")

  Sgen = bus_vars[elemid]
  constraints = SortedDict{String, Constraint}()

  constraints["Cstr_$(busid)_$(elemid)_real"] = real(bounds[1][1]) << Sgen << real(bounds[2][1])

  a, θ = bounds[1][2], angle(bounds[2][2] - bounds[1][2])
  p =  real((Sgen-a)*exp(im*(π/2-θ)))
  constraints["Cstr_$(busid)_$(elemid)_upper"] = 0 << p << Inf

  a, θ = bounds[1][1], angle(bounds[2][1] - bounds[1][1])
  p = real((Sgen-a)*exp(im*(-π/2-θ)))
  constraints["Cstr_$(busid)_$(elemid)_lower"] = 0 << p << Inf

  return constraints
end

# function constraint(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T <: MyType
#   return SortedDict{String, Constraint}()
# end


## 4. Generator cost
function cost(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_elems_var::SortedDict{String, Variable}, Snode::Polynomial, lb, ub) where T <: IIDMGenerator
  Sgen = gencost = Polynomial()
  if elem_formulation == :NbMinVar
    _, S, lb_gen, ub_gen = Snodal(element, bus, elemid, elem_formulation, bus_elems_var)
    abs((lb-lb_gen) - (ub-ub_gen)) ≤ 5e-5 || warn("Power missing at $bus... $(abs((lb-lb_gen) - (ub-ub_gen)))")
    Sgen = add(Snode, -S - (ub-ub_gen))
  elseif elem_formulation == :NewVar
    Sgen = Polynomial(bus_elems_var[elemid])
  else
    warn("No cost contribution for gen $elemid, $bus")
  end

  add!(Sgen, conj(Sgen))
  gencost = (Sgen)*0.5
  return gencost
end
