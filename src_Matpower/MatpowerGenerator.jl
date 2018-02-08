"""
    mutable struct MatpowerGenerator <: AbstractNodeLabel

Mutable structure descendant of AbstractNodeLabel

# Fields
- `power_min::Complex128`
- `power_max::Complex128`
- `cost_degree::Int64`
- `cost_coeffs::Array{Float64}`
- `gen_status::Bool`: generator on/off

"""
mutable struct MatpowerGenerator <: AbstractNodeLabel
  power_min::Complex128 #Smin = Pmin + i Qmin
  power_max::Complex128 #Smax = Pmax + i Qmax
  cost_degree::Int64 # degree of the cost polynomial
  cost_coeffs::Array{Float64} # cost polynomial coeffs from higher to lower degree
  gen_status::Bool
end


"""
    create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, scenario::String) where T <: MatpowerGenerator

Create power variables n for generator `elemid` at `bus` in `bus_vars` if `elem_formulation == :NewVar`\n
Return nothing

"""
function create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, scenario::String) where T <: MatpowerGenerator
  if elem_formulation == :NewVar
    bus_vars[elemid] = Variable(variable_name("Sgen", bus, elemid, scenario), Complex)
  end
  return
end

"""
[cstrname, polynom, lb, ub] = Snodal(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}) where T <: MatpowerGenerator

Return the polynomial contribution in power of generator `elemid` at `bus`(name, value, lower bound, upper bound). Will be used to construct power balance constraints in polynomial problem.\n
If `elem_formulation == :NbMinVar`, return generator bounds `["UNIT", Polynomial() ,Smin, Smax]`\n
If `elem_formulation == :NewVar`, return -Sgen `["UNIT", -Sgen, 0, 0]`\n
Return no contribution `["", Polynomial(), 0, 0]` otherwise
"""
function Snodal(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}) where T <: MatpowerGenerator
  cstrname = "UNIT"
  if elem_formulation == :NbMinVar
    lb = element.power_min
    ub = element.power_max
    if !element.gen_status ## generator is off
      lb = ub = 0
    end
    return [cstrname, Polynomial(), lb, ub]
  elseif elem_formulation == :NewVar
    return [cstrname, Polynomial(-bus_vars[elemid]), 0, 0]
  end
  return ["", Polynomial(), 0, 0] # Raise error here ? Keep doing nothing ? One would expect a generator exists if this is called.
end


"""
    constraint(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T <: MatpowerGenerator

Return all the constraints defined by generator `elemid` at `bus`. Will be used to construct constraints in polynomial problem.\n
If `elem_formulation == :NewVar`, return generator bounds : "Genbounds" => Smin <= Sgen <= Smax\n
Return empty dictionary otherwise
"""
function constraint(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T <: MatpowerGenerator
  if elem_formulation == :NewVar
    lb = element.power_min
    ub = element.power_max
    if !element.gen_status
      ub = lb = 0
    end
    cstrname = get_GenBounds_cstrname()
    p = Polynomial(-bus_vars[elemid])
    return Dict{String, Constraint}(cstrname => lb << p << ub)
  end
  return Dict{String, Constraint}()
end


"""
    cost(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_elems_vars::Dict{String, Variable}, Snode::Polynomial, lb, ub) where T <: MatpowerGenerator

Return the polynomial contribution in objective generator `elemid` at `bus`.
"""
function cost(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_elems_var::Dict{String, Variable}, Snode::Polynomial, lb, ub) where T <: MatpowerGenerator
  Sgen = Polynomial()
  gencost = Polynomial()
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
  Sactive = (Sgen)*0.5
  Sexp = 1 + Polynomial()
  d = element.cost_degree
  degrees = element.cost_coeffs
  for i=1:2                                               ## NOTE: NOT normal behaviour, required for linear cost...
    add!(gencost, degrees[d-i+1] * Sexp)
    Sexp = Sexp*Sactive
  end
  if !element.gen_status
    gencost = 0
  end
  return gencost
end
