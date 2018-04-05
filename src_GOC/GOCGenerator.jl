"""
    mutable struct GOCGenerator <: AbstractNodeLabel

Mutable structure descendant of AbstractNodeLabel

# Fields
- `busname::String`: identify the bus associated to the generator
- `id::String`: identify the generator for the bus `busname`
- `power_min::Complex128`: Smin = Pmin + im Qmin
- `power_max::Complex128`: Smax = Pmax + im Qmax
- `participation_factor::Float64`: coefficient to define active power coupling constraints
- `dict_obj_coeffs::Dict{Int64,Float64}`: dictionary for objective coefficients: degree => coeff
"""
mutable struct GOCGenerator <: AbstractNodeLabel
  busid::Int64
  id::Any
  power_min::Complex128 #Smin = Pmin + i Qmin
  power_max::Complex128 #Smax = Pmax + i Qmax
  participation_factor::Float64
  dict_obj_coeffs::Dict{Int64,Float64} # dict : degree => coeff
end



"""
    create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, scenario::String) where T <: GOCGenerator

Create variables for generator `elemid` at `bus`  in `bus_vars` if necessary.\n
If `elem_formulation == :NewVar`, create Sgen variable for the generator\n
If `elem_formulation == :GOCcoupling`, create Sgen variable for the generator, binary variables for comparison of voltages and create delta variable\n
Return nothing
"""
function create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, scenario::String) where T <: GOCGenerator
    if elem_formulation == :NewVar
        bus_vars[elemid] = Variable(variable_name("Sgen", bus, elemid, scenario), Complex)
    elseif elem_formulation == :GOCcoupling
        bus_vars[elemid] = Variable(variable_name("Sgen", bus, elemid, scenario), Complex)

        # Binary variables for coupling constraints
        Vbc_inf_Vsc = get_binInf_varname(basecase_scenario_name(), scenario, bus)
        bus_vars[Vbc_inf_Vsc] = Variable(Vbc_inf_Vsc, Bool)

        Vsc_inf_Vbc = get_binInf_varname(scenario, basecase_scenario_name(), bus)
        bus_vars[Vsc_inf_Vbc] = Variable(Vsc_inf_Vbc, Bool)

        Veq = get_binEq_varname(scenario, basecase_scenario_name(), bus)
        bus_vars[Veq] = Variable(Veq, Bool)

        # Scaling factor for current scenario
        delta_varname = get_delta_varname(scenario)
        bus_vars[delta_varname] = Variable(delta_varname, Real)
    end

    return
end

"""
    [cstrname, polynom, lb, ub] = Snodal(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}) where T <: GOCGenerator

Return the polynomial contribution in power of generator `elemid` at `bus`(name, value, lower bound, upper bound). Will be used to construct power balance constraints in polynomial problem.\n
If `elem_formulation == :NbMinVar`, return generator bounds `["UNIT", Polynomial() ,Smin, Smax]`\n
If `elem_formulation == :NewVar` or `:GOCCoupling`, return -Sgen `["UNIT", -Sgen, 0, 0]`\n
Return no contribution `["", Polynomial(), 0, 0]` otherwise
"""
function Snodal(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}) where T <: GOCGenerator
    cstrname = "UNIT"
    if elem_formulation == :NbMinVar
      lb = element.power_min
      ub = element.power_max
      return [cstrname, Polynomial(), lb, ub]
  elseif elem_formulation == :NewVar || elem_formulation == :GOCcoupling
      return [cstrname, Polynomial(-bus_vars[elemid]), 0, 0]
    end
    return ["", Polynomial(), 0, 0]
end

"""
    constraint(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T <: GOCGenerator

Return all the constraints defined by generator `elemid` at `bus`. Will be used to construct constraints in polynomial problem.\n
If `elem_formulation == :NewVar`, return generator bounds : "Genbounds" => Smin <= Sgen <= Smax\n
If `elem_formulation == :GOCCoupling`, return generator bounds defined by coupling constraints\n
Return empty dictionary otherwise
"""
function constraint(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T <: GOCGenerator
    cstrs = Dict{String, Constraint}()
    if elem_formulation == :NewVar
      lb = element.power_min
      ub = element.power_max
      cstrname = get_GenBounds_cstrname()
      cstrs[cstrname] = lb << bus_vars[elemid] << ub
    end

    ## Active / Reactive formulation
    if elem_formulation == :GOCcoupling
        # Problem variables and parameters
        Sgen = bus_vars[elemid]

        Volt_sc = bus_vars[volt_name()]
        Volt_bc = OPFpbs[basecase_scenario_name()].mp.node_vars[bus][volt_name()]

        Vbc_inf_Vsc = bus_vars[get_binInf_varname(basecase_scenario_name(), scenario, bus)]
        Vsc_inf_Vbc = bus_vars[get_binInf_varname(scenario, basecase_scenario_name(), bus)]
        Veq = bus_vars[get_binEq_varname(scenario, basecase_scenario_name(), bus)]

        Pmin, Qmin = real(element.power_min), imag(element.power_min)
        Pmax, Qmax = real(element.power_max), imag(element.power_max)
        ϵ, M = get_GOC_Volt_ϵ(), get_GOC_BigM()

        ## Definition of binary variables, upper constraint
        ccname_bindef_upper = get_VoltBinDef_upper()
        cstrs[ccname_bindef_upper] = (abs2(Volt_bc) - abs2(Volt_sc) - (Vbc_inf_Vsc * (-ϵ) + (1-Vbc_inf_Vsc+Veq) * M - 2*M*Veq)) << 0
        ## Definition of binary variables, lower constraint
        ccname_bindef_lower = get_VoltBinDef_lower()
        ##TODO: verify big M constraints
        cstrs[ccname_bindef_lower] = (abs2(Volt_bc) - abs2(Volt_sc) - (Vsc_inf_Vbc * ϵ + (1-Vsc_inf_Vbc+Veq) * (-M) + 2*M*Veq)) >> 0
        ####
        # If Vsc < Vbc, lhs = abs2(Volt_bc) - abs2(Volt_sc) > 0
        ## If Vbc_inf_Vsc = 1, positive lhs < -ϵ => upper constraint not satisfied => NOT POSSIBLE
        ## If Veq = 1, positive lhs <  0 => upper constraint not satisfied => NOT POSSIBLE
        ## If Vsc_inf_Vbc = 1, positive lhs < M, upper constraint satisfied AND lhs > ϵ, lower constraint satisfied

        # If Vbc < Vsc, lhs = abs2(Volt_bc) - abs2(Volt_sc) < 0
        ## If Vbc_inf_Vsc = 1, lhs < -ϵ, upper constraint satisfied AND lhs > -M, lower constraint satisfied
        ## If Veq = 1,  lhs < 0, upper constraint satisfied AND lhs > 0, lower constraint not satisfied => NOT POSSIBLE
        ## If Vsc_inf_Vbc = 1, lhs < M, upper constraint satisfied AND lhs > ϵ, lower constraint not satisfied => NOT POSSIBLE

        # If |Vbc-Vsc|< ϵ, -ϵ < abs2(Volt_bc) - abs2(Volt_sc) < ϵ
        ## If Vbc_inf_Vsc = 1, lhs < -ϵ, upper constraint not satisfied => NOT POSSIBLE
        ## If Veq = 1,  lhs < 0, upper constraint satisfied AND lhs > 0, lower constraint satisfied
        ## If Vsc_inf_Vbc = 1, lhs < M, upper constraint satisfied AND lhs > ϵ, lower constraint not satisfied => NOT POSSIBLE

        ####
        ## Definition of binary variables, complementarity
        ccname_compl = get_VoltBinDef_complement()
        ## Three separate possibilities => only one binary variable == 1
        cstrs[ccname_compl] = (Vbc_inf_Vsc + Vsc_inf_Vbc + Veq == 1)

        P0 = real(OPFpbs[basecase_scenario_name()].mp.node_vars[bus][elemid])
        Pgen, Qgen = real(Sgen), imag(Sgen)
        deltavar = bus_vars[get_delta_varname(scenario)]
        α = element.participation_factor

        ## Coupling constraint active :
        cstrname = get_CC_active_cstrname()
        cstrs[cstrname] = (Pgen - (P0 + α * deltavar)) == 0

        ## Coupling constraint reactive : upper bound
        ccname_upper = get_CC_reactiveupper_cstrname()
        cstrs[ccname_upper] = ((Qgen - Vbc_inf_Vsc * (Qmin - Qmax) - Qmax)) << 0

        ## Coupling constraint reactive : lower bound
        ccname_lower = get_CC_reactivelower_cstrname()
        cstrs[ccname_lower] = ((Qgen - Vsc_inf_Vbc * (Qmax - Qmin) - Qmin)) >> 0
    end

    return cstrs
end


"""
    cost(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, Snode::Polynomial, lb, ub) where T <: GOCGenerator

Return the polynomial contribution in objective generator `elemid` at `bus`.
"""
function cost(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::Dict{String, Variable}, Snode::Polynomial, lb, ub) where T <: GOCGenerator
  Sgen = gencost = Polynomial()
  if elem_formulation == :NbMinVar
    _, S, lb_gen, ub_gen = Snodal(element, bus, elemid, elem_formulation, bus_vars)
    abs((lb-lb_gen) - (ub-ub_gen)) ≤ 5e-5 || warn("Power missing at $bus... $(abs((lb-lb_gen) - (ub-ub_gen)))")
    Sgen = Snode - S - (ub-ub_gen)
    # Sgen = add(Snode, -S - (ub-ub_gen))
  elseif elem_formulation == :NewVar
    Sgen::Polynomial = bus_vars[elemid]
    # Sgen = Polynomial(bus_vars[elemid])
  end
  add!(Sgen, conj(Sgen))
  Sactive = (Sgen)*0.5
  Sexp = 1 + Polynomial()
  dict_coeffs = element.dict_obj_coeffs
  degrees_sorted = sort(collect(keys(dict_coeffs)))
  #quartic objective
  imax = 3
  #quadratic objective :TODO
  #imax = 2
  if imax != 3
    warning("Objective is troncated to $imax degree")
  end
  for i=1:imax
        degree = degrees_sorted[i] ## NOTE: NOT normal behaviour, required for linear cost...
        add!(gencost, dict_coeffs[degree]*Sexp)
        Sexp = Sexp*Sactive
  end
  return gencost
end
