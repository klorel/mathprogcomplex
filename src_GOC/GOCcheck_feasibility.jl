"""
    check_feasibility(elem::T, bus::String, elemid::String, elem_formulation::Symbol, scenario::String, point::Point) where T<:AbstractNodeLabel

Evaluate the constraints defined by `elemid` of type `elem` at `bus` for `scenario` (depending on `elem_formulation`) on `point` and return a dictionary of infeasible constraints SortedDict(ctrname => message with information)
Return nothing if there is not violated constraints.

"""


###Genbounds constraints
function check_feasibility(element::T, bus::String, elemid::String, elem_formulation::Symbol, scenario::String, point::Point,epsilon::Float64) where T<:GOCGenerator
    if elem_formulation == :NewVar
      lb = element.power_min
      ub = element.power_max
      cstrname = get_cstrname(scenario, bus, elemid, "Genbounds")
      varname = variable_name("Sgen", bus, elemid, scenario)
      value = point[Variable(varname,Complex)]
      #println(lb, " <= ", value, " <= ", ub)
      slack_left = value - lb
      slack_right = ub - value
      min_slack = min(real(slack_left),real(slack_right),imag(slack_left),imag(slack_right))
      if min_slack < -epsilon
          #infeasibility
          message = "Power generated not in diagram PQ."
          if real(slack_left) < - epsilon
              message = message*" Pgen < Pmin"
          elseif imag(slack_left) < - epsilon
              message = message*" Qgen < Qmin"
          elseif real(slack_right) < - epsilon
              message = message*" Pgen > Pmax"
          elseif imag(slack_left) < - epsilon
              message = message*" Qgen > Qmax"
          end
          return SortedDict(cstrname => message)
      end
    end
end

# ##test
# scenario = "BaseCase"
# bus = "BUS_1"
# ds_bus = OPFpbs[scenario].ds.bus[bus]
# elemid = "Generator_1"
# elem = ds_bus[elemid]
# elem_formulation = OPFpbs[scenario].mp.node_formulations[bus][elemid]
# point = global_point
# point_infeasible = copy(global_point)
# point_infeasible[Variable(variable_name("Sgen", bus, elemid, scenario), Complex)] += 1000
#
# check_feasibility(elem, bus, elemid, elem_formulation, scenario, point,1e-6)
# check_feasibility(elem, bus, elemid, elem_formulation, scenario, point_infeasible,1e-6)

###VOLTM constraints
function check_feasibility(element::T, bus::String, elemid::String, elem_formulation::Symbol, scenario::String, point::Point,epsilon::Float64) where T<:GOCVolt
    cstrname = get_cstrname(scenario, bus, elemid, "VOLTM")
    Vmin = element.voltage_magnitude_min
    Vmax = element.voltage_magnitude_max
    varname = variable_name("VOLT", bus, "", scenario)
    value = point[Variable(varname,Complex)]
    module_value = value * conj(value)
    if abs(imag(module_value)) > 1e-15
        warn("Imaginary part of $(varname) module is not zero")
    end
    #println(Vmin^2, " <= ", module_value, " <= ", Vmax^2)
    slack_left = module_value - Vmin^2
    slack_right = Vmax^2 - module_value
    if real(slack_left) > - epsilon
        if real(slack_right) > - epsilon
            ##feasibility at epsilon OK
            return
        else
            return SortedDict(cstrname => "Voltage magnitude above Vmax")
        end
    else
        return SortedDict( cstrname => "Voltage magnitude under Vmin")
    end
end

##test
# scenario = "BaseCase"
# bus = "BUS_1"
# ds_bus = OPFpbs[scenario].ds.bus[bus]
# elemid = "VoltGOC"
# elem = ds_bus[elemid]
# elem_formulation = OPFpbs[scenario].mp.node_formulations[bus][elemid]
# point = global_point
#
# check_feasibility(elem, bus, elemid, elem_formulation, scenario, point, 1e-6)


# function check_feasibility(elem::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}, scenario::String, point::Point) where T<:AbstractLinkLabel
#     return
# end
#
"""
    check_feasibility(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}, scenario::String, point::Point, epsilon::Float64) where T<:AbstractLinkLabel

Evaluate the constraints defined by `elemid` of type `elem` at `link` for `scenario` (depending on `elem_formulation`) on `point` and return a dictionary of infeasible constraints SortedDict(ctrname => message with information)
Return nothing if there is not violated constraints.

"""

##Smax constraints TODO: define for all types of abstract link ???
function check_feasibility(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}, scenario::String, point::Point, epsilon::Float64) where T<:GOCLineπ_notransformer
    origid, destid = link.orig, link.dest
    ctr_name_orig = get_cstrname(scenario, link, elemid, "Smaxorig")
    ctr_name_dest = get_cstrname(scenario, link, elemid, "Smaxdest")
    Sor = Sorig(element, link, elemid, elem_formulation, link_vars)
    Sde = Sdest(element, link, elemid, elem_formulation, link_vars)
    Smax = element.power_magnitude_max
    slack_Sorig = Smax^2 - evaluate(abs2(Sor), point)
    slack_Sdest = Smax^2 - evaluate(abs2(Sde), point)
    if abs(imag(slack_Sorig)) > 1e-15 || abs(imag(slack_Sdest)) > 1e-15
        warn("Imaginary part of Sorig or Sdest module is not zero")
    end
    ctr_slacks = SortedDict(ctr_name_orig => real(slack_Sorig), ctr_name_dest => real(slack_Sdest))
    min_violation_ctr = minimum(ctr_slacks)
    if min_violation_ctr[2] > - epsilon
        return
    else
        if real(slack_Sorig) < - epsilon
            if real(slack_Sdest) < - epsilon
                return SortedDict(ctr_name_orig => "Power magnitude at origin above Smax", ctr_name_dest => "Power magnitude at destination")
            else
                return SortedDict(ctr_name_orig => "Power magnitude at origin is above Smax")
            end
        else
            return SortedDict(ctr_name_dest => "Power magnitude at destination Smax")
        end
    end
end

function check_feasibility(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}, scenario::String, point::Point, epsilon::Float64) where T<:GOCLineπ_withtransformer
    origid, destid = link.orig, link.dest
    ctr_name_orig = get_cstrname(scenario, link, elemid, "Smaxorig")
    ctr_name_dest = get_cstrname(scenario, link, elemid, "Smaxdest")
    Sor = Sorig(element, link, elemid, elem_formulation, link_vars)
    Sde = Sdest(element, link, elemid, elem_formulation, link_vars)
    Smax = element.power_magnitude_max
    slack_Sorig = Smax^2 - evaluate(abs2(Sor), point)
    slack_Sdest = Smax^2 - evaluate(abs2(Sde), point)
    if abs(imag(slack_Sorig)) > 1e-15 || abs(imag(slack_Sdest)) > 1e-15
        warn("Imaginary part of Sorig or Sdest module is not zero")
    end
    ctr_slacks = SortedDict(ctr_name_orig => real(slack_Sorig), ctr_name_dest => real(slack_Sdest))
    min_violation_ctr = minimum(ctr_slacks)
    if min_violation_ctr[2] > - epsilon
        return
    else
        if real(slack_Sorig) < - epsilon
            if real(slack_Sdest) < - epsilon
                return SortedDict(ctr_name_orig => "Power magnitude at origin above Smax", ctr_name_dest => "Power magnitude at destination")
            else
                return SortedDict(ctr_name_orig => "Power magnitude at origin is above Smax")
            end
        else
            return SortedDict(ctr_name_dest => "Power magnitude at destination Smax")
        end
    end
end

#test
# scenario = "BaseCase"
# link = collect(keys(OPFpbs[scenario].ds.link))[1]
# ds_link = OPFpbs[scenario].ds.link[link]
# elemid = linepi_notransformer_name("BL")
# elem = ds_link[elemid]
# link_vars = OPFpbs[scenario].mp.link_vars[link]
# elem_formulation = OPFpbs[scenario].mp.link_formulations[link][elemid]
# point = global_point
# check_feasibility(elem, link, elemid, elem_formulation, link_vars, scenario, point, 1e-6)

"""
    check_feasibility_cc_active_power(elem::T, bus::String, elemid::String, elem_formulation::Symbol, scenario::String, point::Point,epsilon::Float64) where T<:AbstractNodeLabel

Evaluate the coupling constraints about active power defined by `elemid` of type `elem` at `bus` for `scenario` (depending on `elem_formulation`) on `point` and return a dictionary of infeasible constraints SortedDict(ctrname => message with information)
Return nothing if the elem is not a generator or if there is not violated constraints.

"""
function check_feasibility_cc_active_power(elem::T, bus::String, elemid::String, elem_formulation::Symbol, scenario::String, point::Point,epsilon::Float64) where T<:AbstractNodeLabel
    return
end

function check_feasibility_cc_active_power(elem::T, bus::String, elemid::String, elem_formulation::Symbol, scenario::String, point::Point,epsilon::Float64) where T<:GOCGenerator
    if !is_basecase_scenario(scenario)
        cc_name = cc_active_power(scenario, bus, elemid)
        part_factor = elem.participation_factor
        delta_name = delta_varname(scenario)
        sgen_bc_varname = variable_name("Sgen", bus, elemid, basecase_scenario_name())
        sgen_scen_varname = variable_name("Sgen", bus, elemid, scenario)
        delta = point[Variable(delta_name, Real)]
        Sgen_basecase = point[Variable(sgen_bc_varname,Complex)]
        Sgen_scenario = point[Variable(sgen_scen_varname,Complex)]
        diff = real(Sgen_scenario - Sgen_basecase - delta*part_factor)
        if abs(diff) < epsilon
            #feasibility
            return
        else
            return SortedDict(cc_name => "Active power generated in contingency not equal to the one in Basecase plus delta * participation factor")
        end
    else
        #warn("No active power coupling constraints for BaseCase scenario")
        return
    end
end

# scenario = "Scenario1"
# check_feasibility_cc_active_power(elem, bus, elemid, elem_formulation, scenario, point,1e-15)
"""
    check_feasibility_cc_reactive_power(element::T, bus::String, elemid::String, elem_formulation::Symbol, scenario::String, point::Point,epsilon::Float64) where T<:AbstractNodeLabel

Evaluate the coupling constraints about reactive power defined by `elemid` of type `element` at `bus` for `scenario` (depending on `elem_formulation`) on `point` and return a dictionary of infeasible constraints SortedDict(ctrname => message with information)
Return nothing if the elem is not a generator or if there is not violated constraints.

"""
function check_feasibility_cc_reactive_power(element::T, bus::String, elemid::String, elem_formulation::Symbol, scenario::String, point::Point,epsilon::Float64) where T<:AbstractNodeLabel
    return
end

function check_feasibility_cc_reactive_power(element::T, bus::String, elemid::String, elem_formulation::Symbol, scenario::String, point::Point,epsilon::Float64) where T<:GOCGenerator
    if !is_basecase_scenario(scenario)
        cc_name = cc_reactive_power(scenario, bus, elemid)
        sgen_scen_varname = variable_name("Sgen", bus, elemid, scenario)
        volt_bc_name = variable_name("VOLT", bus, "", basecase_scenario_name())
        volt_scen_name = variable_name("VOLT", bus, "", scenario)
        Sgen_scenario = point[Variable(sgen_scen_varname,Complex)]
        V_basecase = point[Variable(volt_bc_name,Complex)]
        V_scenario = point[Variable(volt_scen_name,Complex)]
        Qmin = imag(element.power_min)
        Qmax = imag(element.power_max)
        #println(abs(V_basecase)," ", abs(V_scenario))

        if  -epsilon <= abs(V_scenario) - abs(V_basecase) <= epsilon
            ##Vscen = Vbc => no constraint
            #println("Vscenario = Vbasecase")
            return
        elseif abs(V_scenario) - abs(V_basecase) < - epsilon
            #Vscen < Vbc => Qscen = Qmax
            println(bus, elemid, " Vscenario < Vbasecase")
            slack = Qmax - imag(Sgen_scenario)
            println(Qmax, imag(Sgen_scenario))
            println(slack)
            if abs(slack) > epsilon
                return SortedDict(cc_name => "Vscenario < Vbasecase, reactive power generated should be equal to Qmax")
            end

        elseif abs(V_scenario) - abs(V_basecase) > epsilon
            #Vscn > Vbc => Qscen = Qmin
            println(bus, elemid, " Vscenario > Vbasecase")
            slack = imag(Sgen_scenario) - Qmin
            println(Qmin, imag(Sgen_scenario))
            println(slack)
            if abs(slack) > epsilon
                return SortedDict(cc_name => "Vscenario > Vbasecase, reactive power generated should be equal to Qmin")
            end
        end

    else
        #warn("No reactive power coupling constraints for BaseCase scenario")
        return
    end
end

# scenario = "Scenario1"
# bus = "BUS_1"
# ds_bus = OPFpbs[scenario].ds.bus[bus]
# elemid = "Generator_1"
# elem = ds_bus[elemid]
# elem_formulation = OPFpbs[scenario].mp.node_formulations[bus][elemid]
# point = global_point
# check_feasibility_cc_reactive_power(elem, bus, elemid, elem_formulation, scenario, point,1e-7)


## TODO::
# function check_feasibility_power_balance(bus::String, ds::DataSource, gs::GridStructure, mp::MathematicalProgramming, scenario::String, point::Point, epsilon::Float64)
#
#     bus_elems_formulations = mp.node_formulations[bus]
#     bus_elems_var = mp.node_vars[bus]
#     cstrnames, Snode, lb, ub = get_Snodal(busid, ds.bus[bus], bus_elems_formulations, bus_elems_var)
#     S_balance = Snode + get_Slinks_in(busid, ds, gs, mp) + get_Slinks_out(busid, ds, gs, mp)
#
#
# end
"""
    test_feasibility_GOC(filenames, epsilon::Float64)

Check feasibility at `epsilon` for instance defined by `filenames`
Return a dictonary of violated constraints by scenario : dict( scenario => dict(ctrname => message))

"""

function test_feasibility_GOC(instance_path, epsilon::Float64)
    OPFpbs = load_OPFproblems(GOCInput, instance_path)
    introduce_Sgenvariables!(OPFpbs)
    pb_global = build_globalpb(OPFpbs)
    global_point = solution_point(instance_path)
    dict_scenario_infeasible_ctr = SortedDict{String, SortedDict{String,String}}()
    for (scenario, OPFpb) in OPFpbs
        dict_scenario_infeasible_ctr[scenario] = check_feasibility(scenario, OPFpb, global_point,epsilon)
    end
    return dict_scenario_infeasible_ctr
end

"""
    check_feasibility(scenario::String, OPFpb::Scenario, point::Point, epsilon::Float64)

Checks feasibility at `epsilon` for each constraint generated for in the scenario `scenario` constructed from `OPFpb` evaluated at `point`.
Returns a dictionary of constraints violated SortedDict: ctrname => message of error
"""
function check_feasibility(scenario::String, OPFpb::Scenario, point::Point, epsilon::Float64)
    ds = OPFpb.ds
    mp = OPFpb.mp

    dict_infeasible_ctr = SortedDict{String, String}()
        for (bus, elems) in ds.bus
            bus_elems_formulations = mp.node_formulations[bus]
            for (elemid, elem) in elems
                elem_formulation = bus_elems_formulations[elemid]
                dict2 = check_feasibility(elem, bus, elemid, elem_formulation, scenario, point, epsilon)
                if typeof(dict2) != Void
                    merge!(dict_infeasible_ctr, dict2)
                end
                dict3 = check_feasibility_cc_active_power(elem, bus, elemid, elem_formulation, scenario, point, epsilon)
                if typeof(dict3) != Void
                    merge!(dict_infeasible_ctr, dict3)
                end
                dict4 = check_feasibility_cc_reactive_power(elem, bus, elemid, elem_formulation, scenario, point, epsilon)
                if typeof(dict4) != Void
                    merge!(dict_infeasible_ctr, dict4)
                end
            end
        end

        for (link, elems) in ds.link
            link_elems_formulations = mp.link_formulations[link]
            link_vars = mp.link_vars[link]
            for (elemid, elem) in elems
                elem_formulation = link_elems_formulations[elemid]
                dict2 = check_feasibility(elem, link, elemid, elem_formulation, link_vars, scenario, point, epsilon)
                if typeof(dict2) != Void
                    merge!(dict_infeasible_ctr, dict2)
                end
            end
        end
    return dict_infeasible_ctr
end


function is_basecase_scenario(scenario::String)
    if scenario == basecase_scenario_name()
        return true
    else
        return false
    end
end
