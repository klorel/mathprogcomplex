"""

Construct a csv file containing data for all instances
"""

@everywhere ROOT = pwd()
@everywhere include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))
@everywhere include("para_fct.jl"))


# NOTE: Launch in console using 'julia -p nbparathreads compute_GOC_all_scenario.jl'
# where the nb of cores is likely a good value for nbparathreads

instances_path = joinpath("instances", "GOC")
folders = readdir(instances_path)

##empty dictionaries which will contain data to export in csv
instances_nb_scenario_data = SortedDict{String,Int64}()
instances_nb_contingency_data = SortedDict{String,Int64}()
instances_nb_bus_data = SortedDict{String,Int64}()
instances_nb_links_data = SortedDict{String,Int64}()
instances_nb_violated_constraints = SortedDict{String, Float64}()
instances_nb_cc_Qgen_active = SortedDict{String, Float64}()


##loop on instances to get data
epsilon = 1e-3
for folder in folders
    instances_nb_scenario_data[folder] = get_nb_scenarios(instances_path,folder)
    scenario = readdir(joinpath(instances_path, folder))[1]
    instances_nb_contingency_data[folder] = get_nb_contingencies(folder,scenario)
    instances_nb_bus_data[folder] = get_nb_bus(folder, scenario)
    instances_nb_links_data[folder] = get_nb_links(folder,scenario)
    instances_nb_cc_Qgen_active[folder] = mean(get_nb_cc_reactive_power_active(folder, scenario, epsilon) for scenario in readdir(joinpath(instances_path, folder)) if scenario !="scorepara.csv")
    #instances_nb_violated_constraints[folder] = mean(get_feasibility(folder, scenario, epsilon) for scenario in readdir(joinpath(instances_path, folder)))
end
##create csv file
touch("GOCdata.csv")
f = open("GOCdata.csv", "w")
write(f, "Instance;Nb_scenarios;Nb_contingencies;Nb_bus;Nb_links;Mean nb of violated constraints by GOC solution at $(epsilon);Mean nb of reactive Qgen contingency constraints at $(epsilon)\n")
for (folder, nb_scenarios) in instances_nb_scenario_data
    nb_cont = instances_nb_contingency_data[folder]
    nb_bus = instances_nb_bus_data[folder]
    nb_links = instances_nb_links_data[folder]
    # nb_vc = instances_nb_violated_constraints[folder]
    nb_vc = "NONE"
    nb_cc_act = instances_nb_cc_Qgen_active[folder]
    write(f, "$(folder);$(nb_scenarios);$(nb_cont);$(nb_bus);$(nb_links);$(nb_vc);$(nb_cc_act)\n")
end
close(f)



##functions to get specific data
"""
    get_nb_scenarios(instances_path,folder)

Returns the number of scenarios for an instance `folder` in `instances_path`.
"""
function get_nb_scenarios(instances_path,folder)
    folder_path = joinpath(instances_path, folder)
    scenarios = sort(filter(x->!ismatch(r"\.", x), readdir(folder_path)), by=x->parse(split(x, "_")[2]))
    nb_scenarios = length(scenarios)
    return nb_scenarios
end

"""
    get_nb_bus(folder, scenario)

Returns the number of bus for a `scenario` in instance `folder`.
"""
function get_nb_bus(folder, scenario)
    power_data = getpowerdata([folder, scenario])
    try
        nb_bus =  Int(power_data["totalbus"])
        return nb_bus
    catch
        warn("power_data coming from .mat file does not have key totalbus")
        return 0.0
    end
end

"""
    get_nb_contingencies(folder,scenario)

Returns the number of contingencies for a `scenario` in instance `folder`.
"""
function get_nb_contingencies(folder,scenario)
    contingency_data = getcontingencydata([folder, scenario])
    nb_contingency = size(contingency_data,1)
    return nb_contingency
end

"""
    get_nb_links(folder, scenario)

Returns the number of links for a `scenario` in instance `folder`.
"""
function get_nb_links(folder, scenario)
    power_data = getpowerdata([folder, scenario])
    try
        nb_links =  Int(power_data["totalbranch"])
        return nb_links
    catch
        warn("power_data coming from .mat file does not have key totalbranch")
        return 0.0
    end
end

"""
    get_feasibility(folder, scenario, epsilon)

Returns the number of violated constraints at `epsilon` (except balance constraints) for a `scenario` in instance `folder`.
"""
function get_feasibility(folder, scenario, epsilon)
    try
        violations = test_feasibility_GOC([folder, scenario],epsilon)
        nb_violated_constraints = 0
        for (scenario, dict_infeas) in violations
            nb_violated_constraints += length(dict_infeas)
        end
        return nb_violated_constraints

    catch
        return -1.0
    end
end


"""
    get_nb_cc_reactive_power_active(folder, scenario,epsilon)

Returns the number of active constraints of type contingency Qgen coupling constraints at `epsilon` for a `scenario` in instance `folder`.
"""
function get_nb_cc_reactive_power_active(folder, scenario,epsilon)
    try
        return nb_active_constraints_in_scenario_cc_Qgen([folder, scenario], epsilon)
    catch
        println([folder, scenario])
        return -1.0
    end
end
