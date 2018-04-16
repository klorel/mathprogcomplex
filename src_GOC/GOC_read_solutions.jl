function solution_point(instance_path::String)
    solution1 = readdlm(joinpath(instance_path,"solution1.txt"), ',')
    solution2 = readdlm(joinpath(instance_path,"solution2.txt"), ',')
    basecase_generator_solution = read_solution1(solution1)
    contingency_generator_solution, bus_solution, delta_solution, line_flow_solution = read_solution2(solution2)
    power_data = getpowerdata(instance_path)
    index = get_bus_index(power_data)
    participation_factors = create_participation_factors_dict(instance_path, index)
    global_point = create_global_point(solution1, solution2, participation_factors, index)
    # Add delta solution
    delta_point = create_delta_solution(solution2)
    binary_point = compute_binary_values(basecase_generator_solution, bus_solution, power_data)
    global_point = merge(global_point, delta_point,binary_point)
    return global_point
end

function read_solution_point_GOC(instance_path::String, solution_path::String)
    solution1 = readdlm(joinpath(solution_path,"solution1.txt"), ',')
    solution2 = readdlm(joinpath(solution_path,"solution2.txt"), ',')
    basecase_generator_solution = read_solution1(solution1)
    contingency_generator_solution, bus_solution, delta_solution, line_flow_solution = read_solution2(solution2)
    power_data = getpowerdata(instance_path)
    index = get_bus_index(power_data)
    participation_factors = create_participation_factors_dict(instance_path, index)
    global_point = create_global_point(solution1, solution2, participation_factors, index)
    # Add delta solution
    delta_point = create_delta_solution(solution2)
    binary_point = compute_binary_values(basecase_generator_solution, bus_solution, power_data)
    println(binary_point)
    global_point = merge(global_point, delta_point,binary_point)
    return global_point
end


function create_voltagesolutionpoint(bus_solution, index)
    point = Point()
    for i in 1:size(bus_solution,1)
        num_bus = index[bus_solution[i,2]]
        module_volt = bus_solution[i,3]
        value = module_volt*exp(im*bus_solution[i,4]*pi/180)
        if bus_solution[i,1] == 0
            add_coord!(point, Variable(variable_name("VOLT", bus_name(num_bus), "",basecase_scenario_name()), Complex), value)
        else
            id_contingency = bus_solution[i,1]
            scenario_name = scenarioname(id_contingency)
            add_coord!(point, Variable(variable_name("VOLT", bus_name(num_bus), "", scenario_name), Complex), value)
        end
    end
    return point
end


function gen_solution_basecase(basecase_generator_solution,index)
    point = Point()
    for i in 1:size(basecase_generator_solution,1)
        num_bus = index[basecase_generator_solution[i,1]]
        gen_id = basecase_generator_solution[i,2]
        if typeof(gen_id)!=Int64
            gen_id = matchall(r"\d+" ,gen_id)[1]
        else
            gen_id = Int(gen_id)
        end
        Pg = basecase_generator_solution[i,3]
        Qg = basecase_generator_solution[i,4]
        add_coord!(point, Variable(variable_name("Sgen", bus_name(num_bus), generator_name(gen_id), basecase_scenario_name()), Complex), Pg + im*Qg)
    end
    return point
end

function gen_solution_contingencies(basecase_generator_solution,contingency_generator_solution,delta_solution, participation_factors, index)
    point = Point()

    ##reactive power
    for i in 1:size(contingency_generator_solution,1)
        id_contingency = contingency_generator_solution[i,1]
        scenario_name = scenarioname(id_contingency)
        gen_id = contingency_generator_solution[i,4]
        num_bus = index[Int(contingency_generator_solution[i,3])]
        if typeof(gen_id)!=Int64
            gen_id = matchall(r"\d+" ,gen_id)[1]
        else
            gen_id = Int(gen_id)
        end
        Qg = contingency_generator_solution[i,5]
        add_coord!(point, Variable(variable_name("Sgen", bus_name(num_bus), generator_name(gen_id), scenario_name), Complex), im*Qg)
    end


    #active power
    for l in 1:size(delta_solution,1)
        id_contingency = delta_solution[l,1]
        scenario_name = scenarioname(id_contingency)
        for i in 1:size(basecase_generator_solution,1)
            num_bus = index[Int(basecase_generator_solution[i,1])]
            gen_id = basecase_generator_solution[i,2]
            if typeof(gen_id)!=Int64
                gen_id = matchall(r"\d+" ,gen_id)[1]
            else
                gen_id = Int(gen_id)
            end
            Pg = basecase_generator_solution[i,3]

            var = Variable(variable_name("Sgen",bus_name(num_bus), generator_name(gen_id), scenario_name), Complex)

            add_coord!(point, var, Pg + delta_solution[l,2] * participation_factors[bus_name(num_bus)][generator_name(gen_id)])
        end
    end
    return point
end



function create_participation_factors_dict(instance_path, index)
    path_folder = instance_path
    if !("generator.csv" ∈ readdir(path_folder))
        error("missing file : generator.csv must be in $path_folder")
    end
    generator_data = readdlm(joinpath(path_folder,"generator.csv"),',',skipstart=1)

    participation_factors = Dict{String,Dict{String,Float64}}()
    for line in 1:size(generator_data,1)
        term_id = generator_data[line,3]
        if term_id == 9
            num_bus = index[Int(generator_data[line,1])]
            gen_id = generator_data[line,2]
            if typeof(gen_id)!=Float64
                gen_id = matchall(r"\d+" ,gen_id)[1]
            else
                gen_id = Int(gen_id)
            end
            value =  generator_data[line,4]
            if !haskey(participation_factors,bus_name(num_bus))
                participation_factors[bus_name(num_bus)] = Dict{String,Float64}()
            end
            participation_factors[bus_name(num_bus)][generator_name(gen_id)] = value
        end
    end
    return participation_factors
end

function create_delta_solution(solution2)
    point = Point()
    contingency_generator_solution, bus_solution, delta_solution, line_flow_solution = read_solution2(solution2)
    for i in 1:size(delta_solution,1)
        scenario = scenarioname(delta_solution[i,1])
        add_coord!(point, Variable(get_delta_varname(scenario),Real), delta_solution[i,2])
    end
    return point
end


function create_global_point(solution1,solution2, participation_factors,index)
    basecase_generator_solution = read_solution1(solution1)
    contingency_generator_solution, bus_solution, delta_solution, line_flow_solution = read_solution2(solution2)
    voltage_point = create_voltagesolutionpoint(bus_solution,index)
    genpoint =  gen_solution_basecase(basecase_generator_solution,index)

    genpoint_contingencies = gen_solution_contingencies(basecase_generator_solution,contingency_generator_solution,delta_solution, participation_factors,index)

    return merge(voltage_point, genpoint, genpoint_contingencies)
end


function compute_binary_values(basecase_generator_solution, bus_solution, power_data)
    index = get_bus_index(power_data)
    generators = collect(basecase_generator_solution[:,1])
    generators = [ index[gen] for gen in generators]
    module_v_basecase = Dict{Int64, Float64}()
    module_v_scenarios = Dict{String,Dict{Int64, Float64}}()
    for i in 1:size(bus_solution,1)
        num_bus = index[bus_solution[i,2]]
            if bus_solution[i,1] == 0
                module_v_basecase[num_bus] = bus_solution[i,3]
            else
                id_contingency = bus_solution[i,1]
                scenario_name = scenarioname(id_contingency)
                if !haskey(module_v_scenarios, scenario_name)
                    module_v_scenarios[scenario_name] = Dict{Int64, Float64}()
                end
                module_v_scenarios[scenario_name][num_bus] = bus_solution[i,3]
            end
    end
    point = Point()
    for (scenario, dict_modules) in module_v_scenarios
        for num_bus in generators
            if abs(dict_modules[num_bus]^2 - module_v_basecase[num_bus]^2) < get_GOC_Volt_ϵ()
                add_coord!(point, Variable(get_binEq_varname(scenario, basecase_scenario_name(), bus_name(num_bus)),Bool), 1.0)
                add_coord!(point, Variable(get_binInf_varname(basecase_scenario_name(),scenario, bus_name(num_bus)),Bool), 0.0)
                add_coord!(point, Variable(get_binInf_varname(scenario, basecase_scenario_name(),bus_name(num_bus)),Bool), 0.0)
            elseif dict_modules[num_bus]^2 - module_v_basecase[num_bus]^2 > get_GOC_Volt_ϵ()
                add_coord!(point, Variable(get_binInf_varname(basecase_scenario_name(),scenario, bus_name(num_bus)),Bool), 1.0)
                add_coord!(point, Variable(get_binInf_varname(scenario, basecase_scenario_name(),bus_name(num_bus)),Bool), 0.0)
                add_coord!(point, Variable(get_binEq_varname(scenario, basecase_scenario_name(), bus_name(num_bus)),Bool), 0.0)
            else
                add_coord!(point, Variable(get_binInf_varname(scenario, basecase_scenario_name(),bus_name(num_bus)),Bool), 1.0)
                add_coord!(point, Variable(get_binInf_varname(basecase_scenario_name(),scenario, bus_name(num_bus)),Bool), 0.0)
                add_coord!(point, Variable(get_binEq_varname(scenario, basecase_scenario_name(), bus_name(num_bus)),Bool), 0.0)
            end
        end
    end
    # println(point)
    return point
end
# filenames = ["Phase_0_IEEE14_1Scenario", "scenario_1"]
# filenames = ["Phase_0_RTS96", "scenario_10"]
# nb_active_constraints_in_scenario_cc_Qgen(filenames, 1e-4)
# folder = "Phase_0_IEEE14"
# scenario = "scenario_10"
# epsilon = 1e-4
# get_nb_cc_reactive_power_active(folder, scenario, epsilon)

##########################################
#READ
#########################################

function read_solution1(solution1::Array{Any,2})
    nb_lines = size(solution1,1)
    #description_line = Array{Any}()
    line_start = 0
    line_end = 0
    for l in 1:nb_lines
        if solution1[l,1] == "--generation dispatch"
            #description_line = solution1[l+1,1:4]
            line_start = l + 2
            break
        end
    end
    i = line_start
    while solution1[i,1] != "--end of generation dispatch"
        if i==nb_lines
            error("no --end of generation dispatch in input solution1")
        else
        i +=1
        end
    end
    line_end = i-1
    basecase_generator_solution = solution1[line_start:line_end,1:4]
    return basecase_generator_solution
end

#
# basecase_generator_solution_description = ["bus id", "unit id", "pg(MW)", "qg(MVar)"]
# contingency_generator_solution_description = ["conID", "genID", "busID", "unitID", "q(MW)"]
# contingency_bus_solution_description = ["contingency id", "bus id", "v(pu)", "theta(deg)"]
# contingency_delta_solution_description = ["contingency id", "Delta(MW)"]
# contingency_line_flow_solution_description = ["contingency id", "line id", "origin bus id", "destination bus id", "circuit id", "p_origin(MW)", "q_origin(MVar)", "p_destination(MW)", "q_destination(MVar)"]
#




function read_solution2(solution2::Array{Any,2})
    nblines = size(solution2,1)

    i = 1
    while solution2[i,1] != "--contingency generator"
        if i==nblines
            error("no --contingency generator in input solution2")
        else
        i +=1
        end
    end
    #println(solution2[i+1,1:5])
    line_start_ctg, line_end_ctg = findnumarray(i+2, solution2)

    i = line_end_ctg + 1
    while solution2[i,1] != "--bus"
        if i==nblines
            error("no --bus in input solution2")
        else
        i +=1
        end
    end
    #println(solution2[i+1,1:4])
    line_start_bus, line_end_bus = findnumarray(i+2, solution2)

    i = line_end_bus + 1
    while solution2[i,1] != "--Delta"
        if i==nblines
            error("no --Delta in input solution2")
        else
        i +=1
        end
    end
    #println(solution2[i+1,1:2])
    line_start_delta, line_end_delta = findnumarray(i+2, solution2)

    i = line_end_delta + 1
    while solution2[i,1] != "--line flow"
        if i==nblines
            error("no --line flow in input solution2")
        else
        i +=1
        end
    end
    #println(solution2[i+1,1:9])
    line_start_line, line_end_line = findnumarray(i+2, solution2)

    if line_end_line + 1 == nblines
        contingency_generator_solution = solution2[line_start_ctg:line_end_ctg,1:5]
        bus_solution = solution2[line_start_bus:line_end_bus,1:4]
        delta_solution = solution2[line_start_delta:line_end_delta,1:2]
        line_flow_solution = solution2[line_start_line:line_end_line,1:9]

        return contingency_generator_solution, bus_solution, delta_solution, line_flow_solution
    else
        error()
    end
end

function findnumarray(line_start, data)
    i = line_start
    name = data[i-2,1][3:end]
    nb_lines = size(data,1)
    while data[i,1]!="--end of "*name
        if i==nb_lines
            error("no --end of $name in input solution2")
        else
        i +=1
        end
    end
    line_end = i-1
    line_start, line_end
end
