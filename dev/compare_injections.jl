

#
function get_Sgen_GOC(instance_path)
    typeofinput = GOCInput
    OPFpbs = load_OPFproblems(typeofinput, instance_path)

    ## Introducing coupling constraints on generator output
    (typeofinput != GOCInput) || introduce_Sgenvariables!(OPFpbs)

    ## Bulding optimization problem
    pb_global = build_globalpb!(OPFpbs)
    pb_global_real = pb_cplx2real(pb_global)

    ## Reading GOC initial point
    init_point = Point()
    (typeofinput != GOCInput) || (init_point = solution_point(instance_path))
    init_point_real = cplx2real(init_point)

    ## Exporting real problem
    folder,scenario = split(instance_path, '\\')[end-1:end]
    folder = String(folder)
    scenario = String(scenario)
    amplexportpath = joinpath("..","knitro_runs", "$(folder[9:end])_$(scenario)")
    pt_knitro, pt_GOC = read_Knitro_output(amplexportpath, pb_global_real)

    obj_GOC = get_objective(pb_global_real, pt_GOC)
    obj_knitro = get_objective(pb_global_real, pt_knitro)
    diff_rel_obj = (obj_GOC-obj_knitro)/obj_GOC
    Sgens = Dict()
    scenario = basecase_scenario_name()
    scenario = "Scen1"
    _, _, bc_prod_vars_knitro = get_splitted_Cpt(pt_knitro, scenario)
    _, _, bc_prod_vars_GOC = get_splitted_Cpt(pt_GOC, scenario)

    return bc_prod_vars_knitro - bc_prod_vars_GOC, diff_rel_obj
end

ROOT = "D:\\repo\\complex-modeler\\dev"
include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))
instances_folder = "data_GOC"
folder = "Phase_0_IEEE14"
scenario = "scenario_2"
instance_path = joinpath(ROOT, "..", instances_folder,folder, scenario)

get_Sgen_GOC(instance_path)

function treat_sol(csv_file, output)
    tab = readdlm(csv_file, ';')
    nb_lines = size(tab,1)
    vars = Dict{String, Tuple{Float64,Float64}}()
    for i in 1:nb_lines
        varname = tab[i,1]
        if ismatch(r"Sgen", varname) || ismatch(r"Delta", varname)
            vars[varname] = (tab[i,2],tab[i,3])
        end
    end
    return vars
end



function gendiff_output(instance,vars,output)
    dict_scenario_var_real = Dict{String,Set{String}}()
    dict_scenario_var_im = Dict{String,Set{String}}()
    dict_scenario_var_delta = Dict{String,String}()

    for (var,values) in vars
        if ismatch(r"Re",var)
            scenario = String(split(var,"_")[1])
            if !haskey(dict_scenario_var_real,scenario)
                dict_scenario_var_real[scenario] = Set{String}()
            end
            push!(dict_scenario_var_real[scenario],var)
        elseif ismatch(r"Im",var)
            scenario = String(split(var,"_")[1])
            if !haskey(dict_scenario_var_im,scenario)
                dict_scenario_var_im[scenario] = Set{String}()
            end
            push!(dict_scenario_var_im[scenario],var)
        else
            scenario = String(split(var,"_")[1])
            dict_scenario_var_delta[scenario] = var
        end
    end

    scenario_diffs_active, ~ = create_dict_diffs(dict_scenario_var_real, vars)
    scenario_diffs_reactive, nb_generator_per_scenario = create_dict_diffs(dict_scenario_var_im, vars)

    scenario_diffs_delta = Dict{String,Float64}()
    for (scenario,varname) in dict_scenario_var_delta
        knitro_val, GOC_val = vars[varname]
        diff_rel = abs((knitro_val-GOC_val)/GOC_val)
        scenario_diffs_delta[scenario] = diff_rel
    end

    return Dict(scenario => (nb_generator_per_scenario[scenario],prop,scenario_diffs_reactive[scenario]) for (scenario,prop) in scenario_diffs_active)
end



function create_dict_diffs(dict_scenario_var, vars)
    scenario_diffs = Dict{String,Float64}()
    nb_generator_per_scenario = Dict{String,Int64}()
    for (scenario, varnames) in dict_scenario_var
        nb_generator = length(varnames)
        nb_generator_per_scenario[scenario] = nb_generator
        nb_diff = 0
        for varname in varnames
            knitro_val, GOC_val = vars[varname]
            diff_rel = abs((knitro_val-GOC_val)/GOC_val)
            if diff_rel >= 0.01
                nb_diff +=1
            end
        end
        scenario_diffs[scenario] = nb_diff/nb_generator
    end
    return scenario_diffs, nb_generator_per_scenario
end


ROOT = pwd()
include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))
instances_folder = "data_GOC"
folder = "Phase_0_OriginalDataset_IEEE14"
scenario = "scenario_1"


results = Dict{String,Dict{String, Tuple{Int64, Float64, Float64}}}()
for i in 1:75
    scenario = "scenario_$i"
    instance_path = joinpath(ROOT, "..", instances_folder,folder, scenario)
    instance = "$(folder[9:end])_$(scenario)"
    csv_file = joinpath(ROOT, "..","knitro_runs","$(folder[9:end])_$(scenario)", "solution_GOC_point.csv")
    output = joinpath(ROOT, "..","knitro_runs","$(folder[9:end])_$(scenario)", "compare_injections.csv")
    vars = treat_sol(csv_file, output)
    output2 = joinpath(ROOT, "$(folder)_compare_injections_by_scenario.csv")
    dict = gendiff_output(instance,vars,output2)
    results[scenario] = dict
end

output = "$(folder)_results_diff_injections.csv"
touch(output)
f = open(output, "w")

write(f, "#Instance;Scenario;Contingency;Nb generators;Percentage of generators with a significant difference in active power;Percentage of generators with a significant difference in reactive power\n")

for (scenario, dict) in results
    for (contingency, values) in dict
        write(f,"$folder;$scenario;$contingency;$(values[1]);$(values[2]);$(values[3])\n")
    end
end
close(f)
