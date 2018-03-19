include(joinpath("..","src_PowSysMod", "PowSysMod_body.jl"))
include("global_test.jl")

typeofinput = GOCInput
folders = ["Phase_0_IEEE14", "Phase_0_RTS96", "Phase_0_Modified_RTS96", "Phase_0_Modified_IEEE14"]

folder = "Phase_0_Modified_RTS96"
path_folder = joinpath(pwd(),"..","data_GOC", folder)
nb_scenarios = length(readdir(path_folder))-1
# scenarios = Set(readdir(joinpath("instances", "GOC", folder)))
# for scenario in scenarios
#     if ismatch(r".csv", scenario)
#         pop!(scenarios, scenario)
#     end
# end
results = pmap(build_and_solve_instance, [typeofinput for i in 1:nb_scenarios], [joinpath(path_folder, scenario) for scenario in readdir(path_folder) if scenario!="scorepara.csv"])

filename = joinpath("resultsROADEFv2_$folder.csv")
touch(filename)

f = open(filename, "w")

write(f, "Folder;Scenario;Nb_variables;Nb_constraints;Objective;Feasibility;Time build export;Time Knitro resolution;Difference with GOC point\n")
for (scenario, data) in results
    nb_variables = data[1]
    nb_constraints = data[2]
    obj = data[3]
    feas = data[4]
    t_buildexport = data[5]
    t_knitro = data[6]
    diff_GOC = data[7]

    write(f, "$(folder);$(scenario);$(nb_variables);$(nb_constraints);$(obj);$(feas);$(t_buildexport);$(t_knitro);$(diff_GOC)\n")
end
close(f)
