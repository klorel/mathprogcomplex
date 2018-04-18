using DataStructures
##read plateforme results for a dataset
submission_number = "297-1523974791"
results_path = joinpath(pwd(),"..", "..", "Plateforme_GOC_results")
logs_path = joinpath(pwd(),"..", "knitro_runs")

lines = readdlm(joinpath(results_path, submission_number, "$(submission_number)_score.csv"), ',')

dataset = lines[3,2]

println(dataset)

scenarios_data = lines[5:end-2,:]

#Case name,Score,Objective Value,Max Violation,Contingency ID of Max Violation,Contingency type of Max Violation,Computation Time (second),
infeas_by_scenario = SortedDict{String,Tuple{Float64, Int64, String}}()

for i in 1:size(scenarios_data,1)
    scenario = scenarios_data[i,1]
    score = scenarios_data[i,2]
    objective = scenarios_data[i,3]
    max_viol = scenarios_data[i,4]
    cont_id_max_viol = scenarios_data[i,5]
    ctr_type_max_viol = scenarios_data[i,6]
    computation_time = scenarios_data[i,7]

    infeas_by_scenario[scenario] = (max_viol, cont_id_max_viol, ctr_type_max_viol)
    if max_viol >= 1e-6
        println("\n$scenario :  $max_viol, $cont_id_max_viol, $ctr_type_max_viol")
        log_folder = joinpath(logs_path, "$(dataset[9:end])_$scenario")
        log_file = filter(x->ismatch(r".log", x), readdir(log_folder))[1]
        f = open(joinpath(log_folder, log_file),"r")
        loglines = readlines(f)
        close(f)
        printline = false
        for line in loglines
            if ismatch(r"EXIT:", line)
                printline = true
            end
            if printline
                println(line)
                if line[end]== '.'
                    printline = false
                end
            end
        end

    end
end
