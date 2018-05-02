ROOT = pwd()
include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))

folder = "Phase_0_RTS96"
scenarios = ["scenario_$i" for i in 1:100]

results_path = joinpath("..", "knitro_runs")

solve_results = Dict{String, Tuple{Any,Any,Any,Any}}()

nb_scenarios_pb = 0


for scenario in scenarios
    log_path = joinpath(results_path, "$(folder[9:end])_$scenario")
    solve_result_1, solve_result_2, solve_result_3, scaling = read_knitro_info_csvfile(log_path)
    solve_results[scenario] = (solve_result_1, solve_result_2, solve_result_3, scaling)
    if solve_result_1 != 0 || solve_result_2 != 0 || solve_result_3 != 0
        println("scenario $scenario not solved: step 1: $solve_result_1, step 2: $solve_result_2, step 3: $solve_result_3\n")
        nb_scenarios_pb +=1
    end
end

println("nb of scenarios not feasible: $nb_scenarios_pb")
