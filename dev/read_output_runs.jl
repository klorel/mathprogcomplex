ROOT = pwd()
include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))
include("parse_knitro_logs.jl")

folder = "Phase_0_RTS96"
scenarios = ["scenario_$i" for i in 51:100]

results_path = joinpath("..", "knitro_runs")

# solve_results = Dict{String, Tuple{Any,Any,Any,Any}}()
#
# nb_scenarios_pb = 0
#
#
# for scenario in scenarios
#     log_path = joinpath(results_path, "$(folder[9:end])_$scenario")
#     solve_result_1, solve_result_2, solve_result_3, scaling = read_knitro_info_csvfile(log_path)
#     solve_results[scenario] = (solve_result_1, solve_result_2, solve_result_3, scaling)
#     if solve_result_1 != 0 || solve_result_2 != 0 || solve_result_3 != 0
#         println("scenario $scenario not solved: step 1: $solve_result_1, step 2: $solve_result_2, step 3: $solve_result_3\n")
#         nb_scenarios_pb +=1
#     end
# end
#
# println("nb of scenarios not feasible: $nb_scenarios_pb")

f = open(joinpath(results_path, "results", "RTS96phase1_comparison_scaling.csv"), "w")
write(f, "scenario ; nb_iter no scaling; nb_iter scaling 1; CPU time no scaling; CPU time scaling 1\n")

for scenario in scenarios
    root = pwd()
    repo_path = joinpath(results_path, "$(folder[9:end])_$scenario")
    outlog = "phase1_noscaling.log"
    cd(repo_path)
    run(`cmd /c ampl real_minlp.run '>' $(outlog)`)
    nb_it, feaserror, opterror, cpu = parse_log(outlog)
    log_scaling1 = ""
    for file in readdir(pwd())
        if ismatch(r"Knitro_18_May_0", file)
            log_scaling1 = file
        end
    end
    nb_it1, feaserror1, opterror1, cpu1 = parse_log(log_scaling1)
    println("$scenario : CPU no scaling : $(cpu[1]) ; CPU scaling 1 : $(cpu1[1])\n")
    write(f, "$scenario;$(nb_it[1]);$(nb_it1[1]);$(cpu[1]);$(cpu1[1])\n")
    cd(root)
end

close(f)
