include(joinpath(ROOT,"src_PowSysMod", "PowSysMod_body.jl"))

# NOTE: Launch in console using 'julia -p nbparathreads compute_GOC_all_scenario.jl'
# where the nb of cores is a likely good value for nbparathreads

folder = "Phase_0_IEEE14_1Scenario"
work_dir = joinpath("..", "data", "data_GOC", folder)
isfolder(work_dir) || error("Not a folder $work_dir")



scenarios = sort(filter(x->!ismatch(r"\.", x), readdir(folder_path)), by=x->parse(split(x, "_")[2]))

date = Dates.format(now(), "yy_u_dd_HH_MM_SS")
save_folder = joinpath("knitro_runs", "para_run_$date")
mkdir(save_folder)

println("----------> Start para jobs")

# r = build_and_solve_GOC(folder, scenarios[1])
r = pmap(build_and_solve_GOC, [folder for i=1:length(scenarios)], scenarios, [save_folder for i=1:length(scenarios)])

println("----------> para jobs done")

println("length scenarios : $(length(r))")


# r = SortedDict(r)
#
# tasks_th = SortedDict()
# for (key, val) in r
#     if !haskey(tasks_th, val["thread_id"])
#         tasks_th[val["thread_id"]] = 0
#     end
#     tasks_th[val["thread_id"]] += 1
# end
# println("id count:")
# for (k, v) in tasks_th
#     println("$k → $v")
# end
