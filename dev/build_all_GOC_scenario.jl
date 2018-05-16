include(joinpath(ROOT,"src_PowSysMod", "PowSysMod_body.jl"))

# NOTE: Launch in console using 'julia -p nbparathreads compute_GOC_all_scenario.jl'
# where the nb of cores is a likely good value for nbparathreads

folder = "Phase_0_IEEE14_1Scenario"
work_dir = abspath("..", "data", "data_GOC", folder)
ispath(work_dir) || error("Not a folder $work_dir")

scenarios = SortedSet(setdiff(readdir(work_dir), ["scorepara.csv"]))

for scenario in scenarios
    cur_dir = joinpath(work_dir, scenario)
    info("Working in $cur_dir now.")

    OPFpbs = load_OPFproblems(GOCInput, cur_dir)
    println(keys(OPFpbs))
end