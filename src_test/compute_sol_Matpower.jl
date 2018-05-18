ROOT = pwd()
include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))

function build_and_solve_matpower_instance(instance_path)
    # typeofinput = MatpowerSimpleInput
    typeofinput = MatpowerInput
    OPFpbs = load_OPFproblems(typeofinput, instance_path)
    ## Bulding optimization problem
    pb_global = build_globalpb!(OPFpbs)
    pb_global_real = pb_cplx2real(pb_global)
    ## Reading GOC initial point
    init_point = Point()
    init_point_real = cplx2real(init_point)
    ## Exporting real problem
    instance_name = splitdir(instance_path)[end]
    amplexportpath = joinpath("..","knitro_runs", "$(instance_name[1:end-2])")
    my_timer = @elapsed export_to_dat(pb_global_real, amplexportpath, init_point_real)
    @printf("%-35s%10.6f s\n", "export_to_dat", my_timer)
    _, t_knitro, _ = @timed run_knitro(amplexportpath, joinpath(pwd(),"..","src_ampl"))
    sol = open(joinpath("..", "..", "data", "data_Matpower","solutions","sol.txt"), "r")
    f = open(joinpath("..", "..", "data", "data_Matpower","solutions","$(instance_name[1:end-2])_sol.txt"), "w")
    for l in readlines(sol)
        write(f, l)
        write(f, "\n")
    end
    close(sol)
    close(f)

end


data_path = joinpath("..","..","data","data_Matpower","matpower")
instances = sort(readdir(data_path))
# instances = ["LMBM3.m", "WB2.m", "WB3.m", "WB5.m", "case118.m", "case1354pegase.m", "case13659pegase.m", "case14.m", "case145.m", "case1888rte.m", "case1951rte.m", "case2383wp.m", "case24_ieee_rts.m", "case2736sp.m", "case2737sop.m", "case2746wop.m", "case2746wp.m", "case2848rte.m", "case2868rte.m", "case2869pegase.m", "case30.m", "case300.m", "case300mod.m", "case3012wp.m", "case3120sp.m", "case3375wp.m", "case33bw.m", "case39.m", "case57.m", "case6468rte.m", "case6470rte.m","case6495rte.m","case6515rte.m", "case6ww.m", "case89pegase.m", "case9.m", "case9241pegase.m", "case9target.m", "case_ieee30.m", "case_illinois200.m"]
instances = ["LMBM3.m", "WB2.m", "WB3.m", "WB5.m", "case118.m", "case1354pegase.m", "case14.m", "case145.m", "case30.m", "case300.m", "case300mod.m",  "case39.m", "case57.m", "case6ww.m", "case89pegase.m", "case9.m", "case9target.m", "case_ieee30.m", "case_illinois200.m"]
# instances = ["case9.m"#=,"case14.m", "case30.m"=#]
instances = ["case9.m","case1888rte.m"]
println(instances)
nb_instances = length(instances)


println("----------> Start para jobs")

results = pmap(build_and_solve_matpower_instance, [joinpath(data_path, instance) for instance in instances])
