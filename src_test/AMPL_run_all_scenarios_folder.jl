ROOT = pwd()
include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))
using MAT

function solve_GOC_via_AMPL(data_path, folder, scenario)
    folder_path = joinpath(data_path, folder)
    instance_path = joinpath(folder_path, scenario)
    raw = "powersystem.raw"
    gen = "generator.csv"
    con = "contingency.csv"
    rawfile = joinpath(instance_path,raw)
    genfile = joinpath(instance_path, gen)
    contfile = joinpath(instance_path, con)
    OPFpbs = load_OPFproblems(rawfile, genfile, contfile)
    introduce_Sgenvariables!(OPFpbs)
    ## Bulding optimization problem
    pb_global = build_globalpb!(OPFpbs)
    pb_global_real = pb_cplx2real(pb_global)

    init_point = Point()
    init_point_real = cplx2real(init_point)

    ## Exporting real problem
    amplexportpath = joinpath("..","knitro_runs", "$(folder[9:end])_$(scenario)")

    my_timer = @elapsed export_to_dat(pb_global_real, amplexportpath, init_point_real)
    @printf("%-35s%10.6f s\n", "export_to_dat", my_timer)

    _, t_knitro, _ = @timed run_knitro(amplexportpath, joinpath(pwd(),"..","src_ampl"))
    pt_knitro, _ = read_Knitro_output(amplexportpath, pb_global_real)

    solve_result_1, solve_result_2 = read_knitro_info_csvfile(amplexportpath)

    outpath = joinpath(pwd(),"..","solutions", "$folder")
    isdir(outpath) || mkpath(outpath)
    outpath = joinpath(outpath, "$scenario")
    isdir(outpath) || mkpath(outpath)
    write_solutions(OPFpbs, pt_knitro, outpath)

    sol_txt = read_solution_point_GOC(instance_path, outpath)
    pt_txt = cplx2real(sol_txt)
    min_slack = get_minslack(pb_global_real, pt_txt)

    feas,ctr = get_minslack(pb_global_real, pt_knitro)
    obj = get_objective(pb_global_real, pt_knitro)


    return scenario => (solve_result_1, solve_result_2, (feas,ctr), min_slack)
end

function read_args(ARGS)
    if length(ARGS)!=2
        error("First argument must be data_path for example ..\..\data\data_GOC
            Second argument must be instance folder for example Phase_0_IEEE14
        ")
    else
        data_path = ARGS[1]
        folder = ARGS[2]
    end
    return data_path, folder
end

data_path, folder = read_args(ARGS)
folder_path = joinpath(data_path, folder)

scenarios = sort(filter(x->!ismatch(r"\.", x), readdir(folder_path)), by=x->parse(split(x, "_")[2]))
nb_scenarios = length(scenarios)

println("----------> Start para jobs")

# r = build_and_solve_GOC(folder, scenarios[1])
results = pmap(solve_GOC_via_AMPL, [data_path for i=1:length(scenarios)], [folder for i=1:length(scenarios)], scenarios)

println("----------> para jobs done\n")
println("######################################################################################")

filename = joinpath("..","knitro_runs","results_$folder.csv")
touch(filename)

f = open(filename, "w")

write(f, "Scenario;solve_result_1; solve_result_2; min slack from knitro point; min_slack from txt files\n")
nb_scenarios_with_pb = 0
for (scenario, data) in results
    solve_result_1 = data[1]
    solve_result_2 = data[2]
    feas,ctr1 = data[3]
    min_slack,ctr2 = data[4]

    if solve_result_1!=0 || solve_result_2!=0 || feas < -1e-6 || min_slack < -1e-6
        nb_scenarios_with_pb +=1
        println("PB: $scenario not feasible")
    end
    write(f, "$(scenario);$(solve_result_1);$(solve_result_2);$(feas);$(min_slack)\n")
end
if nb_scenarios_with_pb > 0
    println("\nNB OF SCENARIOS NOT FEASIBLE : $nb_scenarios_with_pb/$nb_scenarios. \nSee results_$folder.csv in knitro_runs for more details")
    write(f, "\n ; NB OF SCENARIOS NOT FEASIBLE;$nb_scenarios_with_pb/$nb_scenarios\n")
else
    println("ALL SCENARIOS ($nb_scenarios scenarios) FEASIBLE.\nSee results_$folder.csv in knitro_runs for more details.")
    write(f, "\n ; NB OF SCENARIOS NOT FEASIBLE;0/$nb_scenarios\n")
end
close(f)
