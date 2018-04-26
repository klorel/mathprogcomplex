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

    solve_result_1, opterror1, solve_result_2, opterror2, solve_result_3, opterror3 = read_knitro_info_csvfile(amplexportpath)

    outpath = joinpath(pwd(),"..","solutions", "$folder")
    isdir(outpath) || mkpath(outpath)
    outpath = joinpath(outpath, "$scenario")
    isdir(outpath) || mkpath(outpath)
    write_solutions(OPFpbs, pt_knitro, outpath)

    sol_txt = read_solution_point_GOC(instance_path, outpath)
    pt_txt = cplx2real(sol_txt)
    # min_slack = get_minslack(pb_global_real, pt_txt)
    # rel_slacks_txt = get_relative_slacks(pb_global_real, pt_txt)
    min_slack = get_relativemaxslack(pb_global_real, pt_txt)
    infeas_ctr_txt = get_nb_infeasible_ctr_by_ctrtype(pb_global_real, pt_txt, 1e-6)

    # feas,ctr = get_minslack(pb_global_real, pt_knitro)
    feas,ctr = get_relativemaxslack(pb_global_real, pt_knitro)
    # rel_slacks_knitro = get_relative_slacks(pb_global_real, pt_knitro)
    infeas_ctr_knitro = get_nb_infeasible_ctr_by_ctrtype(pb_global_real, pt_knitro, 1e-6)
    # println(infeas_ctr_knitro)
    obj = get_objective(pb_global_real, pt_knitro)


    return scenario => (solve_result_1, solve_result_2, solve_result_3, (feas,ctr), min_slack, opterror1, opterror2,infeas_ctr_knitro)
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

solve_result_num_msg = Dict{Int64, String}()
solve_result_num_msg[0] = "Locally optimal or satisfactory solution."
solve_result_num_msg[200] = "Convergence to an infeasible point. Problem may be locally infeasible."
solve_result_num_msg[201] = "Relative change in infeasible solution estimate < xtol."
solve_result_num_msg[202] = "Current infeasible solution estimate cannot be improved."
solve_result_num_msg[410] = "Iteration limit reached. Current point is infeasible."

nb_scenarios_per_num = Dict( i => (0,0) for i in keys(solve_result_num_msg))

date = Dates.format(now(), "yy_u_dd_HH_MM_SS")
filename = joinpath("..","knitro_runs","results_$(folder)_$(date).csv")
touch(filename)

f = open(filename, "w")

write(f, "Scenario;solve_result_1; solve_result_2; solve_result_3;max relative slack from knitro point; ctr associated; max relative slack from txt files; ctr associated; opterror1 ; opterror2\n")
filename = joinpath("..","knitro_runs","infeas_by_ctr_$(folder)_$(date).csv")
f2 = open(filename, "w")

ctrtypes = [("BALANCE", "Re"),
            ("BALANCE", "Im"),
          (get_VoltM_cstrname(),"Re"),
          (get_GenBounds_cstrname(),"Re"),
          (get_GenBounds_cstrname(),"Im"),
          (get_NullImpVolt_cstrname(),"Re"),
          (get_VoltBinDef_upper(),"Re"),
          (get_VoltBinDef_lower(),"Re"),
          (get_VoltBinDef_complement(),"Re"),
          (get_CC_active_cstrname(),"Re"),
          (get_CC_reactiveupper_cstrname(),"Re"),
          (get_CC_reactivelower_cstrname(),"Re"),
          (get_Smax_orig_cstrname(),"Re"),
          (get_Smax_dest_cstrname(), "Re")]


write(f2, "Scenario;BALANCE,Re;BALANCE,Im;$(get_VoltM_cstrname()),Re;$(get_GenBounds_cstrname()),Re;$(get_GenBounds_cstrname()),Im;$(get_NullImpVolt_cstrname()),Re;$(get_VoltBinDef_upper()),Re;$(get_VoltBinDef_lower()),Re;$(get_VoltBinDef_complement()),Re;$(get_CC_active_cstrname()),Re;$(get_CC_reactiveupper_cstrname()),Re;$(get_CC_reactivelower_cstrname()),Re;$(get_Smax_orig_cstrname()),Re;$(get_Smax_dest_cstrname()), Re\n")


nb_scenarios_with_pb = 0
for (scenario, data) in results
    solve_result_1 = data[1]
    nb1, nb2 = nb_scenarios_per_num[solve_result_1]
    nb_scenarios_per_num[solve_result_1] = (nb1 + 1, nb2)
    solve_result_2 = data[2]
    nb1, nb2 = nb_scenarios_per_num[solve_result_2]
    nb_scenarios_per_num[solve_result_2] = (nb1, nb2+1)
    solve_result_3 = data[3]
    feas,ctr1 = data[4]
    min_slack,ctr2 = data[5]
    opterror1 = data[6]
    opterror2 = data[7]
    infeas_ctr_knitro = data[8]

    if solve_result_1!=0 || solve_result_2!=0 || solve_result_3!=0 || feas > 1e-6 || min_slack > 1e-6
        nb_scenarios_with_pb +=1
        println("PB: $scenario not feasible")
    end
    write(f, "$(scenario);$(solve_result_1);$(solve_result_2);$(solve_result_3);$(feas);$(ctr1);$(min_slack);$(ctr2);$(opterror1);$(opterror2)\n")
    write(f2, "$(scenario)")

    for ctr in ctrtypes
        if !haskey(infeas_ctr_knitro, ctr)
            value = 0
        else
            value = infeas_ctr_knitro[ctr]
        end
        write(f2, ";$value")
    end
    write(f2, "\n")
end
close(f2)
if nb_scenarios_with_pb > 0
    println("\nNB OF SCENARIOS NOT FEASIBLE : $nb_scenarios_with_pb/$nb_scenarios. \nSee results_$folder.csv in knitro_runs for more details")
    write(f, "\n ; NB OF SCENARIOS NOT FEASIBLE;$nb_scenarios_with_pb/$nb_scenarios\n")
    write(f, "\nNB OF SCENARIOS WITH CODE:; Code ;Phase 1 ;Phase2; Message\n ")
    for (solve_result, (nb1,nb2)) in nb_scenarios_per_num
        write(f, " ;$solve_result;$nb1;$nb2;$(solve_result_num_msg[solve_result])\n ")
    end
else
    println("ALL SCENARIOS ($nb_scenarios scenarios) FEASIBLE.\nSee results_$folder.csv in knitro_runs for more details.")
    write(f, "\n ; NB OF SCENARIOS NOT FEASIBLE;0\n")
end
close(f)
