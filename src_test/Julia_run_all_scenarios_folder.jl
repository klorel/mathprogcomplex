ROOT = pwd()
include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))
using MAT

function solve_GOC_via_Julia(data_path, folder, scenario)
    originalSTDOUT = STDOUT
    date = Dates.format(now(), "yy_u_dd_HH_MM_SS")
    outpath = joinpath("..","JuMP_runs","$(folder)_$(scenario)")
    isdir(outpath) || mkpath(outpath)
    outlog = open(joinpath(outpath, "MyJulia1_$(date).log"), "w")
    redirect_stdout(outlog)

    println("Julia/JuMP test")
    folder_path = joinpath(data_path, folder)
    instance_path = joinpath(folder_path, scenario)
    outpath = joinpath("..","JuMP_runs","$(folder)_$(scenario)")
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
    ##convert to JuMP model
    ##phase 1 : continuous relaxation
    println("
  ##############################################\n
  Phase 1 : resolution  of continuous relaxation\n
  ##############################################")

    mysolver = KnitroSolver(KTR_PARAM_OUTLEV=3,
                            KTR_PARAM_MAXIT=600,
                            KTR_PARAM_SCALE=0,
                            KTR_PARAM_FEASTOL=1e-6,
                            KTR_PARAM_OPTTOL=1e-3,
                            KTR_PARAM_FEASTOLABS=1.0,
                            KTR_PARAM_OPTTOLABS=1.0,
                            KTR_PARAM_BAR_INITPT=2,
                            KTR_PARAM_PRESOLVE=0,
                            KTR_PARAM_HONORBNDS=0,
                            KTR_PARAM_MIP_INTVAR_STRATEGY=1)
    tic()
    my_timer = @elapsed m, variables_jump, ctr_jump, ctr_exp = get_JuMP_cartesian_model(pb_global_real, mysolver)
    @printf("%-35s%10.6f s\n", "get_JuMP_cartesian_model", my_timer)
    toc()
    #resolution
    status1 = solve(m)
    #
    # if status1 == :Optimal
    #   solve_result_1 = 0
    # else
    #   solve_result_1 = 1000
    # end



    ##phase 2 : resolution with complementary constraints + initial point = solution continuous relaxation
    println("
  ##############################################################################################\n
  Phase 2 : resolution with complementary constraints from the solution of continuous relaxation\n
  ##############################################################################################")

    mysolver2 = KnitroSolver(KTR_PARAM_OUTLEV=3,
                            KTR_PARAM_MAXIT=600,
                            KTR_PARAM_SCALE=0,
                            KTR_PARAM_FEASTOL=1e-6,
                            KTR_PARAM_OPTTOL=1e-3,
                            KTR_PARAM_FEASTOLABS=1.0,
                            KTR_PARAM_OPTTOLABS=1.0,
                            KTR_PARAM_BAR_INITPT=2,
                            KTR_PARAM_PRESOLVE=0,
                            KTR_PARAM_HONORBNDS=0,
                            KTR_PARAM_MIP_INTVAR_STRATEGY=2)

    for (varname, varjump) in variables_jump
      setvalue(varjump, getvalue(varjump))
    end
    setsolver(m, mysolver2)
    status2 = solve(m)

    # if status2 == :Optimal
    #   solve_result_2 = 0
    # else
    #   solve_result_2 = 1000
    # end

    # for (varname, varjump) in variables_jump
    #   if getcategory(varjump) == :Bin
    #     println(varname, " : ", getvalue(varjump), "arrondi : ",(getvalue(varjump) > 0.5 ? 1 : 0))
    #     JuMP.fix(varjump,(getvalue(varjump) > 0.5 ? 1 : 0))
    #   end
    # end
    #
    # mysolver3 = KnitroSolver(KTR_PARAM_OUTLEV=3,
    #                         KTR_PARAM_MAXIT=600,
    #                         KTR_PARAM_SCALE=0,
    #                         KTR_PARAM_FEASTOL=1e-6,
    #                         KTR_PARAM_OPTTOL=1e-3,
    #                         KTR_PARAM_FEASTOLABS=1.0,
    #                         KTR_PARAM_OPTTOLABS=1.0,
    #                         KTR_PARAM_BAR_INITPT=2,
    #                         KTR_PARAM_PRESOLVE=0,
    #                         KTR_PARAM_HONORBNDS=0)
    # setsolver(m, mysolver3)
    # solve(m)


    ##get values
    println("Objective value : ", getobjectivevalue(m),"\n")

    # println("----Solution csv writing")
    # f = open(joinpath(outpath,"JuMP_solution.csv"),"w")
    # write(f, "Varname ; Value\n")
    # for (varname, var) in variables_jump
    #   value = getvalue(var)
    #   write(f, "$varname; $value\n")
    # end
    # close(f)

    minslack = +Inf
    ctr_minslack = ""
    # f = open(joinpath(outpath,"KnitroJuMP_constraints_eval.csv"),"w")
    # write(f, "Ctrname ; Value ; LB ; UB\n")
    for (ctrname, (exp,lb,ub)) in ctr_exp
      body_value = getvalue(exp)
      slack = min(body_value-lb, ub-body_value)
      if slack <= minslack
        minslack = slack
        ctr_minslack = ctrname
      end
      # write(f, "$ctrname; $body_value ; $lb ; $ub ; $slack\n")
    end
    # close(f)
    # println("----End solution csv writing\n")
    feas = minslack
    ctr = ctr_minslack
    println("min slack for constraints in problem JuMP Phase 2 : ($minslack,$ctr_minslack)")

    ##create solution1.txt and solution2.txt
    println("--Solution txt writing")
    write_solutions(OPFpbs, variables_jump, outpath)
    println("--End solution txt writing\n")

    println("--Reading solution1.txt and solution2.txt")
    sol_txt = read_solution_point_GOC(instance_path, outpath)
    pt_txt = cplx2real(sol_txt)
    min_slack = get_minslack(pb_global_real, pt_txt)
    infeas_ctr_txt = get_nb_infeasible_ctr_by_ctrtype(pb_global_real, pt_txt, 1e-6)
    println("get_minslack point from txt files:", min_slack)
    close(outlog)
    redirect_stdout(originalSTDOUT)
    return scenario => (status1, status2, (feas,ctr), min_slack, infeas_ctr_txt)
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
results = pmap(solve_GOC_via_Julia, [data_path for i=1:length(scenarios)], [folder for i=1:length(scenarios)], scenarios)

println("----------> para jobs done\n")
println("######################################################################################")



date = Dates.format(now(), "yy_u_dd_HH_MM_SS")
filename = joinpath("..","JuMP_runs","results_$(folder)_$(date).csv")
touch(filename)

f = open(filename, "w")

write(f, "Scenario;solve_result_1; solve_result_2;max relative slack from knitro point; ctr associated; max relative slack from txt files; ctr associated; opterror1 ; opterror2\n")
filename = joinpath("..","JuMP_runs","infeas_by_ctr_$(folder)_$(date).csv")
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
    solve_result_2 = data[2]
    feas,ctr1 = data[3]
    min_slack,ctr2 = data[4]
    infeas_ctr_txt = data[5]

    if solve_result_1!=0 || solve_result_2!=0 || feas > 1e-6 || min_slack > 1e-6
        nb_scenarios_with_pb +=1
        println("PB: $scenario not feasible")
    end
    write(f, "$(scenario);$(solve_result_1);$(solve_result_2);$(feas);$(ctr1);$(min_slack);$(ctr2)\n")
    write(f2, "$(scenario)")

    for ctr in ctrtypes
        if !haskey(infeas_ctr_txt, ctr)
            value = 0
        else
            value = infeas_ctr_txt[ctr]
        end
        write(f2, ";$value")
    end
    write(f2, "\n")
end
close(f2)
if nb_scenarios_with_pb > 0
    println("\nNB OF SCENARIOS NOT FEASIBLE : $nb_scenarios_with_pb/$nb_scenarios. \nSee results_$folder.csv in JuMP_runs for more details")
    write(f, "\n ; NB OF SCENARIOS NOT FEASIBLE;$nb_scenarios_with_pb/$nb_scenarios\n")
else
    println("ALL SCENARIOS ($nb_scenarios scenarios) FEASIBLE.\nSee results_$folder.csv in JuMP_runs for more details.")
    write(f, "\n ; NB OF SCENARIOS NOT FEASIBLE;0\n")
end
close(f)
