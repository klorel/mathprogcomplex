include(joinpath(pwd(),"src_PowSysMod", "PowSysMod_body.jl"))
# push!(LOAD_PATH, "D:\\repo\\GOC_Gimaju\\src_PowSysMod")
# using PowSysMod, JuMP, KNITRO

function MyJulia1(rawFile, genFile, contFile)
  println(rawFile, genFile, contFile)
  #
  # rGOC = "/state/partition1/data/141982//Phase_0_IEEE14_1Scenario/scenario_1/powersystem.raw"
  # gGOC = "/state/partition1/data/141982//Phase_0_IEEE14_1Scenario/scenario_1/generator.csv"
  # cGOC = "/state/partition1/data/141982//Phase_0_IEEE14_1Scenario/scenario_1/contingency.csv"
  #
  # # if rawFile == rGOC && genFile == gGOC && contFile == cGOC
  #   f = open("solution1_AMPL.jl","r")
  #   lines1 = readlines(f)
  #   close(f)
  #
  #   s1 = open("solution1.txt","w")
  #   for line in lines1
  #     write(s1, line,"\n")
  #   end
  #   close(s1)
  #
  #   f = open("solution2_AMPL.jl","r")
  #   lines2 = readlines(f)
  #   close(f)
  #
  #   s2 = open("solution2.txt","w")
  #   for line in lines2
  #     write(s2, line, "\n")
  #   end
  #   close(s2)
  # # end

  folder_path, scenario = splitdir(splitdir(rawFile)[1])
  folder = splitdir(folder_path)[2]
  outpath = joinpath("JuMP_runs","$(folder)_$(scenario)")
  ##read and load files
  OPFpbs = load_OPFproblems(rawFile, genFile, contFile)
  introduce_Sgenvariables!(OPFpbs)
  ## Building optimization problem
  pb_global = build_globalpb!(OPFpbs)

  ## convert into real problem
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
                          KTR_PARAM_FEASTOL=1.0,
                          KTR_PARAM_OPTTOL=1.0,
                          KTR_PARAM_FEASTOLABS=1.001e-6,
                          KTR_PARAM_OPTTOLABS=1e-3,
                          KTR_PARAM_BAR_INITPT=2,
                          KTR_PARAM_PRESOLVE=0,
                          KTR_PARAM_HONORBNDS=0,
                          KTR_PARAM_MIP_INTVAR_STRATEGY=1)
  tic()
  my_timer = @elapsed m, variables_jump, ctr_jump, ctr_exp = get_JuMP_cartesian_model(pb_global_real, mysolver)
  @printf("%-35s%10.6f s\n", "get_JuMP_cartesian_model", my_timer)
  toc()
  #resolution
  solve(m)

  minslack = +Inf
  ctr_minslack = ""
  for (ctrname, (exp,lb,ub)) in ctr_exp
    body_value = getvalue(exp)
    slack = min(body_value-lb, ub-body_value)
    if slack <= minslack
      minslack = slack
      ctr_minslack = ctrname
    end
  end
  println("min slack for constraints in problem JuMP Phase 1 : ($minslack,$ctr_minslack)")

  ##phase 2 : resolution with complementary constraints + initial point = solution continuous relaxation
  println("
##############################################################################################\n
Phase 2 : resolution with complementary constraints from the solution of continuous relaxation\n
##############################################################################################")

  mysolver2 = KnitroSolver(KTR_PARAM_OUTLEV=3,
                          KTR_PARAM_MAXIT=600,
                          KTR_PARAM_SCALE=0,
                          KTR_PARAM_FEASTOL=1.0,
                          KTR_PARAM_OPTTOL=1.0,
                          KTR_PARAM_FEASTOLABS=1.001e-6,
                          KTR_PARAM_OPTTOLABS=1e-3,
                          KTR_PARAM_BAR_INITPT=2,
                          KTR_PARAM_PRESOLVE=0,
                          KTR_PARAM_HONORBNDS=0,
                          KTR_PARAM_MIP_INTVAR_STRATEGY=2)


  setsolver(m, mysolver2)
  # tic()
  # my_timer = @elapsed m2, variables_jump2 = get_JuMP_cartesian_model(pb_global_real, mysolver)
  # @printf("%-35s%10.6f s\n", "get_JuMP_cartesian_model", my_timer)
  # f = open(joinpath(outpath,"JuMP_solution_phase1.csv"),"w")
  # write(f, "Varname ; Value\n")
  # for (varname, varjump) in variables_jump
  #   write(f, "$varname; $(getvalue(varjump))\n")
  #   setvalue(variables_jump[varname], getvalue(varjump))
  # end
  # close(f)
  # toc()
  #resolution
  solve(m)

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
  println("min slack for constraints in problem JuMP Phase 2 : ($minslack,$ctr_minslack)")

  ##create solution1.txt and solution2.txt
  println("--Solution txt writing")
  write_solutions(OPFpbs, variables_jump, pwd())
  # write_solutions(OPFpbs, variables_jump, outpath)
  println("--End solution txt writing\n")

  return pb_global_real

end
