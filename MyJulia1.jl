include(joinpath(pwd(),"src_PowSysMod", "PowSysMod_body.jl"))
# push!(LOAD_PATH, "D:\\repo\\GOC_Gimaju\\src_PowSysMod")
# using PowSysMod, JuMP, KNITRO

function MyJulia1(rawFile, genFile, contFile)
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
                          KTR_PARAM_FEASTOLABS=1e-6,
                          KTR_PARAM_OPTTOLABS=1e-3,
                          KTR_PARAM_BAR_INITPT=2,
                          KTR_PARAM_PRESOLVE=0,
                          KTR_PARAM_HONORBNDS=0,
                          KTR_PARAM_MIP_INTVAR_STRATEGY=1)
  tic()
  my_timer = @elapsed m, variables_jump = get_JuMP_cartesian_model(pb_global_real, mysolver)
  @printf("%-35s%10.6f s\n", "get_JuMP_cartesian_model", my_timer)
  toc()
  #resolution
  solve(m)

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
                          KTR_PARAM_FEASTOLABS=1e-6,
                          KTR_PARAM_OPTTOLABS=1e-3,
                          KTR_PARAM_BAR_INITPT=2,
                          KTR_PARAM_PRESOLVE=0,
                          KTR_PARAM_HONORBNDS=0,
                          KTR_PARAM_MIP_INTVAR_STRATEGY=2)


  setsolver(m, mysolver2)
  # tic()
  # my_timer = @elapsed m2, variables_jump2 = get_JuMP_cartesian_model(pb_global_real, mysolver)
  # @printf("%-35s%10.6f s\n", "get_JuMP_cartesian_model", my_timer)
  for (varname, varjump) in variables_jump
    setvalue(variables_jump[varname], getvalue(varjump))
  end
  # toc()
  #resolution
  solve(m)

  ##get values
  println("Objective value : ", getobjectivevalue(m),"\n")

  println("----Solution csv writing")
  f = open("JuMP_solution.csv","w")
  write(f, "Varname ; Value\n")
  for (varname, var) in variables_jump
    value = getvalue(var)
    write(f, "$varname; $value\n")
  end
  close(f)
  println("----End solution csv writing\n")

  ##create solution1.txt and solution2.txt
  println("--Solution txt writing")
  write_solutions(OPFpbs, variables_jump)
  println("--End solution txt writing\n")
   return pb_global_real

end
