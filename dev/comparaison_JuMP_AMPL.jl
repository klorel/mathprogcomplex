include(joinpath(pwd(),"..","src_PowSysMod", "PowSysMod_body.jl"))
include("polyproblem_to_jump.jl")
instances = ["case9.m","case30.m","case118.m","case300.m","case1354pegase.m", "case2869pegase.m","case13659pegase.m"]

for instance in instances
  println(instance)
  instance_path = joinpath(pwd(),"..", "data_Matpower", "matpower",instance)

  ##read and load files
  OPFpbs = load_OPFproblems(MatpowerInput, instance_path)
  ## Bulding optimization problem
  pb_global = build_globalpb!(OPFpbs)

  tic()
  pb_global_real = pb_cplx2real(pb_global)
  t_cplx2real = toq()

  println("t_cplx2real : ", t_cplx2real)

  ## Exporting real problem
  tic()
  ## initial point
  init_point = Point()
  init_point_real = cplx2real(init_point)
  instance_name  = String(split(instance_path,'\\')[end])
  amplexportpath = joinpath("..","knitro_runs", "$(instance_name[1:end-2])")

  my_timer = @elapsed export_to_dat(pb_global_real, amplexportpath, init_point_real)
  @printf("%-35s%10.6f s\n", "export_to_dat", my_timer)
  t_buildexport = toq()

  println("time build export : ", t_buildexport)

  tic()
  mysolver = KnitroSolver(KTR_PARAM_OUTLEV=3,
                          KTR_PARAM_MAXIT=600,
                          KTR_PARAM_SCALE=0,
                          KTR_PARAM_FEASTOL=1.0,
                          KTR_PARAM_OPTTOL=1.0,
                          KTR_PARAM_FEASTOLABS=1e-6,
                          KTR_PARAM_OPTTOLABS=1e-3,
                          KTR_PARAM_BAR_INITPT=2,
                          KTR_PARAM_PRESOLVE=0,
                          KTR_PARAM_HONORBNDS=0)

  # my_timer = @elapsed (m, variables_jump) = get_JuMP_cartesian_model(pb_global, mysolver)
  # @printf("%-35s%10.6f s\n", "get_JuMP_cartesian_model", my_timer)
  (m, variables_jump) = get_JuMP_cartesian_model(pb_global_real, mysolver)
  t_buildmodel = toq()

  println("t_buildmodel : ", t_buildmodel)

end
