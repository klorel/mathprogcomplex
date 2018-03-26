if length(ARGS)==0
  error("Argument must be a matpower instance as case9.m for example")
end

instance = ARGS[1]

instance_path = joinpath(pwd(),"..", "data_Matpower", "matpower",instance)

include(joinpath(pwd(),"..","src_PowSysMod", "PowSysMod_body.jl"))
include("polyproblem_to_jump.jl")

function MyJulia_matpower(instance_path)
  ##read and load files
  OPFpbs = load_OPFproblems(MatpowerInput, instance_path)
  ## Bulding optimization problem
  pb_global = build_globalpb!(OPFpbs)

  ##convert to JuMP model
  ##convert to real variables

  tic()
  pb_poly_real = pb_cplx2real(pb_global)

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

  println(mysolver)

  # my_timer = @elapsed (m, variables_jump) = get_JuMP_cartesian_model(pb_global, mysolver)
  # @printf("%-35s%10.6f s\n", "get_JuMP_cartesian_model", my_timer)
  (m, variables_jump) = get_JuMP_cartesian_model(pb_poly_real, mysolver)
  t_buildmodel = toq()

  println("t_buildmodel : ", t_buildmodel)



  solve(m)

  ##create solution1.txt and solution2.txt
end


MyJulia_matpower(instance_path)
