include(joinpath(pwd(),"..","src_PowSysMod", "PowSysMod_body.jl"))
include("polyproblem_to_jump.jl")

function MyJulia1(rawfile,genfile,contfile)
  ##read and load files
  OPFpbs = load_OPFproblems(rawfile, genfile, contfile)
  introduce_Sgenvariables!(OPFpbs)
  ## Bulding optimization problem
  pb_global = build_globalpb!(OPFpbs)

  ##convert to JuMP model

  mysolver = KnitroSolver(KTR_PARAM_OUTLEV=3,
                          KTR_PARAM_MAXIT=600,
                          KTR_PARAM_SCALE=0,
                          KTR_PARAM_HESSOPT=0,
                          KTR_PARAM_FEASTOL=1.0,
                          KTR_PARAM_OPTTOL=1.0,
                          KTR_PARAM_FEASTOLABS=1e-6,
                          KTR_PARAM_OPTTOLABS=1e-3,
                          KTR_PARAM_BAR_INITPT=2,
                          KTR_PARAM_PRESOLVE=0,
                          KTR_PARAM_HONORBNDS=0,
                          KTR_PARAM_MIP_INTVAR_STRATEGY=0) 
  tic()
  m = get_JuMP_cartesian_model(pb_global, mysolver)
  toc()

  solve(m)


  ##create solution1.txt and solution2.txt
end

raw = "powersystem.raw"
gen = "generator.csv"
con = "contingency.csv"

instance_path = joinpath(pwd(),"..", "data_GOC", "Phase_0_IEEE14_1Scenario","scenario_1")
# instance_path = joinpath(pwd(), "data_GOC", "Phase_0_RTS96","scenario_1")
# instance_path = joinpath(pwd(), "data_GOC", "Phase_0_Feas179","scenario_1")


rawfile = joinpath(instance_path,raw)
genfile = joinpath(instance_path, gen)
contfile = joinpath(instance_path, con)

MyJulia1(rawfile,genfile,contfile)
