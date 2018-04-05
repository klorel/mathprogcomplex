##Test connexion with JuMP


ROOT = pwd()
include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))
include("polyproblem_to_jump.jl")
# typeofinput, instance_path, eps = read_arguments(ARGS)

typeofinput = MatpowerInput
data_path = joinpath(ROOT,"..","data_Matpower","matpower")
instance = "case9.m"


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
                        KTR_PARAM_MIP_INTVAR_STRATEGY=0) #0 nothing, 1 relax, 2 complementary

instance_path = joinpath(data_path, instance)

OPFpbs = load_OPFproblems(typeofinput, instance_path)

## Introducing coupling constraints on generator output
(typeofinput != GOCInput) || introduce_Sgenvariables!(OPFpbs)

## Bulding optimization problem
pb_global = build_globalpb!(OPFpbs)

print(pb_global)

m = get_JuMP_cartesian_model(pb_global, mysolver)

print(m)


solve(m)
