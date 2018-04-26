#Test connexion with JuMP
ROOT = pwd()
include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))


typeofinput = MatpowerInput
data_path = joinpath(ROOT,"..","..","data","data_Matpower","matpower")
instance = "case9.m"


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
                        KTR_PARAM_MIP_INTVAR_STRATEGY=0) #0 nothing, 1 relax, 2 complementary

instance_path = joinpath(data_path, instance)
OPFpbs = load_OPFproblems(typeofinput, instance_path)
## Bulding optimization problem
pb_global = build_globalpb!(OPFpbs)
## convert into real problem
pb_global_real = pb_cplx2real(pb_global)

m, variables_jump, ctr_jump, ctr_exp = get_JuMP_cartesian_model(pb_global_real, mysolver)

solve(m)
