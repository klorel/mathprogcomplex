ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))
include(joinpath(ROOT, "src_PowSysMod", "PowSysMod_body.jl"))

using BenchmarkTools

OPFpbs = load_OPFproblems(MatpowerInput, joinpath("..", "data", "data_Matpower", "matpower", "case89pegase.m"))
# OPFpbs = load_OPFproblems(MatpowerInput, joinpath("..", "data", "data_Matpower", "matpower", "case30.m"))
const problem_c = build_globalpb!(OPFpbs)
const problem = pb_cplx2real(problem_c)

d = 1

relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                # symmetries=[PhaseInvariance],
                                d = d)

max_cliques = get_maxcliques(relax_ctx, problem)

momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem, max_cliques)

mmtrel_pb = MomentRelaxation{Float64}(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)

@benchmark (sdpinstance = build_SDPInstance(relax_ctx, mmtrel_pb))
