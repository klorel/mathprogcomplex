ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))

using Base.Profile
using ProfileView
using BenchmarkTools

function toprofile(n, problem, d)
    # problem = buildPOP_WB2(v2max=1.022, setnetworkphase=false)
    # problem = buildPOP_WB5()

    for i=1:n
        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        # symmetries=[PhaseInvariance],
                                        d = d)

        max_cliques = get_maxcliques(relax_ctx, problem)

        momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem, max_cliques)

        mmtrel_pb = MomentRelaxation{Float64}(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)

        # sdpinstance = build_SOSrelaxation(relax_ctx, mmtrel_pb)

    end
end


OPFpbs = load_OPFproblems(MatpowerInput, joinpath("..", "data", "data_Matpower", "matpower", "case89pegase.m"))
const problem_c = build_globalpb!(OPFpbs)
const problem = pb_cplx2real(problem_c)

d = 1

toprofile(1, problem, d)

Profile.clear()
@profile toprofile(10, problem, d)
ProfileView.view()


## Benchmark
relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                # symmetries=[PhaseInvariance],
                                d = d)

max_cliques = get_maxcliques(relax_ctx, problem)

momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem, max_cliques)

@benchmark (mmtrel_pb = MomentRelaxation{Float64}(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques))