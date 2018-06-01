ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))

using Base.Profile
using ProfileView

function toprofile(n, problem, d)
    # problem = buildPOP_WB2(v2max=1.022, setnetworkphase=false)
    # problem = buildPOP_WB5()

    for i=1:n
        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        # symmetries=[PhaseInvariance],
                                        d = d)

        sdpinstance = build_relaxation(problem, relax_ctx)
    end
end


function main()

    OPFpbs = load_OPFproblems(MatpowerInput, joinpath("..", "data", "data_Matpower", "matpower", "case30.m"))
    problem_c = build_globalpb!(OPFpbs)
    problem = pb_cplx2real(problem_c)

    toprofile(1, problem, 1)

    Profile.clear()
    @profile toprofile(3, problem, 1)
    ProfileView.view()
end

# main()
