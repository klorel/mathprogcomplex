ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))

using Base.Profile
using ProfileView

function toprofile(n, d)
    # problem = buildPOP_WB2(v2max=1.022, setnetworkphase=false)
    problem = buildPOP_WB5()

    for i=1:n
        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        # symmetries=[PhaseInvariance],
                                        d = d)

        sdpinstance = build_relaxation(problem, relax_ctx)
    end
end


function main()

    toprofile(1, 1)

    Profile.clear()
    @profile toprofile(10, 2)
    ProfileView.view()
end

# main()
