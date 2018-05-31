ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))

using Base.Profile
using ProfileView

function toprofile(n)
    for i=1:n
        problem = buildPOP_WB2(v2max=1.022, setnetworkphase=false)
        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        # symmetries=[PhaseInvariance],
                                        d = 1)

        sdpinstance = build_relaxation(problem, relax_ctx)
    end
end


function main()

    toprofile(1)

    Profile.clear()
    @profile toprofile(100)
    ProfileView.view()
end

main()