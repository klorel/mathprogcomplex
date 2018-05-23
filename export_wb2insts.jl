ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))


function main()

    problem = buildPOP_WB2(rmeqs=false) # v2max=0.976


    ## Rank relaxation
    logpath = joinpath(pwd(), "WB2_ordre1")
    mkpath(logpath)

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        d = 2)

    primobj, dualobj = run_hierarchy(problem, relax_ctx, logpath, indentedprint=true)


    ## Order 2 relaxation
    logpath = joinpath(pwd(), "WB2_ordre2")
    mkpath(logpath)

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                    d = 4)

    primobj, dualobj = run_hierarchy(problem, relax_ctx, logpath)

end

main()
