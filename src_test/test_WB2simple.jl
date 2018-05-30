ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))

function main()
    repo = LibGit2.GitRepo(pwd()); branch = LibGit2.shortname(LibGit2.head(repo))
    date = String(Dates.format(now(), "mm_dd-HHhMM"))
    testfolder = joinpath("Mosek_runs", branch, "testWB2simple", date)
    ispath(testfolder) && rm(testfolder, recursive=true)

    params = OrderedSet([(1.028, 1, true, 885.71),
                        (1.028, 2, true, NaN),
                        (1.028, 3, true, 905.73),
                        (1.028, 1, false, 885.71),
                        (1.028, 2, false, NaN),
                        (1.028, 3, false, 905.73)])


    outstream = @sprintf("%15s  %6s  %15s  %15s  %15s  | %15s  %20s  %15s\n", "v2max", "d", "objective", "CJ obj", "Δ", "elapsed time", "tot. bytes alloc (MB)", "gctime")

    ## Rank relaxation
    for (v2max, d, rmeqs, CJobj) in params
        problem = buildPOP_WB2(v2max=1.028, rmeqs=rmeqs, fixvariable=false)

        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            symmetries=[PhaseInvariance],
                                            d = d)

        logpath = joinpath(testfolder, "v2max_$(v2max)_d_$(d)_"*(rmeqs?"noeqs":"eqs"))
        !ispath(logpath) || rm(logpath, recursive=true)
        mkpath(logpath)

        (cur_obj, dualobj), t, bytes, gctime, memallocs = @timed run_hierarchy(problem, relax_ctx, logpath)

        outstream = outstream * @sprintf("%15f  %6i  %15f  %15f  %15f  | %15f  %20f  %15f\n", v2max, d, cur_obj, CJobj, abs(cur_obj-CJobj), t, bytes/10^6, gctime)
    end
    print(outstream)
end

main()