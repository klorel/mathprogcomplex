ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))

function main()

    params = OrderedSet([(1.028, 2, true, 885.71),
                         (1.028, 4, true, NaN),
                         (1.028, 6, true, 905.73),
                         (1.028, 2, false, 885.71),
                         (1.028, 4, false, NaN),
                         (1.028, 6, false, 905.73)])


    outstream = @sprintf("%15s  %6s  %15s  %15s  %15s  | %15s  %20s  %15s\n", "v2max", "d", "objective", "CJ obj", "Δ", "elapsed time", "tot. bytes alloc (MB)", "gctime")

    ## Rank relaxation
    for (v2max, d, rmeqs, CJobj) in params
        problem = buildPOP_WB2(v2max=1.028, rmeqs=rmeqs)

        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        d = d)

        logpath = joinpath("Mosek_runs", "test_WB2simple", "v2max_$(v2max)_d_$(d)_"*(rmeqs?"noeqs":"eqs"))
        !ispath(logpath) || rm(logpath, recursive=true)
        mkpath(logpath)

        (cur_obj, dualobj), t, bytes, gctime, memallocs = @timed run_hierarchy(problem, relax_ctx, logpath)

        outstream = outstream * @sprintf("%15f  %6i  %15f  %15f  %15f  | %15f  %20f  %15f\n", v2max, d, cur_obj, CJobj, abs(cur_obj-CJobj), t, bytes/10^6, gctime)
    end
    print(outstream)
end

main()