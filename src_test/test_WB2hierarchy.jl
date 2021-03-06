ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))

function print_primobj(io::IO, primobjectives, sols)
    @printf(io, "%15s  %6s  %15s  %15s  %15s  | %15s  %20s  %15s\n", "v2max", "d", "objective", "CJ obj", "Δobjs", "elapsed time", "tot. bytes alloc (MB)", "gctime")
    for ((v2max, d), (cur_obj, t, bytes, gctime)) in primobjectives
        target_obj = -1
        if d == 1
            target_obj = sols[v2max][3] # Rank relaxation value
        elseif d == sols[v2max][1]
            target_obj = sols[v2max][2] # Optimal value
        end

        CJ_obj = target_obj!=-1?target_obj:NaN
        Δobj = target_obj!=-1?abs(cur_obj-target_obj):NaN
        @printf(io, "%15f  %6i  %15f  %15f  %15f  | %15f  %20f  %15f\n", v2max, d, cur_obj, CJ_obj, Δobj, t, bytes, gctime)
    end
end

function main()
    repo = LibGit2.GitRepo(pwd()); branch = LibGit2.shortname(LibGit2.head(repo))
    date = String(Dates.format(now(), "mm_dd-HHhMM"))
    testfolder = joinpath("Mosek_runs", branch, "testWB2extensive", date)
    ispath(testfolder) && rm(testfolder, recursive=true)

    sols = SortedDict(0.976 => (2, 905.76, 905.76),
                      0.983 => (2, 905.73, 903.12),
                      0.989 => (2, 905.73, 900.84),
                      0.996 => (2, 905.73, 898.17),
                      1.002 => (2, 905.73, 895.86),
                      1.009 => (2, 905.73, 893.16),
                      1.015 => (2, 905.73, 890.82),
                      1.022 => (3, 905.73, 888.08),
                      1.028 => (3, 905.73, 885.71),
                      1.035 => (2, 882.97, 882.97))

    primobjectives_noeqs = SortedDict()
    for (v2max, vals) in sols
        dmax = vals[1]

        ## Rank relaxation of WB2 vanilla
        problem = buildPOP_WB2(v2max=v2max, setnetworkphase=false)

        for d = 1:dmax
            info("Working on WB2, no eqs, v2max=$v2max, d=$d")
            logpath = joinpath(testfolder, "WB2_v2max_$(v2max)_d_$(d)_noeq")
            ispath(logpath) && rm(logpath, recursive=true); mkpath(logpath)
            println("Saving file at $logpath")

            relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                                symmetries=[PhaseInvariance],
                                                d = d)

            (primobj, dualobj), t, bytes, gctime, memallocs = @timed run_hierarchy(problem, relax_ctx, logpath, save_pbs=true);
            primobjectives_noeqs[(v2max, d)] = (primobj, t, bytes / 10^6, gctime)
        end
        # ## Minimal CV order of Phase-fixed OPF according to Josz
        # problem = buildPOP_WB2(v2max=v2max, rmeqs=true, setnetworkphase=true)
        # d = dmax

        # info("Working on WB2, no eqs, v2max=$v2max, d=$d")
        # logpath = joinpath(testfolder, "WB2_v2max_$(v2max)_d_$(d)_noeq")
        # ispath(logpath) && rm(logpath, recursive=true); mkpath(logpath)
        # println("Saving file at $logpath")

        # relax_ctx = set_relaxation(problem; hierarchykind=:Real,
        #                                     d = d)

        # (primobj, dualobj), t, bytes, gctime, memallocs = @timed run_hierarchy(problem, relax_ctx, logpath, save_pbs=true);
        # primobjectives_noeqs[(v2max, d)] = (primobj, t, bytes / 10^6, gctime)
    end

    # primobjectives_eqs = SortedDict()
    # for (v2max, vals) in sols
    #     dmax = vals[1]

    #     problem = buildPOP_WB2(v2max=v2max, rmeqs=false, setnetworkphase=true)

    #     for d=1:dmax  # Complex d definition, twice real one...
    #         info("Working on WB2, eqs, v2max=$v2max, d=$d")
    #         logpath = joinpath(testfolder, "WB2_v2max_$(v2max)_d_$(d)_eq")
    #         ispath(logpath) && rm(logpath, recursive=true); mkpath(logpath)
    #         println("Saving file at $logpath")

    #         relax_ctx = set_relaxation(problem; hierarchykind=:Real,
    #                                             d = d)

    #         (primobj, dualobj), t, bytes, gctime, memallocs = @timed run_hierarchy(problem, relax_ctx, logpath);
    #         primobjectives_eqs[(v2max, d)] = (primobj, t, bytes / 10^6, gctime)
    #     end
    # end

    println("- No equality constraints:")
    print_primobj(STDOUT, primobjectives_noeqs, sols)

    # println("- Equality constraints:")
    # print_primobj(STDOUT, primobjectives_eqs, sols)

    logfile = joinpath(testfolder, "test_WB2.txt")

    isfile(logfile) && rm(logfile)

    info("Writing results to $logfile")
    open(logfile, "w") do fout
        println(fout, "- No equality constraints:")
        print_primobj(fout, primobjectives_noeqs, sols)

        # println(fout, "- Equality constraints:")
        # print_primobj(fout, primobjectives_eqs, sols)
    end

    return primobjectives_noeqs#, primobjectives_eqs
end

main()