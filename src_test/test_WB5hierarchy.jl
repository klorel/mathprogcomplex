ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))

function print_primobj(io::IO, primobjectives, sols)
    @printf(io, "%15s  %6s  %15s  %15s  %15s  | %15s  %20s  %15s\n", "q5min", "d", "objective", "CJ obj", "Δobjs", "elapsed time", "tot. bytes alloc (MB)", "gctime")
    for ((q5min, d), (cur_obj, t, bytes, gctime)) in primobjectives
        target_obj = -1
        if d == 2
            target_obj = sols[q5min][3] # Rank relaxation value
        elseif d == 2*sols[q5min][1]
            target_obj = sols[q5min][2] # Optimal value
        end

        CJ_obj = target_obj!=-1?target_obj:NaN
        Δobj = target_obj!=-1?abs(cur_obj-target_obj):NaN
        @printf(io, "%15f  %6i  %15f  %15f  %15f  | %15f  %20f  %15f\n", q5min, d/2, cur_obj, CJ_obj, Δobj, t, bytes, gctime)
        if target_obj!=-1 && !isapprox(cur_obj, target_obj, atol=1e-2)
            warn("Test WB5: (q5min, d) = ($q5min, $d) found objective is $cur_obj, expected $target_obj.")
        end
    end
end

function main()
    repo = LibGit2.GitRepo(pwd()); branch = LibGit2.shortname(LibGit2.head(repo))
    date = String(Dates.format(now(), "mm_dd-HHhMM"))
    testfolder = joinpath("Mosek_runs", branch, "testWB5extensive", date)
    ispath(testfolder) && rm(testfolder, recursive=true)

    sols = SortedDict(-30.80 => (2, 945.83, 945.83),
                       -20.51 => (2, 1146.48, 954.82),
                       -10.22 => (2, 1209.11, 963.83),
                        00.07 => (2, 1267.79, 972.85),
                        10.36 => (2, 1323.86, 981.89),
                        20.65 => (2, 1377.97, 990.95),
                        30.94 => (2, 1430.54, 1005.13),
                        41.23 => (2, 1481.81, 1033.07),
                        51.52 => (2, 1531.97, 1070.39),
                        61.81 => (1, 1114.90, 1114.90))

    primobjectives_noeqs = SortedDict()
    for (q5min, vals) in sols
        dmax = vals[1]

        problem = buildPOP_WB5(q5min=q5min, rmeqs=true)
        dmax = 1 ## To Be Removed
        for d=2:2:2*dmax  # Complex d definition, twice real one...
            info("Working on WB5, no eqs, q5min=$q5min, d=$d")
            logpath = joinpath(testfolder, "WB5_v2max_$(q5min)_d_$(d)_noeq")
            ispath(logpath) && rm(logpath, recursive=true); mkpath(logpath)
            println("Saving file at $logpath")

            relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                                d = d,
                                                symmetries = [PhaseInvariance])

            println("-------> $(now())")
            (primobj, dualobj), t, bytes, gctime, memallocs = @timed run_hierarchy(problem, relax_ctx, logpath);
            primobjectives_noeqs[(q5min, d)] = (primobj, t, bytes / (10^6), gctime)
        end
    end

    primobjectives_eqs = SortedDict()
    for (q5min, vals) in sols
        dmax = vals[1]

        problem = buildPOP_WB5(q5min=q5min, rmeqs=false)
        dmax = 1 ## To Be Removed
        for d=2:2:2*dmax  # Complex d definition, twice real one...
            info("Working on WB5, eqs, q5min=$q5min, d=$d")
            logpath = joinpath(testfolder, "WB5_v2max_$(q5min)_d_$(d)_eq")
            ispath(logpath) && rm(logpath, recursive=true); mkpath(logpath)
            println("Saving file at $logpath")

            relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                                d = d,
                                                symmetries = [PhaseInvariance])

            println("-------> $(now())")
            (primobj, dualobj), t, bytes, gctime, memallocs = @timed run_hierarchy(problem, relax_ctx, logpath);
            primobjectives_eqs[(q5min, d)] = (primobj, t, bytes / (10^6), gctime)
        end
    end

    println("- No equality constraints:")
    print_primobj(STDOUT, primobjectives_noeqs, sols)

    println("- Equality constraints:")
    print_primobj(STDOUT, primobjectives_eqs, sols)

    logfile = joinpath(testfolder, "test_WB5.txt")

    isfile(logfile) && rm(logfile)

    info("Writing results to $logfile")
    open(logfile, "w") do fout
        println(fout, "- No equality constraints:")
        print_primobj(fout, primobjectives_noeqs, sols)

        println(fout, "- Equality constraints:")
        print_primobj(fout, primobjectives_eqs, sols)
    end

    return primobjectives_noeqs, primobjectives_eqs
end

main()