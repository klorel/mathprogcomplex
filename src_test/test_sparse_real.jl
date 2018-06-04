ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))

function print_primobj(io::IO, primobjectives, sols)
    @printf(io, "%15s  %6s  %15s  %15s  %15s  | %15s  %20s  %15s\n", "q5min", "d", "objective", "CJ obj", "Δobjs", "elapsed time", "tot. bytes alloc (MB)", "gctime")
    for ((q5min, d), (cur_obj, t, bytes, gctime)) in primobjectives

        target_obj = sols[q5min] # Rank relaxation value

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
    testfolder = joinpath("Mosek_runs", branch, "testWB5sparse", date)
    ispath(testfolder) && rm(testfolder, recursive=true)

    # Testing all rank relaxations
    sols = SortedDict(-30.80 => 945.83,
                      -20.51 => 954.82,
                      -10.22 => 963.83,
                      00.07 => 972.85,
                      10.36 => 981.89,
                      20.65 => 990.95,
                      30.94 => 1005.13,
                      41.23 => 1033.07,
                      51.52 => 1070.39,
                      61.81 => 1114.90)

    primobjectives_noeqs = SortedDict()
    for (q5min, vals) in sols

        problem = buildPOP_WB5(q5min=q5min, rmeqs=true)
        d = 2

        info("Working on WB5, no eqs, q5min=$q5min, d=$d")
        logpath = joinpath(testfolder, "WB5_q5min_$(q5min)_d_$(d)_noeq")
        ispath(logpath) && rm(logpath, recursive=true); mkpath(logpath)
        println("Saving file at $logpath")

        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                # symmetries=[PhaseInvariance],
                                issparse=true,
                                d = 1)

        max_cliques = get_WB5cliques(relax_ctx, problem)

        println("-------> $(now())")
        (primobj, dualobj), t, bytes, gctime, memallocs = @timed run_hierarchy(problem, relax_ctx, logpath, max_cliques=max_cliques);
        primobjectives_noeqs[(q5min, d)] = (primobj, t, bytes / (10^6), gctime)

        # Saving max_cliques
        open(joinpath(logpath, "maxcliques_relaxctx.txt"), "w") do fcliques
            println(fcliques, "max_cliques are:")
            println(fcliques, max_cliques)
            println(fcliques, "relaxation_ctx is:")
            println(fcliques, relax_ctx)
        end
    end

    primobjectives_eqs = SortedDict()

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
end

main()
