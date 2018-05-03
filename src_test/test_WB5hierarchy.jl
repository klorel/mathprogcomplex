ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))

function print_primobj(primobjectives, sols)
    @printf("%15s  %6s  %15s  %15s  | %15s  %20s  %15s\n", "q5min", "d", "objective", "CJ obj", "elapsed time", "tot. bytes alloc (MB)", "gctime")
    for ((q5min, d), (cur_obj, t, bytes, gctime)) in primobjectives
        target_obj = -1
        if d == 2
            target_obj = sols[q5min][3] # Rank relaxation value
        elseif d == 2*sols[q5min][1]
            target_obj = sols[q5min][2] # Optimal value
        end

        if target_obj!=-1 && !isapprox(cur_obj, target_obj, atol=1e-2)
            warn("Test WB5: (q5min, d) = ($q5min, $d) found objective is $cur_obj, expected $target_obj.")
        end

        @printf("%15f  %6i  %15f  %15f  | %15f  %20f  %15f\n", q5min, d, cur_obj, target_obj!=-1?abs(cur_obj-target_obj):NaN, t, bytes, gctime)
    end
end

function main()
    sols = SortedDict([-30.80 => (2, 945.83, 945.83),
                       -20.51 => (2, 1146.48, 954.82),
                       -10.22 => (2, 1209.11, 963.83),
                        00.07 => (2, 1267.79, 972.85),
                        10.36 => (2, 1323.86, 981.89),
                        20.65 => (2, 1377.97, 990.95),
                        30.94 => (2, 1430.54, 1005.13),
                        41.23 => (2, 1481.81, 1033.07),
                        51.52 => (2, 1531.97, 1070.39),
                        61.81 => (1, 1114.90, 1114.90)])

    primobjectives_noeqs = SortedDict()
    for (q5min, vals) in sols
        dmax = vals[1]
        problem = buildPOP_WB5(q5min=q5min, rmineqs=true)

        dmax = 1 ## To Be Removed
        for d=2:2:2*dmax  # Complex d definition, twice real one...
            relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                                d = d,
                                                symmetries = [PhaseInvariance])
            logpath = joinpath("Mosek_runs", "WB5_q5min_$(q5min)_d_$(d)")
            ispath(logpath) && rm(logpath, recursive=true); mkpath(logpath)

            println("-------> $(now())")
            (primobj, dualobj), t, bytes, gctime, memallocs = @timed run_hierarchy(problem, relax_ctx, logpath);
            primobjectives_noeqs[(q5min, d)] = (primobj, t, bytes / (10^6), gctime)
        end
    end

    primobjectives_eqs = SortedDict()
    for (q5min, vals) in sols
        dmax = vals[1]
        problem = buildPOP_WB5(q5min=q5min, rmineqs=false)

        dmax = 1 ## To Be Removed
        for d=2:2:2*dmax  # Complex d definition, twice real one...
            relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                                d = d,
                                                symmetries = [PhaseInvariance])
            logpath = joinpath("Mosek_runs", "WB5_q5min_$(q5min)_d_$(d)")
            ispath(logpath) && rm(logpath, recursive=true); mkpath(logpath)

            println("-------> $(now())")
            (primobj, dualobj), t, bytes, gctime, memallocs = @timed run_hierarchy(problem, relax_ctx, logpath);
            primobjectives_eqs[(q5min, d)] = (primobj, t, bytes / (10^6), gctime)
        end
    end

    println("- No equality constraints:")
    print_primobj(primobjectives_noeqs, sols)

    println("- Equality constraints:")
    print_primobj(primobjectives_eqs, sols)

    return primobjectives_noeqs, primobjectives_eqs
end

main()
