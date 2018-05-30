ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))

function main()
    repo = LibGit2.GitRepo(pwd()); branch = LibGit2.shortname(LibGit2.head(repo))
    date = String(Dates.format(now(), "mm_dd-HHhMM"))
    testfolder = joinpath("Mosek_runs", branch, "testWB2simple", date)
    ispath(testfolder) && rm(testfolder, recursive=true)

    params = OrderedSet([(1.022, 1, false, 861.51),
                        (1.022, 2, false, 901.38),
                        (1.022, 3, false, 905.73)])
                        # (1.022, 1, true,  861.51),
                        # (1.022, 2, true,  901.38),
                        # (1.022, 3, true,  905.73)])


    outstream = @sprintf("-> relaxations of | unchanged WB2 OPF | no equality contraints | ball constraint\n")
    outstream = outstream * @sprintf("%15s  %6s  %15s  %15s  %15s  | %15s  %20s  %15s\n", "v2max", "d", "objective", "CJ obj", "Δ", "elapsed time", "tot. bytes alloc (MB)", "gctime")

    ## Rank relaxation of unchanged WB2 (+ ball constraint)
    for (v2max, d, phaseset, CJobj) in params
        problem = buildPOP_WB2(v2max=v2max, rmeqs=true, setnetworkphase=phaseset)

        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            # symmetries=[PhaseInvariance],
                                            d = d)

        logpath = joinpath(testfolder, "v2max_$(v2max)_d_$(d)_"*(phaseset?"phaseset":"phasefree"))
        !ispath(logpath) || rm(logpath, recursive=true)
        mkpath(logpath)

        (cur_obj, dualobj), t, bytes, gctime, memallocs = @timed run_hierarchy(problem, relax_ctx, logpath, save_mmtpb=true)

        outstream = outstream * @sprintf("%15f  %6i  %15f  %15f  %15f  | %15f  %20f  %15f\n", v2max, d, cur_obj, CJobj, abs(cur_obj-CJobj), t, bytes/10^6, gctime)
    end

    outstream = outstream * @sprintf("\n-> relaxations of | fixed phase WB2 OPF | no equality contraints | ball constraint\n")
    outstream = outstream * @sprintf("%15s  %6s  %15s  %15s  %15s  | %15s  %20s  %15s\n", "v2max", "d", "objective", "CJ obj", "Δ", "elapsed time", "tot. bytes alloc (MB)", "gctime")

    ## Rank relaxation of fixed phase WB2 (+ ball constraint)
    for (v2max, d, phaseset, CJobj) in params
        phaseset=true
        problem = buildPOP_WB2(v2max=v2max, rmeqs=true, setnetworkphase=phaseset)

        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            # symmetries=[PhaseInvariance],
                                            d = d)

        logpath = joinpath(testfolder, "v2max_$(v2max)_d_$(d)_"*(phaseset?"phaseset":"phasefree"))
        !ispath(logpath) || rm(logpath, recursive=true)
        mkpath(logpath)

        (cur_obj, dualobj), t, bytes, gctime, memallocs = @timed run_hierarchy(problem, relax_ctx, logpath, save_mmtpb=true)

        outstream = outstream * @sprintf("%15f  %6i  %15f  %15f  %15f  | %15f  %20f  %15f\n", v2max, d, cur_obj, CJobj, abs(cur_obj-CJobj), t, bytes/10^6, gctime)
    end

    outstream = outstream * @sprintf("\n-> relaxations of | unchanged WB2 OPF | equality contraints | ball constraint\n")
    outstream = outstream * @sprintf("%15s  %6s  %15s  %15s  %15s  | %15s  %20s  %15s\n", "v2max", "d", "objective", "CJ obj", "Δ", "elapsed time", "tot. bytes alloc (MB)", "gctime")

    ## Rank relaxation of unchanged WB2 (+ ball constraint)
    for (v2max, d, phaseset, CJobj) in params
        problem = buildPOP_WB2(v2max=v2max, rmeqs=false, setnetworkphase=phaseset)

        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            # symmetries=[PhaseInvariance],
                                            d = d)

        logpath = joinpath(testfolder, "v2max_$(v2max)_d_$(d)_"*(phaseset?"phaseset":"phasefree"))
        !ispath(logpath) || rm(logpath, recursive=true)
        mkpath(logpath)

        (cur_obj, dualobj), t, bytes, gctime, memallocs = @timed run_hierarchy(problem, relax_ctx, logpath, save_mmtpb=true)

        outstream = outstream * @sprintf("%15f  %6i  %15f  %15f  %15f  | %15f  %20f  %15f\n", v2max, d, cur_obj, CJobj, abs(cur_obj-CJobj), t, bytes/10^6, gctime)
    end

    outstream = outstream * @sprintf("\n-> relaxations of | fixed phase WB2 OPF | equality contraints | ball constraint\n")
    outstream = outstream * @sprintf("%15s  %6s  %15s  %15s  %15s  | %15s  %20s  %15s\n", "v2max", "d", "objective", "CJ obj", "Δ", "elapsed time", "tot. bytes alloc (MB)", "gctime")

    ## Rank relaxation of fixed phase WB2 (+ ball constraint)
    for (v2max, d, phaseset, CJobj) in params
        phaseset=true
        problem = buildPOP_WB2(v2max=v2max, rmeqs=false, setnetworkphase=phaseset)

        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            # symmetries=[PhaseInvariance],
                                            d = d)

        logpath = joinpath(testfolder, "v2max_$(v2max)_d_$(d)_"*(phaseset?"phaseset":"phasefree"))
        !ispath(logpath) || rm(logpath, recursive=true)
        mkpath(logpath)

        (cur_obj, dualobj), t, bytes, gctime, memallocs = @timed run_hierarchy(problem, relax_ctx, logpath, save_mmtpb=true)

        outstream = outstream * @sprintf("%15f  %6i  %15f  %15f  %15f  | %15f  %20f  %15f\n", v2max, d, cur_obj, CJobj, abs(cur_obj-CJobj), t, bytes/10^6, gctime)
    end

    print(outstream)
    open(joinpath(testfolder, "summary.log"), "w") do fout
        print(fout, outstream)
    end
end

main()