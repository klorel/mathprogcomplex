ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))

function main()
    v2max = 1.028
    problem = buildPOP_WB2(v2max=1.028, rmineqs=true)

    ## Rank relaxation
    d2=2  # Complex d definition, twice real one...
    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        d = d2)

    (cur_obj2, dualobj), t2, bytes2, gctime2, memallocs = @timed run_hierarchy(problem, relax_ctx)

    ## Exact relaxation
    d = 6
    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        d = d)

    (cur_obj, dualobj), t, bytes, gctime, memallocs = @timed run_hierarchy(problem, relax_ctx)
    info("Rank relaxation")
    @printf("%15s  %6s  %15s  %15s  | %15s  %20s  %15s\n", "v2max", "d", "objective", "CJ obj", "elapsed time", "tot. bytes alloc (MB)", "gctime")
    @printf("%15f  %6i  %15f  %15f  | %15f  %20f  %15f\n", v2max, d2, cur_obj2, 885.71, t2, bytes2/10^6, gctime2)
    @printf("%15f  %6i  %15f  %15f  | %15f  %20f  %15f\n", v2max, d, cur_obj, 905.73, t, bytes/10^6, gctime)
end

main()
