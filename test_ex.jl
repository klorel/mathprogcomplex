ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))


function main()

    x1 = Variable("x1", Real)
    x2 = Variable("x2", Real)
    x3 = Variable("x3", Real)

    problem = Problem()
    add_variable!(problem, x1)
    add_variable!(problem, x2)
    add_variable!(problem, x3)
    set_objective!(problem, -2*x1+x2-x3)
    add_constraint!(problem, "ctr1", (x1*(4*x1-4*x2+4*x3-20)+x2*(2*x2-2*x3+9)+x3*(2*x3-13)+24) >> 0)
    add_constraint!(problem, "ctr2", (x1+x2+x3) << 4)
    add_constraint!(problem, "ctr3", (3*x2+x3) << 6)
    add_constraint!(problem, "def_x1", 0 << x1 << 2)
    add_constraint!(problem, "def_x2", 0 << x2)
    add_constraint!(problem, "def_x3", 0 << x3 << 3)

    println("Problem is:\n$problem")

    data = SortedDict()
    for d=1:4
        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            d = d)

        logpath = joinpath("Mosek_runs", "Lasserreex", "d_$d")
        ispath(logpath) && rm(logpath, recursive=true)
        mkpath(logpath)
        cur_obj, dualobj = run_hierarchy(problem, relax_ctx, logpath, indentedprint=true, save_pbs=true)

        @show d, cur_obj, dualobj
        data[d] = (cur_obj, dualobj)
        # run_hierarchy(problem::Problem, relax_ctx::RelaxationContext, logpath; indentedprint=false,
        #                                                                         max_cliques::SortedDict{String, SortedSet{Variable}}=SortedDict{String, SortedSet{Variable}}(),
        #                                                                         save_pbs=false)
    end

    for (d, (cur_obj, dualobj)) in data
        println("$d  -> $cur_obj, $dualobj")
    end

end

main()