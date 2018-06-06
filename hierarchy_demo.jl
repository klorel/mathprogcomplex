ROOT = pwd()
include(joinpath("src_SOShierarchy", "SOShierarchy.jl"))


## GLOBAL OPTIMIZATION WITH POLYNOMIALS AND THE PROBLEM OF MOMENTS
# JEAN B. LASSERRE, 2001
# http://www.ii.uib.no/~lennart/drgrad/Lasserre2001.pdf

function main()
    problem = Problem()
    x1 = Variable("x1", Real); add_variable!(problem, x1)
    x2 = Variable("x2", Real); add_variable!(problem, x2)

    set_objective!(problem, -(x1-1)^2 -(x1-x2)^2 -(x2-3)^2)
    add_constraint!(problem, "crt1", (1-(x1-1)^2) >> 0)
    add_constraint!(problem, "crt2", (1-(x1-x2)^2) >> 0)
    add_constraint!(problem, "crt3", (1-(x2-3)^2) >> 0)

    workdir = joinpath("Mosek_runs", "testfolder")
    ispath(workdir) && rm(workdir, recursive=true); mkpath(workdir)

    println(problem)
    readline()

    for d=1:2
        logpath = joinpath(workdir, "exmaple2_d_$d"); mkpath(logpath)

        info("Solving real problem relaxation at order $d")

        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            d = d)

        println("\nRelaxation parameters:")
        println(relax_ctx)
        readline()

        max_cliques = get_maxcliques(relax_ctx, problem)
        println("Maximal cliques:")
        println(max_cliques)
        readline()

        ########################################
        # Compute moment matrices parameters: order et variables
        momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem, max_cliques)

        ########################################
        # Compute moment and localization matrices
        moment_relaxation_pb = MomentRelaxation(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)
        println("Moment relaxation problem:")
        println(moment_relaxation_pb)
        readline()

        ########################################
        # Convert to a primal SDP problem
        SOS_pb = build_SDPInstance(relax_ctx, moment_relaxation_pb)
        println("SOS problem:")
        println(SOS_pb)
        readline()

        println("Mosek solve:")
        mosekpb = build_mosekpb(SOS_pb, logpath)

        primal = SortedDict{Tuple{String,String,String}, Float64}()
        dual = SortedDict{Tuple{String, String, String}, Float64}()
        primobj, dualobj = solve_mosek(mosekpb::SDP_Problem, primal, dual; logname = joinpath(logpath, "Mosek_run.log"))

        readline()
    end
end

main()