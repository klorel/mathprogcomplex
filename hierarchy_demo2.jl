ROOT = pwd()
include(joinpath("src_SOShierarchy", "SOShierarchy.jl"))


## GLOBAL OPTIMIZATION WITH POLYNOMIALS AND THE PROBLEM OF MOMENTS
# JEAN B. LASSERRE, 2001
# http://www.ii.uib.no/~lennart/drgrad/Lasserre2001.pdf

function main()

    workdir = joinpath("Mosek_runs", "testfolder")
    ispath(workdir) && rm(workdir, recursive=true); mkpath(workdir)

    OPFpbs = load_OPFproblems(MatpowerInput, joinpath("..", "data", "data_Matpower", "matpower", "WB2.m"))
    problem = build_globalpb!(OPFpbs)

    relax_ctx = set_relaxation(problem; hierarchykind=:Complex,
                                        # symmetries=[PhaseInvariance],
                                        d = 1)

    println(problem)
    readline()

    d = 1

    logpath = joinpath(workdir, "exmaple2_d_$d"); mkpath(logpath)

    info("Solving complex problem relaxation at order $d")

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

    println("Mosek solve: Coming soon")

    readline()
end

main()
