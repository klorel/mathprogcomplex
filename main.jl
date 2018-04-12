ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))


function main()

    ########################################
    # Construction du problème type
    # rawproblem = buildPOP_1v1c()
    # rawproblem = buildPOPR_2v1c()
    # rawproblem = buildPOP_1v2c()
    # rawproblem = buildPOP_2v3c()
    # rawproblem = buildPOP_WB2()
    # rawproblem = buildPOP_WB2_expl()

    ########################################
    # Normalizing pb and setting relaxation order by constraint
    # problem = normalize_problem(rawproblem)
    # relax_ctx = set_relaxation(problem, hierarchykind=:Complex, d = 1)
    
    real_pb = false
    if real_pb
        # Construction of the initial problem
        rawproblem = buildPOPR_2v2cbis()

        # Transform the problem to canonical form and set relaxation parameters
        problem = normalize_problem(rawproblem)
        relax_ctx = set_relaxation(problem, hierarchykind=:Real, d = 2)
    else
        # Construction of the initial problem
        rawproblem = buildPOP_WB2_expl()

        # Transform the problem to canonical form and set relaxation parameters
        problem = normalize_problem(rawproblem)
        relax_ctx = set_relaxation(problem; hierarchykind=:Complex,
                                            d = 2,
                                            symmetries = [PhaseInvariance])
        relax_ctx.di["moment_cstr"] = 2
    end

    println("\n--------------------------------------------------------")
    println("relax_ctx = \n$relax_ctx")
    
    println("\n--------------------------------------------------------")
    println("problem = \n$problem")
    
    ########################################
    # Build sparsity pattern, compute maximal cliques
    max_cliques = get_maxcliques(relax_ctx, problem)
    
    println("\n--------------------------------------------------------")
    println("max cliques = \n$max_cliques")
    
    ########################################
    # Compute moment matrices parameters: order et variables
    moments_params = build_sparsity(relax_ctx, problem, max_cliques)
    println("\n--------------------------------------------------------")
    println("Matrix moment parameters =")
    for (key, (val1, val2)) in moments_params
        println("$key \t -> $val1, di-ki=$val2")
    end
    
    ########################################
    # Compute partial moment hierarchy
    mmtrel_pb = MomentRelaxationPb(relax_ctx, problem, moments_params, max_cliques)
    println("\n--------------------------------------------------------")
    println("mmtrel_pb = $mmtrel_pb")
    
    ########################################
    # Convert to a primal SDP problem
    SDP_body, SDP_rhs = build_SDP(relax_ctx, mmtrel_pb)
    println("\n--------------------------------------------------------")
    println("SDP_body = \n$SDP_body")

    println("\n--------------------------------------------------------")
    println("SDP_rhs = \n$SDP_rhs")


    ########################################
    # Calcul d'une solution par un solveur
    # m, Zi, yα_re, yα_im, expo2int, int2expo = make_JuMPproblem(SDP_SOS, SCSSolver(max_iters=5000000, eps=1e-3, verbose=true), relax_ctx)

    # println("-----> SDP_SOS problem size: ", Base.summarysize(m)/1024, " ko")
    # println("-----> JuMP problem size: ", Base.summarysize(m)/1024, " ko")

    ########################################
    # Résolution du SDP par un solveur
    # println("-----> Starting solve")
    # solve(m)

    # println("\n-----> Objective value: ", getobjectivevalue(m), "\n")

    # # for (cstrname, mmb) in B_i
    # #     println("$cstrname \t= ", getvalue(Zi[cstrname]), "\n")
    # # end

    # println("\n\n----->Lagrange multipliers : yα =")
    # yα = - getdual(yα_re) - im*getdual(yα_im)
    # print_cmat(yα)

    return
end

main()