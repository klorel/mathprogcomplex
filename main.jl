ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))


function main()

    ########################################
    # Construction du problème type
    # problem = buildPOP_1v1c()
    # problem = buildPOPR_2v1c()
    # problem = buildPOP_1v2c()
    # problem = buildPOP_2v3c()
    # problem = buildPOP_WB2()
    # problem = buildPOP_WB2_expl()

    ########################################
    # Normalizing pb and setting relaxation order by constraint
    # relax_ctx = set_relaxation(problem, hierarchykind=:Complex, d = 1)
    
    real_pb = true
    if real_pb
        # Build the init problem and set relaxation parameters
        problem = buildPOPR_2v2cbis()
        relax_ctx = set_relaxation(problem; hierarchykind=:Real, 
                                            d = 1,
                                            symmetries = [PhaseInvariance])
    else
        # Build the init problem and set relaxation parameters
        problem = buildPOP_WB2_expl()
        relax_ctx = set_relaxation(problem; hierarchykind=:Complex,
                                            d = 2,
                                            symmetries = [PhaseInvariance])
        relax_ctx.di[get_momentcstrname()] = 2
    end

    println("\n--------------------------------------------------------")
    println("problem = \n$problem")
    
    println("\n--------------------------------------------------------")
    println("relax_ctx = \n$relax_ctx")
    
    ########################################
    # Build sparsity pattern, compute maximal cliques
    max_cliques = get_maxcliques(relax_ctx, problem)
    
    println("\n--------------------------------------------------------")
    println("max cliques =")
    for (cliquename, vars) in max_cliques
        print("$cliquename = ")
        for var in vars print("$var, ") end
        @printf("\b\b \n")
    end
    
    ########################################
    # Compute moment matrices parameters: order et variables
    moments_params = build_sparsity(relax_ctx, problem, max_cliques)
    println("\n--------------------------------------------------------")
    println("Matrix moment parameters =")
    for (key, (val1, val2)) in moments_params
        print("$key \t -> di-ki = $val2, \tcliques = ")
        for clique in val1 print("$clique, ") end
        @printf("\b\b \n")
    end
    
    ########################################
    # Compute partial moment hierarchy
    mmtrel_pb = MomentRelaxationPb(relax_ctx, problem, moments_params, max_cliques)
    println("\n--------------------------------------------------------")
    println("mmtrel_pb = $mmtrel_pb")
    
    ########################################
    # Convert to a primal SDP problem
    sdp = build_SDP(relax_ctx, mmtrel_pb)
    println("\n--------------------------------------------------------")
    println("sdp = \n$sdp")

    export_SDP(relax_ctx, sdp, pwd())

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