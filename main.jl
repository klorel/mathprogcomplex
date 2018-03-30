ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))


function main()

    ########################################
    # Construction du problème type
    rawproblem = buildPOP_1v1c()
    # rawproblem = buildPOP_1v2c()
    # rawproblem = buildPOP_2v3c()
    # rawproblem = buildPOP_WB2()

    rawproblem = buildPOPR_2v2c()

    ########################################
    # Normalizing pb and setting relaxation order by constraint
    problem = normalize_problem(rawproblem)
    relax_ctx = set_relaxation(problem, ismultiordered = false, hierarchykind=:Real, d = 1)

    println("\n")
    println(relax_ctx)

    println("\n")

    println(problem)

    ########################################
    # Construction du sparsity pattern, extension chordale, cliques maximales.

    max_cliques = get_allvars(relax_ctx, problem)
    println("\n--------------------------------------------------------")
    println("max cliques = $max_cliques")

    moments_params = build_sparsity(relax_ctx, problem, max_cliques)
    println("\n--------------------------------------------------------")
    println("moment params =")
    for (key, (val1, val2)) in moments_params
        println("$key \t -> $val2, $val1")
    end
    
    ########################################
    # Calcul des matrices de moment

    momentmatrices = compute_momentmat(relax_ctx, problem, moments_params, max_cliques)
    println("\n--------------------------------------------------------")
    println("moment matrices =")
    for (cstr, mm) in momentmatrices
        println("Constraint $cstr :")
        println(mm)
    end
    println("--------------------")
    


    # B_i = compute_Bibycstr(problem, momentmatrices, max_cliques, cliquevarsbycstr, orderbyclique, relax_ctx)
    
    # SDP_SOS = build_SDP_SOS(problem, max_cliques, B_i, cliquevarsbycstr, orderbyclique, relax_ctx)
    
    # ########################################
    # # Calcul d'une solution par un solveur
    # m, Zi, yα_re, yα_im, expo2int, int2expo = make_JuMPproblem(SDP_SOS, SCSSolver(max_iters=5000000, eps=1e-3, verbose=true), relax_ctx)

    # println("-----> SDP_SOS problem size: ", Base.summarysize(m)/1024, " ko")
    # println("-----> JuMP problem size: ", Base.summarysize(m)/1024, " ko")

    # ########################################
    # # Résolution du SDP par un solveur
    # println("-----> Starting solve")
    # solve(m)

    # println("\n-----> Objective value: ", getobjectivevalue(m), "\n")

    # # for (cstrname, mmb) in B_i
    # #     println("$cstrname \t= ", getvalue(Zi[cstrname]), "\n")
    # # end

    # println("\n\n----->Lagrange multipliers : yα =")
    # yα = - getdual(yα_re) - im*getdual(yα_im)
    # print_cmat(yα)

end

main()