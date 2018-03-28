using JuMP, SCS

ROOT = pwd()
include(joinpath("set_problems.jl"))
include(joinpath("func_definitions.jl"))

function main()

    ########################################
    # Construction du problème type
    # rawproblem = buildPOP_1v1c()
    rawproblem = buildPOP_1v2c()
    # rawproblem = buildPOP_2v3c()
    # rawproblem = buildPOP_WB2()

    ########################################
    # Normalizing pb and setting relaxation order by constraint
    problem = normalize_problem(rawproblem)
    relax_ctx = set_relaxation(problem, issparse = false, ismultiordered = false, d = 2)

    println(problem)

    ########################################
    # Construction du sparsity pattern
    sparsity_pattern = compute_sparsitypattern(problem, relax_ctx)

    # Extension chordale et détection des cliques maximales
    compute_chordalextension!(sparsity_pattern)
    max_cliques = compute_maxcliques(sparsity_pattern)

    ########################################
    # Relaxation degree par clique and variables par constrainte
    varsbycstr = compute_varsbycstr(problem)
    cliquevarsbycstr = compute_varsbycstr(sparsity_pattern, max_cliques, varsbycstr)
    orderbyclique = compute_cliqueorders(sparsity_pattern, varsbycstr, max_cliques, relax_ctx)

    ########################################
    # Calcul des matrices B_i et pose du probleme
    momentmatrices = compute_momentmat(problem, max_cliques, cliquevarsbycstr, orderbyclique, relax_ctx)

    println("-------------------- momentmatrices:")
    for (cstr, mm) in momentmatrices
        println("$cstr :")
        println(mm)
    end
    println("--------------------")
    
    B_i = compute_Bibycstr(problem, momentmatrices, max_cliques, cliquevarsbycstr, orderbyclique, relax_ctx)
    
    println("-------------------- B_i:")
        for (cstr, mmb) in B_i
        println("$cstr :")
        println(mmb)
    end
    println("--------------------")
    
    
    SDP_SOS = build_SDP_SOS(problem, max_cliques, B_i, cliquevarsbycstr, orderbyclique, relax_ctx)
    
    println("-------------------- SDP_SOS:")
    println(SDP_SOS)
    println("--------------------")
    
    ########################################
    # Calcul d'une solution par un solveur
    m, Zi, yα_re, yα_im, expo2int, int2expo = make_JuMPproblem(SDP_SOS, SCSSolver(max_iters=5000000, eps=1e-3, verbose=true), relax_ctx)

    println("-----> SDP_SOS problem size: ", Base.summarysize(m)/1024, " ko")
    println("-----> JuMP problem size: ", Base.summarysize(m)/1024, " ko")

    ########################################
    # Résolution du SDP par un solveur
    println("\n-----> Starting solve")
    solve(m)

    println("\n\n-----> Objective value: ", getobjectivevalue(m), "\n")

    # for (cstrname, mmb) in B_i
    #     println("$cstrname \t= ", getvalue(Zi[cstrname]), "\n")
    # end

    println("\n\n----->Lagrange multipliers : yα =")
    yα = - getdual(yα_re) - im*getdual(yα_im)
    print_cmat(yα)

    # ########################################
    # # Reconstruction de la solution du POP
    # zsol = recover_zsol(yα, expo2int, problem)
    # println("\n\n----->Z solution : z =$zsol")

    # println("\nMoment matrix:")
    # println(evaluate(momentmatrices["moment_cstr"]), zsol)

    # for cstrname in setdiff(keys(momentmatrices), Set(["moment_cstr"]))
    #     println("Localizing matrix : $cstrname ")
    #     println(evaluate(momentmatrices[cstrname], zsol))
    # end

end

main()