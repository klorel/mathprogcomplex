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

    # problem, relax_ctx = lasserre_ex1()
    # problem, relax_ctx = lasserre_ex2()
    # problem, relax_ctx = lasserre_ex3()
    # problem, relax_ctx = lasserre_ex5(d=4)

    # problem = buildPOP_1v2()

    # problem = buildPOP_WB2(v2max=0.983, rmineqs=true)
    problem = buildPOP_WB5()
    # problem = buildPOP_case9()

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        d = 2,
                                        # symmetries=[PhaseInvariance],
                                        issparse=true)


    println("\n--------------------------------------------------------")
    println("problem = \n$problem")

    println("\n--------------------------------------------------------")
    println("relax_ctx = \n$relax_ctx")

    ########################################
    # Construction du sparsity pattern, extension chordale, cliques maximales.
    # max_cliques = get_maxcliques(relax_ctx, problem)
    max_cliques = get_WB5cliques(relax_ctx, problem)
    # max_cliques = get_case9cliques(relax_ctx, problem)

    println("\n--------------------------------------------------------")
    println("max cliques =")
    for (cliquename, vars) in max_cliques
        print("$cliquename = ")
        for var in vars print("$var, ") end
        @printf("\b\b \n")
    end

    ########################################
    # Compute moment and localizing matrices parameters: order et variables
    momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem, max_cliques)
    println("\n--------------------------------------------------------")
    println("moment params =")
    for (cliquename, dcl) in momentmat_param
        println("Moment matrix, $cliquename \t -> dcl = $dcl")
    end
    for (key, (val1, val2)) in localizingmat_param
        @printf("%15s \t -> di-ki = %i, \tcliques = ", key, val2)
        for clique in val1 print("$clique, ") end
        @printf("\b\b \n")
    end

    ########################################
    # Build the moment relaxation problem
    mmtrel_pb = MomentRelaxationPb(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques; verbose=true)
    println("\n--------------------------------------------------------")
    println("mmtrel_pb = $mmtrel_pb")

    # ########################################
    # # Convert to a primal SDP problem
    sdpinstance = build_SDPInstance(relax_ctx, mmtrel_pb)
    # println("\n--------------------------------------------------------")
    # println("sdpinstance = \n$sdpinstance")

    export_SDP(relax_ctx, sdpinstance, pwd())
    sdp_instance = read_SDPInstance(pwd())

    println("VAR_TYPES size:     $(size(sdp_instance.VAR_TYPES))")
    println("BLOCKS size:        $(size(sdp_instance.BLOCKS))")
    println("LINEAR size:        $(size(sdp_instance.LINEAR))")
    println("CONST size:         $(size(sdp_instance.CONST))")

    sdp = SDP_Problem()

    set_constraints!(sdp, sdp_instance)
    set_vartypes!(sdp, sdp_instance)
    set_blocks!(sdp, sdp_instance)

    # # warn("sdpvars : ")
    # # for (blockname, block) in sdp.name_to_sdpblock
    # #     println("$blockname \t $block")
    # # end

    # # warn("symvars : ")
    # # for (blockname, block) in sdp.name_to_symblock
    # #     println("$blockname \t $block")
    # # end
    # # @show sdp.n_scalvarsym


    set_matrices!(sdp, sdp_instance)

    # warn("sdp.matrices")
    # for ((ctrname, blockname, var1, var2), coeff) in sdp.matrices
    #     @printf("%45s %30s %20s %20s -> %f\n", ctrname, blockname, var1, var2, coeff)
    # end

    warn("sdp.lin_matsym")
    for ((ctrname, blockname, var1, var2), coeff) in sdp.lin_matsym
        @printf("%45s %30s %20s %20s -> %f\n", ctrname, blockname, var1, var2, coeff)
    end

    set_linear!(sdp, sdp_instance)

    warn("sdp.lin")
    for ((ctrname, var), coeff) in sdp.linear
        @printf("%45s %30s -> %f\n", ctrname, var, coeff)
    end

    warn("sdp.scalvar_to_id")
    for (varname, id) in sdp.scalvar_to_id
        @printf("%45s -> %i\n", varname, id)
    end

    set_const!(sdp, sdp_instance)

    # println(sdp)

    primal = SortedDict{Tuple{String,String,String}, Float64}()
    dual = SortedDict{Tuple{String, String, String}, Float64}()

    primobj, dualobj = solve_mosek(sdp::SDP_Problem, primal, dual, debug=false)

    # println("Primal solution")
    # for ((blockname, var1, var2), val) in primal
    # @printf("%15s %5s %5s %f\n", blockname, var1, var2, val)
    # end

    # println("\nDual solution NEGATED")
    # for var in problem.variables
    #     ctrname = get_momentcstrname()
    #     var1 = var[1]
    #     var2 = "1"
    #     val = dual[(ctrname, var1, var2)]
    #     println("($(ctrname), $(var1), $(var2)) = $(-val)")
    # end

    # println("Objectives : $primobj, $dualobj")
end

main()
