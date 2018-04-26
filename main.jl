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

    # pb_ind = 2
    # if pb_ind == 0
    #     # Build the init problem and set relaxation parameters
    #     # problem = buildPOPR_2v2cbis()
    #     problem = buildPOPR_2v2c()
    #     change_eq_to_ineq!(problem)

    #     relax_ctx = set_relaxation(problem; hierarchykind=:Real,
    #                                         d = 1)
    #                                         # symmetries = [PhaseInvariance])
    # elseif pb_ind == 1
    #     # Build the init problem and set relaxation parameters
    #     problem = buildPOP_WB2_expl()
    #     relax_ctx = set_relaxation(problem; hierarchykind=:Complex,
    #                                         d = 2,
    #                                         symmetries = [PhaseInvariance])
    #     relax_ctx.di[get_momentcstrname()] = 2
    # else
    #     # WB2 problem converted to real
    #     WB2_C = buildPOP_WB2_expl()
    #     change_eq_to_ineq!(WB2_C)

    #     problem = pb_cplx2real(WB2_C)

    #     relax_ctx = set_relaxation(problem; hierarchykind=:Real,
    #                                         d = 2)
    #                                         # symmetries = [PhaseInvariance])
    # end

    # problem, relax_ctx = lasserre_ex1()
    # problem, relax_ctx = lasserre_ex2()
    # problem, relax_ctx = lasserre_ex3()
    # problem, relax_ctx = lasserre_ex5(d=3)

    # problem_C = buildPOP_EllJoszMolc()
    # change_eq_to_ineq!(problem_C)
    # problem = pb_cplx2real(problem_C)

    # relax_ctx = set_relaxation(problem; hierarchykind=:Real,
    #                                     d = 4)

    # problem_C = buildPOP_WB2_expl()
    # change_eq_to_ineq!(problem_C)
    # problem = pb_cplx2real(problem_C)

    # relax_ctx = set_relaxation(problem; hierarchykind=:Real,
    #                                     d = 2,
    #                                     symmetries = [PhaseInvariance])

    # problem = buildPOP_WB2R_expl()
    # change_eq_to_ineq!(problem)
    # relax_ctx = set_relaxation(problem; hierarchykind=:Real,
    #                                     d = 3)
    #                                     # symmetries = [PhaseInvariance])


    println("\n--------------------------------------------------------")
    println("problem = \n$problem")

    println("\n--------------------------------------------------------")
    println("relax_ctx = \n$relax_ctx")

    ########################################
    # Construction du sparsity pattern, extension chordale, cliques maximales.
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
    println("moment params =")
    for (key, (val1, val2)) in moments_params
        print("$key \t -> di-ki = $val2, \tcliques = ")
        for clique in val1 print("$clique, ") end
        @printf("\b\b \n")
    end

    ########################################
    # Calcul des matrices de moment

    mmtrel_pb = MomentRelaxationPb(relax_ctx, problem, moments_params, max_cliques)
    println("\n--------------------------------------------------------")
    # println("mmtrel_pb = $mmtrel_pb")

    ########################################
    # Convert to a primal SDP problem
    sdpinstance = build_SDPInstance(relax_ctx, mmtrel_pb)
    println("\n--------------------------------------------------------")
    # println("sdpinstance = \n$sdpinstance")
    export_SDP(relax_ctx, sdpinstance, pwd())

    sdp_instance = read_SDPInstance(pwd())

    println("VAR_TYPES size:     $(size(sdp_instance.VAR_TYPES))")
    println("BLOCKS size:        $(size(sdp_instance.BLOCKS))")
    println("LINEAR size:        $(size(sdp_instance.LINEAR))")
    println("CONST size:         $(size(sdp_instance.CONST))")

    sdp = SDP_Problem()

    set_constraints!(sdp, sdp_instance)
    set_blocks!(sdp, sdp_instance)
    set_matrices!(sdp, sdp_instance)
    set_linear!(sdp, sdp_instance)
    set_const!(sdp, sdp_instance)

    # println("SDP_Problem :\n$sdp")

    primal=SortedDict{Tuple{String,String,String}, Float64}()
    dual=SortedDict{String, Float64}()

    solve_mosek(sdp::SDP_Problem, primal, dual, debug=false, print_sol=false)
end

main()
