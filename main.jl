ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))
include(joinpath(ROOT, "dev", "get_cliques.jl"))


function main()

    # problem = buildPOP_WB2(v2max=1.022, setnetworkphase=false)
    # relax_ctx = set_relaxation(problem; hierarchykind=:Real,
    #                                     # symmetries=[PhaseInvariance],
    #                                     d = 1)
    #
    # problem = buildPOP_WB2(v2max=1.022, setnetworkphase=true)
    # relax_ctx = set_relaxation(problem; hierarchykind=:Real,
    #                                     # symmetries=[PhaseInvariance],
    #                                     d = 1)

    ###GOC
    data_path = joinpath("..", "data", "data_GOC")
    folder = "Phase_0_IEEE14_1Scenario"
    scenario = "scenario_1"
    folder_path = joinpath(data_path, folder)
    instance_path = joinpath(folder_path, scenario)
    raw = "powersystem.raw"
    gen = "generator.csv"
    con = "contingency.csv"
    rawfile = joinpath(instance_path,raw)
    genfile = joinpath(instance_path, gen)
    contfile = joinpath(instance_path, con)
    OPFpbs = load_OPFproblems(rawfile, genfile, contfile)
    introduce_Sgenvariables!(OPFpbs)
    ## Bulding optimization problem
    pb_global = build_globalpb!(OPFpbs)
    pb_global_real = pb_cplx2real(pb_global)
    problem = pb_global_real
    # problem = convert_mipb_to_pb(pb_global_real)
    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        # symmetries=[PhaseInvariance],
                                        issparse=true,
                                        d = 2)



    println("\n--------------------------------------------------------")
    println("problem = \n$problem")

    println("\n--------------------------------------------------------")
    println("relax_ctx = \n$relax_ctx")

    ########################################
    # Construction du sparsity pattern, extension chordale, cliques maximales.
    # max_cliques = get_maxcliques(relax_ctx, problem)
    max_cliques = get_cliques(problem)

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
    mmtrel_pb = MomentRelaxation(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)

    println("\n--------------------------------------------------------")
    println("mmtrel_pb = $mmtrel_pb")

    ########################################
    # Convert to a primal SDP problem
    sdpinstance = build_SDPInstance(relax_ctx, mmtrel_pb)

    # println("\n--------------------------------------------------------")
    # println("sdpinstance = \n$sdpinstance")

    path = joinpath(pwd(), "Mosek_runs", "worksdp")
    mkpath(path)
    export_SDP(sdpinstance, path)
    sdp_instance = read_SDPInstance(path)

    println("VAR_TYPES size:     $(size(sdp_instance.VAR_TYPES))")
    println("BLOCKS size:        $(size(sdp_instance.BLOCKS))")
    println("LINEAR size:        $(size(sdp_instance.LINEAR))")
    println("CONST size:         $(size(sdp_instance.CONST))")

    sdp = SDP_Problem()

    set_constraints!(sdp, sdp_instance)
    set_vartypes!(sdp, sdp_instance)
    set_blocks!(sdp, sdp_instance)
    set_linvars!(sdp, sdp_instance)

    set_matrices!(sdp, sdp_instance)
    set_linear!(sdp, sdp_instance)
    set_const!(sdp, sdp_instance)

    # println(sdp)

    primal = SortedDict{Tuple{String,String,String}, Float64}()
    dual = SortedDict{Tuple{String, String, String}, Float64}()

    primobj, dualobj = solve_mosek(sdp::SDP_Problem, primal, dual)

    # # println("Primal solution")
    # # for ((blockname, var1, var2), val) in primal
    # # @printf("%15s %5s %5s %f\n", blockname, var1, var2, val)
    # # end

    # # println("\nDual solution NEGATED")
    # # for var in problem.variables
    # #     ctrname = get_momentcstrname()
    # #     var1 = var[1]
    # #     var2 = "1"
    # #     val = dual[(ctrname, var1, var2)]
    # #     println("($(ctrname), $(var1), $(var2)) = $(-val)")
    # # end

    # println("Objectives : $primobj, $dualobj")
end

main()
