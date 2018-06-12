ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))


function main()

    problem = buildPOP_WB5(q5min=-20.51, rmeqs=false) #v2max=0.983

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        # symmetries=[PhaseInvariance],
                                        issparse=true,
                                        d = 1)

    println("\n--------------------------------------------------------")
    println("problem = \n$problem")

    println("\n--------------------------------------------------------")
    println("relax_ctx = \n$relax_ctx")

    ########################################
    # Construction du sparsity pattern, extension chordale, cliques maximales.
    # max_cliques = get_maxcliques(relax_ctx, problem)
    max_cliques = get_WB5cliques(relax_ctx, problem)

    println("\n--------------------------------------------------------")
    println("max cliques =")
    print(max_cliques)

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
    mmtrel_pb = MomentRelaxation{Float64}(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)

    println("\n--------------------------------------------------------")
    println("mmtrel_pb = $mmtrel_pb")

    ########################################
    # Convert to a primal SDP problem
    sdpinstance = build_SOSrelaxation(relax_ctx, mmtrel_pb)

    println("\n--------------------------------------------------------")
    println("sdpinstance = \n$sdpinstance")

    path = joinpath(pwd(), "Mosek_runs", "worksdp")
    ispath(path) && rm(path, recursive=true)
    mkpath(path)
    export_SDP(sdpinstance, path, renamemoments=false)

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

end

main()
