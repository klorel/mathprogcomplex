
function run_hierarchy(problem::Problem, relax_ctx::RelaxationContext, logpath; indentedprint=false,
                                                                                max_cliques::SortedDict{String, SortedSet{Variable}}=SortedDict{String, SortedSet{Variable}}())

    ########################################
    # Construction du sparsity pattern, extension chordale, cliques maximales.
    if max_cliques == Dict{String, SortedSet{Variable}}()
        max_cliques = get_maxcliques(relax_ctx, problem)
    end

    ########################################
    # Compute moment matrices parameters: order et variables
    momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem, max_cliques)

    ########################################
    # Calcul des matrices de moment
    mmtrel_pb = MomentRelaxationPb(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)

    ########################################
    # Convert to a primal SDP problem
    sdpinstance = build_SDPInstance(relax_ctx, mmtrel_pb)

    export_SDP(relax_ctx, sdpinstance, logpath, indentedprint=indentedprint)
    sdp_instance = read_SDPInstance(logpath)

    ########################################
    # Build Mosek-like structure for SDP problem
    sdp = SDP_Problem()

    set_constraints!(sdp, sdp_instance)
    set_vartypes!(sdp, sdp_instance)
    set_blocks!(sdp, sdp_instance)
    set_linvars!(sdp, sdp_instance)

    set_matrices!(sdp, sdp_instance)
    set_linear!(sdp, sdp_instance)
    set_const!(sdp, sdp_instance)

    primal = SortedDict{Tuple{String,String,String}, Float64}()
    dual = SortedDict{Tuple{String, String, String}, Float64}()

    primobj, dualobj = solve_mosek(sdp::SDP_Problem, primal, dual; logname = joinpath(logpath, "Mosek_run.log"))

    info("Writing results to $logpath")
    params_file = joinpath(logpath, "maxcliques_relaxctx.txt")
    isfile(params_file) && rm(params_file)
    open(params_file, "w") do fcliques
        println(fcliques, "# max_cliques are:")
        println(fcliques, max_cliques)
        println(fcliques, "# relaxation_ctx is:")
        println(fcliques, relax_ctx)
    end

    return primobj, dualobj
end