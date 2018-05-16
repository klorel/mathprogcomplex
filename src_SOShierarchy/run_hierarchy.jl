
function run_hierarchy(problem::Problem, relax_ctx::RelaxationContext; logpath = "")
    ########################################
    # Construction du sparsity pattern, extension chordale, cliques maximales.
    max_cliques = get_maxcliques(relax_ctx, problem)
    # max_cliques = get_WB5cliques(relax_ctx, problem)

    ########################################
    # Compute moment matrices parameters: order et variables
    momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem, max_cliques)

    ########################################
    # Calcul des matrices de moment
    mmtrel_pb = MomentRelaxationPb(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)

    ########################################
    # Convert to a primal SDP problem
    sdpinstance = build_SDPInstance(relax_ctx, mmtrel_pb)

    export_SDP(relax_ctx, sdpinstance, logpath)
    sdp_instance = read_SDPInstance(logpath)

    ########################################
    # Build Mosek-like structure for SDP problem
    sdp = SDP_Problem()

    set_constraints!(sdp, sdp_instance)
    set_vartypes!(sdp, sdp_instance)
    set_blocks!(sdp, sdp_instance)

    set_matrices!(sdp, sdp_instance)
    set_linear!(sdp, sdp_instance)
    set_const!(sdp, sdp_instance)

    primal = SortedDict{Tuple{String,String,String}, Float64}()
    dual = SortedDict{Tuple{String, String, String}, Float64}()

    if logpath != ""
        primobj, dualobj = solve_mosek(sdp::SDP_Problem, primal, dual; logname = joinpath(logpath, "Mosek_run.log"))
    else
        primobj, dualobj = solve_mosek(sdp::SDP_Problem, primal, dual)
    end

    return primobj, dualobj
end