ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))


function run_hierarchy(problem::Problem, relax_ctx::RelaxationContext, logname)
    ########################################
    # Construction du sparsity pattern, extension chordale, cliques maximales.
    max_cliques = get_maxcliques(relax_ctx, problem)

    ########################################
    # Compute moment matrices parameters: order et variables
    moments_params = build_sparsity(relax_ctx, problem, max_cliques)

    ########################################
    # Calcul des matrices de moment
    mmtrel_pb = MomentRelaxationPb(relax_ctx, problem, moments_params, max_cliques)

    ########################################
    # Convert to a primal SDP problem
    sdpinstance = build_SDPInstance(relax_ctx, mmtrel_pb)

    export_SDP(relax_ctx, sdpinstance, pwd())
    sdp_instance = read_SDPInstance(pwd())

    ########################################
    # Build Mosek-like structure for SDP problem
    sdp = SDP_Problem()

    set_constraints!(sdp, sdp_instance)
    set_blocks!(sdp, sdp_instance)
    set_matrices!(sdp, sdp_instance)
    set_linear!(sdp, sdp_instance)
    set_const!(sdp, sdp_instance)

    primal = SortedDict{Tuple{String,String,String}, Float64}()
    dual = SortedDict{Tuple{String, String, String}, Float64}()

    primobj, dualobj = solve_mosek(sdp::SDP_Problem, primal, dual; logname = logname)

    return primobj, dualobj
end


function main()

    params = SortedSet([(-30.80, 2)])
                        # (-20.51, 2),
                        # (-10.22, 2),
                        # (0.07, 2),
                        # (10.36, 2),
                        # (20.65, 2),
                        # (30.94, 2),
                        # (41.23, 2),
                        # (51.52, 2),
                        # (61.81, 2)])

    primobjectives = SortedDict()

    for (q5min, dmax) in params
        problem = buildPOP_WB5(q5min=q5min)

        for d=2:2:2*dmax  # Complex d definition, twice real one...
            relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                                d = d)

            primobj, dualobj = run_hierarchy(problem, relax_ctx, "Mosek_WB5_v5max_$(q5min)_(d)_$(d).log")
            primobjectives[(q5min, d)] = primobj
        end
    end

    @printf("%5s  %4s  %s\n", "q5min", "dmax", "objective")
    for ((q5min, data), obj) in primobjectives
        @printf("%4f  %4i  %f\n", q5min, data, obj)
    end

    return primobjectives
end

main()
