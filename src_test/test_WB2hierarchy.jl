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

    params = SortedSet([(0.976, 2),
                        (0.983, 2),
                        (0.989, 2),
                        (0.996, 2),
                        (1.002, 2),
                        (1.009, 2),
                        (1.015, 2),
                        (1.022, 3),
                        (1.028, 3),
                        (1.035, 2)])

    primobjectives = SortedDict()

    for (v2max, dmax) in params
        problem = buildPOP_WB2(v2max=v2max)

        for d=2:2:2*dmax  # Complex d definition, twice real one...
            relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                                d = d)

            primobj, dualobj = run_hierarchy(problem, relax_ctx, "Mosek_WB2_v2max_$(v2max)_(d)_$(d).log")
            primobjectives[(v2max, d)] = primobj
        end
    end

    @printf("%5s  %4s  %s\n", "v2max", "dmax", "objective")
    for ((v2max, data), obj) in primobjectives
        @printf("%4f  %4i  %f\n", v2max, data, obj)
    end

    return primobjectives
end

main()
