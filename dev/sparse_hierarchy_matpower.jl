include(joinpath(pwd(),"src_PowSysMod", "PowSysMod_body.jl"))

function build_and_solve_matpower_relaxation(instance::String, sparse_param::Bool, order::Int64, QCQP_param::Bool)
    instance_path = joinpath(pwd(),"..","data", "data_Matpower", "matpower", instance)
    if QCQP_param
        problem_c, point = import_from_dat(joinpath("..", "data", "data_Matpower", "matpower_QCQP", instance[1:end-2]*".dat"))
        problem = pb_cplx2real(problem_c)
    else
        OPFpbs = load_OPFproblems(MatpowerInput, instance_path)
        ## Bulding optimization problem
        pb_global = build_globalpb!(OPFpbs)
        problem = pb_cplx2real(pb_global)
    end

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        # symmetries=[PhaseInvariance],
                                        issparse=sparse_param,
                                        d = order)

    println("\n--------------------------------------------------------")
    println("problem = \n$problem")

    println("\n--------------------------------------------------------")
    println("relax_ctx = \n$relax_ctx")

    ########################################
    # Construction du sparsity pattern, extension chordale, cliques maximales.
    if sparse_param
        # max_cliques = get_cliques(problem)
        max_cliques = get_cliques_matpower(instance_path)
    else
        max_cliques = get_maxcliques(relax_ctx, problem)
    end

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

end


instances = readir(joinpath("..","data", "data_Matpower", "blocks_AMD_clique"))
println(instances)
