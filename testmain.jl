ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))
include(joinpath("dev","get_cliques.jl"))

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
    input = MatpowerInput
    if input == GOCInput
        folder = "Phase_0_IEEE14_1Scenario"
        # folder = "Phase_0_IEEE14"
        # folder = "Phase_0_RTS96"
        scenario = "scenario_1"
        instance_path = joinpath(pwd(),"..","data", "data_GOC", folder, scenario)
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
        problem = pb_cplx2real(pb_global)
    elseif input == MatpowerInput
        typeofinput = MatpowerInput
        data_path = joinpath(ROOT,"..","data","data_Matpower","matpower")
        instance = "WB2.m"
        instance_path = joinpath(data_path, instance)
        OPFpbs = load_OPFproblems(typeofinput, instance_path)
        ## Bulding optimization problem
        pb_global = build_globalpb!(OPFpbs)
        problem = pb_cplx2real(pb_global)
    end

    change_eq_to_ineq!(problem)

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            d = 1)
                                            # symmetries = [PhaseInvariance])


    println("\n--------------------------------------------------------")
    println("problem = \n$problem")

    println("\n--------------------------------------------------------")
    println("relax_ctx = \n$relax_ctx")

    ########################################
    # Construction du sparsity pattern, extension chordale, cliques maximales.
    max_cliques = get_maxcliques(relax_ctx, problem)
    # max_cliques = get_cliques(problem)

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
    println("mmtrel_pb = $mmtrel_pb")

    ########################################
    # Convert to a primal SDP problem
    sdpinstance = build_SDPInstance(relax_ctx, mmtrel_pb)
    println("\n--------------------------------------------------------")
    println("sdpinstance = \n$sdpinstance")
    export_SDP(relax_ctx, sdpinstance, pwd())

    sdp_instance = read_SDPInstance(pwd())

    println("VAR_TYPES size:     $(size(sdp_instance.VAR_TYPES))")
    println("BLOCKS size:        $(size(sdp_instance.BLOCKS))")
    println("LINEAR size:        $(size(sdp_instance.LINEAR))")
    println("CONST size:         $(size(sdp_instance.CONST))")

    sdp = SDP_Problem()

    set_constraints!(sdp, sdp_instance, debug=true)
    set_blocks!(sdp, sdp_instance, debug=true)
    set_matrices!(sdp, sdp_instance, debug=true)
    set_linear!(sdp, sdp_instance, debug=true)
    set_const!(sdp, sdp_instance, debug=true)


    primal=SortedDict{Tuple{String,String,String}, Float64}()
    dual=SortedDict{String, Float64}()

    solve_mosek(sdp::SDP_Problem, primal::SortedDict{Tuple{String,String,String}, Float64}, dual::SortedDict{String, Float64}, debug=false)
end

main()
