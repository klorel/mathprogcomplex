ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))

using Base.Profile
using ProfileView

function toprofile(n, problem, d)
    # problem = buildPOP_WB2(v2max=1.022, setnetworkphase=false)
    # problem = buildPOP_WB5()

    for i=1:n
        sparse_param = true
        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            # symmetries=[PhaseInvariance],
                                            issparse=sparse_param,
                                            d = 1)



        println("\n--------------------------------------------------------")
        println("problem = \n$problem")

        println("\n--------------------------------------------------------")
        println("relax_ctx = \n$relax_ctx")

        ########################################
        # Construction du sparsity pattern, extension chordale, cliques maximales.
        if sparse_param
            # max_cliques = get_cliques(problem)
            # max_cliques = get_cliques_matpower(instance_path)
            # max_cliques = get_cliques_matpower_forQCQP(instance_path)
            max_cliques = SortedDict{String, SortedSet{Variable}}()
            max_cliques["clique1"] = SortedSet{Variable}()
            for var in problem.variables
                push!(max_cliques["clique1"], Variable(var[1], var[2]))
            end
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
end


function main()

    OPFpbs = load_OPFproblems(MatpowerInput, joinpath("..", "data", "data_Matpower", "matpower", "case14.m"))
    problem_c = build_globalpb!(OPFpbs)
    problem = pb_cplx2real(problem_c)

    toprofile(1, problem, 1)

    Profile.clear()
    @profile toprofile(3, problem, 1)
    ProfileView.view()
end


main()
