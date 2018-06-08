ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))

using Convex, Mosek

function set_variables!(matrix_variables, scalar_variables, primalSDP)
    for (names, block) in primalSDP.name_to_sdpblock
        n = length(block.var_to_id)
        matrix_variables[names] = Convex.Variable(n, n)
    end

    for varname in keys(primalSDP.scalvar_to_id)
        scalar_variables[varname] = Convex.Variable()
    end
end

function build_Convex_SDPPrimal(primalSDP::SDP_Problem{Float64}, mysolver)

    matrix_variables = Dict{String, Any}()
    scalar_variables = Dict{String, Any}()

    set_variables!(matrix_variables, scalar_variables, primalSDP)

    objective = 0
    constraints = Dict{SDP_Moment, Any}()
    for ((sdpmoment, matname, γ, δ), λ) in primalSDP.matrices
        i = primalSDP.name_to_sdpblock[matname].var_to_id[γ]
        j = primalSDP.name_to_sdpblock[matname].var_to_id[δ]
        n = size(matrix_variables[matname], 1)

        if sdpmoment in primalSDP.obj_keys

            objective += vecdot(matrix_variables[matname], sparse([i], [j], [λ], n, n))

            # objective += matrix_variables[matname][primalSDP.name_to_sdpblock[matname].var_to_id[γ],
            #                                        primalSDP.name_to_sdpblock[matname].var_to_id[δ]] * λ
            println(" Adding to objective, Zi $matname, i,j = $i, $j  *\t $λ")
        else
            !haskey(constraints, sdpmoment) && (constraints[sdpmoment] = 0)
            # constraints[sdpmoment] += matrix_variables[matname] * sparse([primalSDP.name_to_sdpblock[matname].var_to_id[γ]],
            #                                                             [primalSDP.name_to_sdpblock[matname].var_to_id[δ]],
            #                                                             [λ],
            #                                                             size(matrix_variables[matname])...)
            # constraints[sdpmoment] += matrix_variables[matname][primalSDP.name_to_sdpblock[matname].var_to_id[γ],
            #                                                     primalSDP.name_to_sdpblock[matname].var_to_id[δ]] * λ

            constraints[sdpmoment] += vecdot(matrix_variables[matname], sparse([i], [j], [λ], n, n))
            println(" Adding to constraint $sdpmoment, Zi $matname, i,j = $i, $j  *\t $λ")

        end
    end

    for ((sdpmoment, scalvar), λ) in primalSDP.linear
        if sdpmoment in primalSDP.obj_keys
            objective += scalar_variables[scalvar] * λ
            println("Ading to objective, $scalvar \t $λ")
        else
            !haskey(constraints, sdpmoment) && (constraints[sdpmoment] = 0)
            constraints[sdpmoment] += scalar_variables[scalvar] * λ
        end
    end

    for (sdpmoment, λ) in primalSDP.cst_ctr
        if sdpmoment in primalSDP.obj_keys && λ!=0
            objective += λ
            println("Adding to objective $λ")
        elseif λ != 0.0
            !haskey(constraints, sdpmoment) && (constraints[sdpmoment] = 0)
            constraints[sdpmoment] += λ
        end
    end

    ctrs = Convex.Constraint[]

    for ctrname in sort(collect(keys(constraints)))
        push!(ctrs, (constraints[ctrname] == 0))
    end
    for (matname, var) in matrix_variables
        push!(ctrs, var in :SDP)
    end

    problem = minimize(objective, ctrs)
    return problem
end



function main()
    problem = buildPOP_WB2(v2max=1.022, setnetworkphase=false, rmeqs=true)
    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        # symmetries=[PhaseInvariance],
                                        d = 1)

    println("\n--------------------------------------------------------")
    println("problem = \n$problem")

    max_cliques = get_maxcliques(relax_ctx, problem)

    println("\n--------------------------------------------------------")
    println("max cliques =\n$max_cliques")

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

    # Convert to a primal SDP problem
    sdpinstance = build_SDPInstance(relax_ctx, mmtrel_pb)

    logpath=joinpath("Mosek_runs", "SCSdev")
    ispath(logpath) && rm(logpath, recursive=true)
    mkpath(logpath)
    export_SDP(sdpinstance, logpath, renamemoments=false)

    main2()
end

function main2()
    logpath=joinpath("Mosek_runs", "SCSdev")
    sdp = build_mosekpb(logpath)

    println("\n\n\n")
    m = build_Convex_SDPPrimal(sdp, MosekSolver())

    solve!(m)
    return m.status
end


# main()
# main2()

function main3()
    problem = buildPOP_WB2(rmeqs=false) # v2max=0.976

    logpath = joinpath("Mosek_runs", "totoa")
    mkpath(logpath)
    export_to_dat(problem, joinpath("Mosek_runs", "totoa"))
end