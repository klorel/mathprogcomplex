ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))

using Convex, SCS

function build_Convex_SDPPrimal(primalSDP::MomentRelaxation{Float64}, mysolver)

    m = Model(solver = mysolver)

    # Collect all moments
    moments = Set{Moment}(keys(primalSDP.objective))
    for (key, mmtmat) in primalSDP.constraints
        for ((γ, δ), p) in mmtmat.mm
            union!(moments, keys(p))
        end
    end

    # Build LMI index mappings
    ctr_2_keymap = Dict()
    for (key, mmtmat) in primalSDP.constraints
        expokeys = Set{Exponent}()
        for (γ, δ) in keys(mmtmat.mm)
            union!(expokeys, [γ, δ])
        end

        expo_to_ind = Dict{Exponent, Int}()
        i=0
        for expo in sort(collect(expokeys))
            i+=1
            expo_to_ind[expo] = i
        end

        ctr_2_keymap[key] = expo_to_ind

        info("- ctr is: ", key)
        for expo in sort(collect(keys(expo_to_ind)))
            ind = expo_to_ind[expo]
            println("    $expo      $ind")
        end
    end

    #Define JuMP variables
    moment_to_JuMPvar = Dict{Moment, JuMP.Variable}()
    for moment in sort(collect(moments))
        moment_to_JuMPvar[moment] = @variable(m, basename=format_string(moment.conj_part)*format_string(moment.expl_part)*moment.clique)
    end

    # Define objective
    objective = 0
    for (mmt, val) in primalSDP.objective
        objective += moment_to_JuMPvar[mmt]*val
    end
    @objective(m, Min, objective)


    # Define constraints
    warn("\nBuilding constraints")
    for ((ctr, clique), mmtmat) in primalSDP.constraints
        println("   -> $ctr, $clique")

        moment_to_scalmat = Dict{Moment, Any}()

        for ((γ, δ), p) in mmtmat.mm
            for (moment, val) in p
                !haskey(moment_to_scalmat, moment) && moment_to_scalmat[moment] = ([],[],[])
                moment_to_scalmat[moment][1] = ctr_2_keymap[(ctr, clique)][γ]
                moment_to_scalmat[moment][2] = ctr_2_keymap[(ctr, clique)][δ]
                moment_to_scalmat[moment][3] = val
            end
        end

        LMIbody = 0
        for (moment, triplets) in moment_to_scalmat
            LMIbody += moment_to_JuMPvar[moment] * sparse(triplets...)
        end
        @constraint(m, SDP, LMIbody)
    end

    return m
end

sparse([1, 2, 3], [1, 2, 3], [0, 2, 0])


function main()
    problem = buildPOP_WB2(v2max=1.022, setnetworkphase=true)
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

    println("\n\n\n")
    m = build_JuMP_SDPDual(mmtrel_pb, SCSSolver())
    show(m)
end

main()