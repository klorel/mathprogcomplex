ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))


function main()

    #######################################################
    # Parameters
    d = 2
    max_nodes = 14



    repo = LibGit2.GitRepo(pwd()); branch = LibGit2.shortname(LibGit2.head(repo))
    date = String(Dates.format(now(), "mm_dd-HHhMM"))
    export_path = joinpath("Mosek_exports", "order$d", branch, date)
    ispath(export_path) && rm(export_path, recursive=true)


    matpower_path = joinpath("..", "data", "data_Matpower", "matpower")
    instances = OrderedSet(sort(readdir(matpower_path), by=x->parse(first(matchall(r"\d+", x)))))

    filter!(x-> parse(first(matchall(r"\d+", x))) ≤ max_nodes, instances)
    setdiff!(instances, Set(["case24_ieee_rts.m"]))
    @show collect(instances)
    i = 0
    for instance in instances
        i += 1
        info("working on $instance ($i/$(length(instances)))")

        OPFpbs = load_OPFproblems(MatpowerInput, joinpath(matpower_path, instance))
        problem_c = build_globalpb!(OPFpbs)
        problem = pb_cplx2real(problem_c)

        POP_exponents = Set{Exponent}()
        for expo in keys(problem.objective)
            push!(POP_exponents, expo)
        end
        for (ctrname, ctr) in problem.constraints
            for expo in keys(ctr.p)
                push!(POP_exponents, expo)
            end
        end
        POP_exponents = OrderedSet{Exponent}(sort(collect(POP_exponents)))

        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            d = d)


        max_cliques = get_maxcliques(relax_ctx, problem)
        momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem, max_cliques)

        mmtrel_pb = MomentRelaxation{Float64}(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)
        sdpinstance = build_SOSrelaxation(relax_ctx, mmtrel_pb)

        logpath = joinpath(export_path, instance[1:end-2])
        mkpath(logpath)
        info(logpath)
        export_SDP(sdpinstance, logpath; POP_exponents=POP_exponents)
    end
end

main()