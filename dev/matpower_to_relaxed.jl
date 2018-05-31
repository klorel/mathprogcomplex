ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))


function main()

    export_path = joinpath("Mosek_exports", "order2")
    ispath(export_path) && rm(export_path, recursive=true)

    matpower_path = joinpath("..", "data", "data_Matpower", "matpower")

    instances = OrderedSet(sort(readdir(matpower_path), by=x->parse(first(matchall(r"\d+", x)))))

    @show collect(instances)
    i = 0
    for instance in instances
        i += 1
        if parse(first(matchall(r"\d+", instance))) < 60
            info("working on $instance ($i/$(length(instances)))")

            OPFpbs = load_OPFproblems(MatpowerInput, joinpath(matpower_path, instance))
            problem_c = build_globalpb!(OPFpbs)
            problem = pb_cplx2real(problem_c)

            relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                                d = 2)

            sdpinstance = build_relaxation(problem, relax_ctx)

            logpath = joinpath(export_path, instance[1:end-2])
            mkpath(logpath)
            export_SDP(relax_ctx, sdpinstance, logpath, indentedprint=false)
        else
            warn("jumping $instance ($i/$(length(instances)))")
        end
    end
end

main()