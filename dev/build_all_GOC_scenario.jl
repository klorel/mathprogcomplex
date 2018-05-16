using MAT
include(joinpath(pwd(),"src_PowSysMod", "PowSysMod_body.jl"))


export_path = abspath(pwd(), "hadr_export")

folder = "Phase_0_IEEE14_1Scenario"
work_dir = abspath("..", "data", "data_GOC", folder)
ispath(work_dir) || error("Not a folder $work_dir")

scenarios = SortedSet(setdiff(readdir(work_dir), ["scorepara.csv"]))

for scenario in scenarios
    cur_dir = joinpath(work_dir, scenario)
    info("working on $cur_dir")

    OPFpbs = load_OPFproblems(GOCInput, cur_dir)

    for contingency in SortedSet(keys(OPFpbs))
        ds = OPFpbs[contingency].ds
        mp = OPFpbs[contingency].mp

        print_with_color(:light_green, "$scenario -> $contingency\n")

        ## Setting problem for simple OPF: no constraint for lines
        for link in get_links(OPFpbs, contingency)
            linkformulations = mp.link_formulations[link]
            for (elemid, elem) in ds.link[link]
                linkformulations[elemid] = :NoCtr
            end
        end

        OPFpb = build_Problem!(OPFpbs, contingency)
        println(OPFpb)

        cur_export = joinpath(export_path, "$(scenario)_$(contingency)")
        ispath(cur_export) && rm(cur_export, recursive=true)
        mkpath(cur_export)
        warn("Exporting to $cur_export")
        export_to_dat(OPFpb, cur_export)
    end

end