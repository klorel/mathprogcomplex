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

        ## Setting problem for simple OPF: linear objetive function
        for bus in get_buses(OPFpbs, contingency)
            for (elemid, elem) in ds.bus[bus]
                if typeof(elem) == GOCGenerator
                    for deg in keys(elem.dict_obj_coeffs)
                        deg ≤ 1 || delete!(elem.dict_obj_coeffs, deg)
                    end
                end
            end
        end

        OPFpb = build_Problem!(OPFpbs, contingency)

        cur_export = joinpath(export_path, "$(scenario)_$(contingency)")
        ispath(cur_export) && rm(cur_export, recursive=true)
        mkpath(cur_export)
        warn("Exporting to $cur_export")
        export_to_dat(OPFpb, cur_export)
    end

end