"""
    build_POP!(typeofinput::T, OPFpbs)

Return the optimization problem depending on `typeofinput`
"""
function build_POP!(typeofinput::T, OPFpbs) where T<:Type{MatpowerInput}
    return build_globalpb!(OPFpbs)
end

function build_POP!(typeofinput::T, OPFpbs) where T<:Type{IIDMInput}
    return build_globalpb!(OPFpbs)
end

function build_POP!(typeofinput::T, OPFpbs) where T<:Type{GOCInput}
    introduce_Sgenvariables!(OPFpbs)
    return build_globalpb!(OPFpbs)
end

function exportPOP(typeofinput::T, instance_path::String, pb_global_real) where T<:Type{GOCInput}
    global_point = solution_point(instance_path)
    folder,scenario = split(instance_path, '\\')[end-1:end]
    folder = String(folder)
    scenario = String(scenario)
    amplexportpath = joinpath("..","knitro_runs", "$(folder[9:end])_$(scenario)_move")
    if !isdir(amplexportpath)
        mkpath(amplexportpath)
    end
    global_point_real = cplx2real(global_point)
    export_to_dat(pb_global_real, amplexportpath, global_point_real)
    return amplexportpath
end

function exportPOP(typeofinput::T, instance_path::String, pb_global_real) where T<:Type{MatpowerInput}
    point = Point()
    instance_name  = String(split(instance_path,'\\')[end])
    amplexportpath = joinpath("..","knitro_runs", "$(instance_name[1:end-2])_move")
    if !isdir(amplexportpath)
        mkpath(amplexportpath)
    end
    point_real = cplx2real(point)
    export_to_dat(pb_global_real, amplexportpath, point_real)
    return amplexportpath
end





"""
    build_and_solve_GOC(instance::String, scenario::String, save_folder::String)

Build and solve GOC instance in `instance`\`scenario`, save logs in `save_folder` and returns information as minslack of solution points.
"""

function build_and_solve_GOC(instance::String, scenario::String, save_folder::String)
    filenames_path = joinpath("data_GOC", instance, scenario)
    filenames = vcat(instance, scenario, readdir(filenames_path))

    ## Build global pb
    OPFpbs = load_OPFproblems(GOCInput, filenames_path)
    introduce_Sgenvariables!(OPFpbs)
    pb_global = build_globalpb!(OPFpbs)
    ## Reading GOC solution
    global_point = solution_point(filenames_path)
    ## Export problem
    date = Dates.format(now(), "yy_u_dd_HH_MM_SS")
    amplexportpath = joinpath(save_folder, "$(date)_$(instance[8:end])_$(scenario)")

    export_dat_GOC_problem(amplexportpath, OPFpbs, pb_global, global_point)

    pb_global_real = pb_cplx2real(pb_global)

    ## Knitro call and solution import
    pt_strat1, pt_strat2 = get_knitro_solutions(amplexportpath)

    ## Results:
    Data = SortedDict()
    Data["minslack_GOC"] = get_minslack(pb_global, global_point)
    Data["minslack_strat1"] = get_minslack(pb_global_real, pt_strat1)[1]
    Data["minslack_strat2"] = get_minslack(pb_global_real, pt_strat1)[2]
    Data["thread_id"] = myid()

    println("Minslack of GOC global solution: $(get_minslack(pb_global, global_point))")
    println("Minslack for strat 1 : $(get_minslack(pb_global_real, pt_strat1))")
    println("Minslack for strat 2 : $(get_minslack(pb_global_real, pt_strat2))")
    println("====X $instance, $scenario done   job id: $(myid())")

    return (instance, scenario) => Data
end
