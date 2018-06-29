ROOT = pwd()
include(joinpath(ROOT,"src_PowSysMod", "PowSysMod_body.jl"))
include(joinpath(ROOT,"src_SOShierarchy", "SOShierarchy.jl"))
include(joinpath(ROOT,"dev", "get_cliques.jl"))



function build_and_solve_SDPrelax_GOC(data_path::String, folder::String, scenario::String, order::Int64,output_dir::String)
    folder_path = joinpath(data_path, folder)
    instance_path = joinpath(folder_path, scenario)
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

    # init_point = Point()
    # init_point_real = cplx2real(init_point)
    #
    # ## Exporting real problem
    amplexportpath = joinpath("knitro_runs", "$(folder[9:end])_$(scenario)")
    #
    # my_timer = @elapsed export_to_dat(pb_global_real, amplexportpath, init_point_real)
    # input_dat1, input_dat2 = joinpath(amplexportpath,"real_minlp_instance_noSmax.dat"), joinpath(amplexportpath,"real_minlp_precond_cstrs_noSmax.dat")
    #
    # !ispath(output_dir) && mkpath(output_dir)
    #
    # ## Reading *.dat
    # problem_C, point = import_from_dat(input_dat1, precondcstrspath=input_dat2)
    #
    # ## Converting problem to real
    # problem = pb_cplx2real(problem_C)
    pt_knitro, pt_GOC = read_Knitro_output(amplexportpath, problem)
    obj = get_objective(problem, pt_knitro)

    ## Build relaxation
    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        issparse=true,
                                        d = order)

    # max_cliques = get_maxcliques(relax_ctx, problem)
    max_cliques = get_cliques(problem)

    ########################################
    # Calcul des matrices de moment
    momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem, max_cliques)
    mmtrel_pb = MomentRelaxation{Float64}(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)

    sdpinstance = build_SDPInstance(relax_ctx, mmtrel_pb)

    # ispath(output_dir) && rm(output_dir, recursive=true)
    # mkpath(output_dir)
    export_SDP(sdpinstance, output_dir)

    sdp_instance = read_SDPInstance(output_dir)

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



    originalSTDOUT = STDOUT
    outlog = open(joinpath(amplexportpath, "relaxSDP_$(folder)_$(scenario).log"), "w")
    redirect_stdout(outlog)

    primobj, dualobj = solve_mosek(sdp::SDP_Problem, primal, dual)

    close(outlog)
    redirect_stdout(originalSTDOUT)

    return scenario => obj, primobj, dualobj
end

function read_args(ARGS)
    if length(ARGS)!=2
        error("First argument must be data_path for example ..\..\data\data_GOC
            Second argument must be instance folder for example Phase_0_IEEE14
        ")
    else
        data_path = ARGS[1]
        folder = ARGS[2]
    end
    return data_path, folder
end

data_path, folder = read_args(ARGS)
folder_path = joinpath(data_path, folder)

# scenarios = ["scenario_$i" for i in 77:100]
scenarios = sort(filter(x->!ismatch(r"\.", x), readdir(folder_path)), by=x->parse(split(x, "_")[2]))
# println(scenarios)
nb_scenarios = length(scenarios)

dat_paths = [joinpath(pwd(),"knitro_runs", "$(folder[9:end])_$scenario", "real_minlp_instance_noSmax.dat") for scenario in scenarios]
output_dirs = [joinpath(pwd(),"knitro_runs", "$(folder[9:end])_$scenario") for scenario in scenarios]
order = 1


println("----------> Start para jobs")

# r = build_and_solve_GOC(folder, scenarios[1])
results = pmap(build_and_solve_SDPrelax_GOC, [data_path for i=1:nb_scenarios],[folder for i=1:nb_scenarios], scenarios, [1 for i=1:nb_scenarios], output_dirs)

f = open("gaps_$(folder[9:end]).csv","w")
write(f, "scenario;Upper bound; Lower bound (primalobj) ; Lowerbound (dualobj) ; Gap\n")

for (scenario, objs) in results
    gap = abs((objs[1]-objs[3])/objs[1])
    write(f,"$scenario;$(objs[1]);$(objs[2]);$(objs[3]);$gap\n")
end

close(f)
