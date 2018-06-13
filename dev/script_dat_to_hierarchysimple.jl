ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))

using ArgParse

function parse_commandline()
  s = ArgParseSettings()
  s.description = "Load the *.sdp files at the given location and start Mosel solve."
  @add_arg_table s begin
    "instance_name"
        help="instance name (without .dat or .m extension)"
        arg_type = String
        required = true
    "d"
        help="relaxation order"
        arg_type = Int
        required = true
    "logpath"
        help = "folder storing all output"
        arg_type = String
        required = true
  end
  return parse_args(s)
end


function main(args)
    input_params = parse_commandline()

    instance_name = input_params["instance_name"]
    d = input_params["d"]
    hierarchykind = :Real
    symmetries = DataType[]
    logpath = input_params["logpath"]


    ## Run funcitons once for precompilation...
    instance_dat = joinpath("..", "data", "data_Matpower", "matpower_QCQP", "WB2.dat")
    problem_C, point = import_from_dat(instance_dat)
    problem = pb_cplx2real(problem_C)
    relax_ctx = set_relaxation(problem; hierarchykind=hierarchykind,
                                        d=1,
                                        params = Dict(:opt_outlev=>0,
                                                      :opt_outmode=>0)
    run_hierarchy(problem, relax_ctx, logpath; save_pbs=false)


    ## Build real problem
    instance_dat = joinpath("..", "data", "data_Matpower", "matpower_QCQP", instance_name*".dat")
    problem_C, point = import_from_dat(instance_dat)
    problem = pb_cplx2real(problem_C)

    ## Build relaxation context
    relax_ctx = set_relaxation(problem; hierarchykind=hierarchykind,
                                        d=d,
                                        symmetries=symmetries,
                                        params = Dict(:opt_outlev=>1,
                                                      :opt_outmode=>1,
                                                      :opt_outcsv=>1,
                                                      :opt_outname=>joinpath(logpath, "momentsos.log"),
                                                      :opt_outcsvname=>joinpath(logpath, "momentsos_solve.csv")))
    relax_ctx.relaxparams[:pb_name] = instance_name

    run_hierarchy(problem, relax_ctx, logpath; save_pbs=false)

    final_output(relax_ctx)

end

main(ARGS)