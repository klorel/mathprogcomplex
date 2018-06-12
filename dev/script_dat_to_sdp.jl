ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))
include(joinpath(ROOT, "src_PolynomialOptim", "PolynomialOptim.jl"))

using ArgParse



function parse_commandline()
  s = ArgParseSettings()
  s.description = "Build a relaxation at the specified order of the real POP described in *.dat into sdp files."
  @add_arg_table s begin
    "dat_path"
        help = ".dat file location"
        arg_type = String
        required = true
    "sdp_path"
        help = "*.sdp output folder"
        arg_type = String
        required = true
    "d"
        help = "hierarchy relaxation order"
        arg_type = Int
        required = true
  end
  return parse_args(s)
end


function main(args)
    input_params = parse_commandline()

    input_dat = input_params["dat_path"]
    isfile(input_dat) || error("Input $input_dat is not a file.")

    output_dir = input_params["sdp_path"]
    !ispath(output_dir) && mkpath(output_dir)

    ## Reading *.dat
    problem_C, point = import_from_dat(input_dat)

    ## Converting problem to real
    problem = pb_cplx2real(problem_C)

    ## Build relaxation
    @show input_params["d"]
    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        d = input_params["d"])

    max_cliques = get_maxcliques(relax_ctx, problem)

    ########################################
    # Calcul des matrices de moment
    momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem, max_cliques)
    mmtrel_pb = MomentRelaxation{Float64}(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)

    sdpinstance = build_SOSrelaxation(relax_ctx, mmtrel_pb)

    export_SDP(sdpinstance, output_dir)
end

main(ARGS)