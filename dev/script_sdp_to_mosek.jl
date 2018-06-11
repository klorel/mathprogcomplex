ROOT = pwd()
include(joinpath(ROOT, "src_SOShierarchy", "SOShierarchy.jl"))


using ArgParse

function parse_commandline()
  s = ArgParseSettings()
  s.description = "Load the *.sdp files at the given location and start Mosel solve."
  @add_arg_table s begin
    "sdp_path"
        help = "*.sdp input folder"
        arg_type = String
        required = true
    "msk_log"
        help = "States whether to save Mosek log - optional"
        arg_type = Bool
        default = false
  end
  return parse_args(s)
end


function main(args)
    input_params = parse_commandline()

    work_dir = input_params["sdp_path"]
    ispath(work_dir) || error("Input $work_dir is not a valid path.")

    ## Reading *.dat
    sdp_instance = read_SDPInstance(work_dir)

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

    primal = SortedDict{Tuple{String,String,String}, Float64}()
    dual = SortedDict{Tuple{String, String, String}, Float64}()

    logname = ""
    input_params["msk_log"] && (logname = joinpath(input_params["sdp_path"], "msk_out.log"))

    primobj, dualobj = solve_mosek(sdp::SDP_Problem, primal, dual, logname=logname)

    input_params["msk_log"] && info("Mosek log saved at $logname.")
end

main(ARGS)