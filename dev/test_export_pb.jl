ROOT = pwd()
include(joinpath(ROOT, "src_PowSysMod", "PowSysMod_body.jl"))

folder = "Phase_0_IEEE14_1Scenario"
scenario = "scenario_1"
filenames = [folder, scenario]
filenames_path = joinpath("data_GOC", folder, scenario)
instance_path = filenames_path

## Build global pb
OPFpbs = load_OPFproblems(GOCInput, instance_path)
introduce_Sgenvariables!(OPFpbs)
pb_global = build_globalpb!(OPFpbs)

## Build solution
global_point = solution_point(instance_path)
println("minslack of GOC global solution: $(get_minslack(pb_global, global_point))")

## Export problem
pb_global_real = pb_cplx2real(pb_global)
pt_global_real = cplx2real(global_point)

# output_file = joinpath("toto", "real_minlp_instance.dat")
export_to_dat(pb_global_real, "src_ampl", pt_global_real)
# pb_global_real = pb_cplx2real(pb_global)
# pt_global_real = cplx2real(global_point)
#
# output_file = joinpath(amplexportpath, "real_minlp_instance.dat")
# export_to_dat(pb_global_real, output_file, pt_global_real)
