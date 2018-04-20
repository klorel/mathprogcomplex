ROOT = pwd()
include(joinpath(ROOT, "src_PowSysMod", "PowSysMod_body.jl"))

# folder = "Phase_0_IEEE14"
# folder = "Phase_0_Feas179"
# folder = "Phase_0_RTS96"
folder = "Phase_0_IEEE14_1Scenario"
scenario = "scenario_1"
filenames = [folder, scenario]
filenames_path = joinpath("data_GOC", folder, scenario)
instance_path = filenames_path

## Build global pb
OPFpbs = load_OPFproblems(GOCInput, instance_path)
introduce_Sgenvariables!(OPFpbs)
pb_global = build_globalpb(OPFpbs)

## Build solution
global_point = solution_point(instance_path)

## Add binary variables value
# bin_point = create_bin_point(OPFpbs, global_point)

println("minslack of GOC global solution: $(get_minslack(pb_global, global_point))")

# evaluate(pb.constraints["Scen1_3_Gen1_CC_Supper"].p, global_point)
# evaluate(pb.constraints["Scen1_3_Gen1_CC_Slower"].p, global_point)

get_minslack(pb_global, global_point)

pb_sansCC = Problem()
pb_sansCC.variables = deepcopy(pb_global.variables)
pb_sansCC.objective = deepcopy(pb_global.objective)
for (cstrname, cstr) in pb_global.constraints
    if !ismatch(r"Pgen|Qgen", cstrname)
        add_constraint!(pb_sansCC, cstrname, cstr)
    end
end
println("minslack of GOC solution without coupling constraints: $(get_minslack(pb_sansCC, global_point))")

## Export problem
amplexportpath = joinpath("knitro_runs", "$(folder[9:end])_$(scenario)_move")

pb_global_real = pb_cplx2real(pb_global)
global_point_real = cplx2real(global_point)

export_to_dat(pb_global_real, amplexportpath, global_point)


## Knitro call and solution import
run_knitro(amplexportpath, joinpath(pwd(), "src_ampl"))
pt_knitro, pt_GOC = read_Knitro_output(amplexportpath, pb_global_real)


length(global_point_real)
length(pt_knitro)

println("Minslack for GOC point : $(get_minslack(pb_global_real, pt_GOC))")
println("Minslack for Knitro point : $(get_minslack(pb_global_real, pt_knitro))")

println("Objective for GOC point :    $(get_objective(pb_global_real, pt_GOC))")
println("Objective for Knitro point : $(get_objective(pb_global_real, pt_knitro))")

println("||pt_GOC - global_point||_2 = ", norm(pt_GOC - global_point, 2))
println("||pt_GOC - global_point||_Inf = ", norm(pt_GOC - global_point, Inf))

# slacks = sort(collect(get_slacks(pb_global_real, pt_GOC)), by=x->real(x[2]))
# violated_cstrs = collect(filter(x->real(x[2])<-1e-4, slacks))
# println("GOC 1e-4 Violated_constraints    : $(length(violated_cstrs))")
#
# slacks = sort(collect(get_slacks(pb_global_real, pt_knitro)), by=x->real(x[2]))
# violated_cstrs = collect(filter(x->real(x[2])<-1e-4, slacks))
# println("Knitro 1e-4 Violated_constraints : $(length(violated_cstrs))")


# bc_volt_vars, bc_bin_vars, bc_prod_vars = get_splitted_Cpt(pt_knitro, "BaseCase")
# sc_volt_vars, sc_bin_vars, sc_prod_vars = get_splitted_Cpt(pt_knitro, "Scen1")
# sc_delta_var = get_delta_var(pt_knitro, "Scen1")
# haskey(pt_knitro, Variable("BaseCase_3_Gen1_Sgen", Complex))
#
#
# using Plots
#
# plot_Volt_vars(pt_knitro, keys(OPFpbs))
# # plot_Sgen_vars(pt_knitro, ["BaseCase"])
# plot_Sgen_vars(pt_knitro, keys(OPFpbs))
#
# plot_ViVj_vars(pt_knitro, ["BaseCase"])
# plot_ViVj_vars(pt_knitro, keys(OPFpbs))

# scatter!(real(ptsVbc), imag(ptsVbc), linealpha=0.1, marker=:circle, series_annotations=[text(matchall(r"", x[1])[1]) for (x,y) in bc_volt_vars])
