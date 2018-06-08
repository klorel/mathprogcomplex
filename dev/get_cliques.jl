###

include(joinpath(pwd(),"src_PowSysMod", "PowSysMod_body.jl"))


function convert_mipb_to_pb!(pb::Problem)
    variables = pb.variables
    for var in variables
        println(var)
        if var[2] == Bool
            variable = Variable(var[1], Real)
            add_constraint!(pb, "Bin_$(var[1])", variable^2 - variable == 0)
            var[2] = Real
        end
    end
    return pb
end

V1 = Variable("VOLT_1",Complex)
V2 = Variable("VOLT_2",Complex)
b = Variable("b", Bool)
p_obj = V1*conj(V2)+2*abs2(V1)+1
p_ctr1 = abs2(V1)
p_ctr2 = V1*conj(V2)+abs2(V1)

problem_poly=Problem()
add_variable!(problem_poly,V1)
add_variable!(problem_poly,V2)
add_variable!(problem_poly,b)

add_constraint!(problem_poly,"ctr1",0.95^2<< p_ctr1 <<1.05^2)
add_constraint!(problem_poly,"ctr2",p_ctr2==0)
add_constraint!(problem_poly,"testb", (b + V2 )==0)

set_objective!(problem_poly,p_obj)

print(problem_poly)

PB_poly_real = pb_cplx2real(problem_poly)

print(convert_mipb_to_pb!(PB_poly_real))


function get_cliques(pb_real::Problem)
    cliques = SortedDict{String, SortedSet{Variable}}()
    vars = pb_real.variables
    for var in vars
        splits = split(var[1],'_')
        scen = splits[1]
        if scen == "BinVolt"
            splits[3] == basecase_scenario_name() ? scen = splits[5] : scen = splits[3]
        end
        clique ="clique-$scen"
        if !haskey(cliques, clique)
            cliques[clique] = SortedSet{Variable}()
        end
        push!(cliques[clique], Variable(var[1], var[2]))
    end
    return cliques
end


# folder = "Phase_0_IEEE14_1Scenario"
# # folder = "Phase_0_IEEE14"
# # folder = "Phase_0_RTS96"
# scenario = "scenario_1"
#
# instance_path = joinpath(pwd(),"..","data", "data_GOC", folder, scenario)
#
#
# raw = "powersystem.raw"
# gen = "generator.csv"
# con = "contingency.csv"
# rawfile = joinpath(instance_path,raw)
# genfile = joinpath(instance_path, gen)
# contfile = joinpath(instance_path, con)
# OPFpbs = load_OPFproblems(rawfile, genfile, contfile)
# introduce_Sgenvariables!(OPFpbs)
# ## Bulding optimization problem
# pb_global = build_globalpb!(OPFpbs)
# pb_global_real = pb_cplx2real(pb_global)
# cliques = get_cliques(pb_global_real)
# println(cliques)
#
# println(length(cliques["C-$(basecase_scenario_name())"]))
#
# println(length(cliques["C-Scen1"]))
