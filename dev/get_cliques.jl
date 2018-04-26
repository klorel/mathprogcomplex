###

include(joinpath(pwd(),"..","src_PowSysMod", "PowSysMod_body.jl"))

folder = "Phase_0_IEEE14_1Scenario"
# folder = "Phase_0_IEEE14"
# folder = "Phase_0_RTS96"
scenario = "scenario_1"

instance_path = joinpath(pwd(),"..","..","data", "data_GOC", folder, scenario)


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
pb_global_real = pb_cplx2real(pb_global)

cliques = Dict{String, Set{Variable}}()

function get_cliques(pb_real::Problem)
    vars = pb_real.variables

    generator_buses = Set()
    for var in vars
        if ismatch(r"Sgen", var[1])
            bus_id = split(var[1],'_')[2]
            push!(generator_buses, bus_id)
        end
    end

    for var in vars
        splits = split(var[1],'_')
        scen = splits[1]
        bus_id = splits[2]
        if bus_id ∈ generator_buses
            clique ="C-bus-$bus_id)"
            if !haskey(cliques, clique)
                cliques[clique] = Set{Variable}()
            end
            push!(cliques[clique], Variable(var[1], var[2]))
        else
            clique ="C-$scen"
            if !haskey(cliques, clique)
                cliques[clique] = Set{Variable}()
            end
            push!(cliques[clique], Variable(var[1], var[2]))
        end
    end
    return cliques
end

get_cliques(pb_global_real)
