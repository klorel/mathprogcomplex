###Test lignes impedance nulle

ROOT = pwd()

cd("dev")
include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))

typeofinput = GOCInput

instances_folder = joinpath(pwd(),"..", "data_GOC")
folder = "Phase_0_Infeas179"
scenario = "scenario_1"

instance_path = joinpath(instances_folder,folder,scenario)


OPFpbs = load_OPFproblems(typeofinput, instance_path)

(typeofinput != GOCInput) || introduce_Sgenvariables!(OPFpbs)

## Bulding optimization problem
pb_global = build_globalpb!(OPFpbs)

(typeofinput != GOCInput) || (init_point = solution_point(instance_path))

feas,ctr = get_minslack(pb_global, init_point)
##Sorig à mettre dans point !
obj = get_objective(pb_global, init_point)


###test sur IEEE14 modifie
folder = "Phase_0_IEEE14"
scenario = "scenario_1"
instance_path = joinpath(instances_folder,folder,scenario)

#################
#READ OPF
##################
power_data,generator_data,contingency_data = getdataGOC(instance_path)
##info in bus
bus,index = read_data_bus(power_data)
##info in generator
gendata_dict = generator_data_to_dict(generator_data,index)
##info in generator
bus = add_generator_data!(power_data,gendata_dict,bus,index)
##info in branch
link = read_branch_data(power_data,index)


####################
#Modif
for (line, dict) in link
    if line.orig == "BUS_4" && line.dest == "BUS_7"
        data = dict["Lineπ_transfo_BL"]
        branch_id = "BL"
        link[line] = SortedDict(nullimpedance_withtransformer_name(branch_id) => GOCNullImpedance_withtransformer(line, data.id, data.susceptance, data.transfo_ratio, data.transfo_phase, data.power_magnitude_max))
    end
end
######################
#info basecase
ds = DataSource(bus,link)
node_linksin, node_linksout = SortedDict{String, SortedSet{Link}}(), SortedDict{String, SortedSet{Link}}()
node_vars = SortedDict{String, SortedDict{String, Variable}}()
link_vars = SortedDict{Link, SortedDict{String, Variable}}()
gs = GridStructure("BaseCase", node_linksin, node_linksout)
node_formulations = SortedDict{String, SortedDict{Tuple{Type, String}, Symbol}}()
link_formulations = SortedDict{Link, SortedDict{Tuple{Type, String}, Symbol}}()
mp = MathematicalProgramming(node_formulations, link_formulations, node_vars,link_vars)
##read scenarios
OPFpbs = scenarios_data(ds, gs, mp, power_data, contingency_data,index)

#################
#LOAD OPF
##################
for scenario in keys(OPFpbs)
  node_linksin, node_linksout = SortedDict{String, SortedSet{Link}}(), SortedDict{String, SortedSet{Link}}()

  for (linkname, _) in OPFpbs[scenario].ds.link
    if !haskey(node_linksin, linkname.dest)
      node_linksin[linkname.dest] = SortedSet{Link}()
    end
    union!(node_linksin[linkname.dest], [linkname])

    if !haskey(node_linksout, linkname.orig)
      node_linksout[linkname.orig] = SortedSet{Link}()
    end
    union!(node_linksout[linkname.orig], [linkname])
  end

  OPFpbs[scenario].gs.node_linksin = node_linksin
  OPFpbs[scenario].gs.node_linksout = node_linksout

  node_vars=SortedDict{String, SortedDict{String, Variable}}(busname => SortedDict{String, Variable}() for busname in keys(OPFpbs[scenario].ds.bus))
  link_vars=SortedDict{Link, SortedDict{String, Polynomial}}(link => SortedDict{String, Polynomial}() for link in keys(OPFpbs[scenario].ds.link))

  node_formulations = SortedDict{String, SortedDict{String, Symbol}}()
  link_formulations = SortedDict{Link, SortedDict{String, Symbol}}()
  for bus in keys(OPFpbs[scenario].ds.bus)
    node_formulations[bus] = SortedDict{String, Symbol}(elem => :NbMinVar for elem in keys(OPFpbs[scenario].ds.bus[bus]))
  end
  for link in keys(OPFpbs[scenario].ds.link)
    link_formulations[link] = SortedDict{String, Symbol}(elem => :NbMinVar for elem in keys(OPFpbs[scenario].ds.link[link]))
  end

  OPFpbs[scenario].mp = MathematicalProgramming(node_formulations, link_formulations, node_vars, link_vars)
end


## Bulding optimization problem
pb_global = build_globalpb!(OPFpbs)

println(keys(pb_global.variables))

println(pb_global.constraints["Scen1_4-7_NullImp_transfo_BL_$(get_NullImpVolt_cstrname())"])
