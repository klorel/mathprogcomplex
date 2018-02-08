# using ComplexModeler, Modeler, Poly,
using LightXML
# using ComplexModeler
ROOT=pwd()
include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))

filename = joinpath("instances", "IIDM", "pfLille_10postes.xml")
OPFpbs = build_OPFproblems(IIDMInput, [filename])

scenario = "BaseCase"
pb = build_Problem!(OPFpbs, scenario)
# print(pb)

println("-----------------------------------------------------------------------------")
voltage = Point()
for (busname, buselems) in OPFpbs["BaseCase"].ds.bus
  voltage[OPFpbs["BaseCase"].mp.node_vars[busname]["Volt"]] = buselems["MeasureVolt"].V
end

## Testing nodal power balance
for cstrname in keys(pb.constraints)
  if ismatch(r"Cstr_", cstrname)
    delete!(pb.constraints, cstrname)
  end
end

cstrs = get_constraints(pb, voltage)

slacks = get_slacks(pb, voltage)
println("minslack : ", get_minslack(pb, voltage))

# ordered_slacks = sort([ [k v] for (k, v) in slacks], by=x->abs(x[2]))
# for k in sort(collect(keys(slacks)), by=x->x[1], rev=true)
#   println("→ $k   \t $(slacks[k]) \t$(pb.constraints[k[1]].lb) < \t $(cstrs[k]) < \t$(pb.constraints[k[1]].ub)")
# end

## Power at orig and dest of lines.
# line_power = Dict{Link, Array{Complex128}}()
# for link in get_links(OPFpbs, "BaseCase")
#   link_elem = OPFpbs["BaseCase"].ds.link[link]
#   length(link_elem) == 1 || warn("link $link has several lines.")
#   for (line_name, line) in link_elem
#     line_power[link] = [line.S_orig line.S_dest]
#   end
# end

## Testing nodal power balance
ds = OPFpbs["BaseCase"].ds
gs = OPFpbs["BaseCase"].gs
mc = OPFpbs["BaseCase"].mc

pb_nodal_balance = Problem()
for (busname, bus_elems) in ds.bus
  bus_elems_formulations = mc.node_formulations[busname]
  bus_elems_var = gs.node_vars[busname]

  cstrname, Snode, lb, ub = get_Snodal(busname, bus_elems, bus_elems_formulations, bus_elems_var)
  S_balance = Snode + get_Slinks_in(busname, ds, gs, mc) + get_Slinks_out(busname, ds, gs, mc)

  if cstrname == "_"
    cstrname = "LOAD_"
  end
  add_constraint!(pb_nodal_balance, cstrname*string(busname), lb << S_balance << ub)
end

get_constraints(pb_nodal_balance, voltage)
slacks = get_slacks(pb_nodal_balance, voltage)
minslack, _ = get_minslack(pb_nodal_balance, voltage)
println("minslack on nodal powerbalance: $minslack")
# for k in sort(collect(keys(slacks)), by=x->x[1], rev=true)
#   println("→ $k   \t $(slacks[k]) \t$(pb_nodal_balance.constraints[k[1]].lb) < \t $(cstrs[k]) < \t$(pb_nodal_balance.constraints[k[1]].ub)")
# end

## Testing line power balance
line_power = Point()
pb_line_balance = Problem()
for link in get_links(OPFpbs, "BaseCase")
  for (elementname, elem) in ds.link[link]
    # elementname = collect(keys(ds.link[link]))[1]
    # element = ds.link[link][elementname]
    elem_formulation = OPFpbs["BaseCase"].mc.link_formulations[link][elementname]
    link_vars = OPFpbs["BaseCase"].gs.link_vars[link]
    Sor = Sorig(element, link, elementname, elem_formulation, link_vars)
    Sde = Sdest(element, link, elementname, elem_formulation, link_vars)

    add_constraint!(pb_line_balance, "", element.S_orig << Sor << element.S_orig)
    evaluate(Sor, voltage)
    evaluate(Sde, voltage)

    element.S_dest
  end
end

link = collect(keys(ds.link))[1]

export_to_dat(pb, "my$instance.dat")
sort_dat("my$instance.dat", "my$(instance)_sorted.dat")
# sort_dat("$instance.dat", "$(instance)_sorted.dat")

err, dict = compare_dat("my$instance.dat", filename_dat)

println("error = $err, $(length(dict))")
