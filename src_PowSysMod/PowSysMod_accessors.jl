## Gridstructure accessors
function get_buses(OPFpbs::OPFProblems, scenario::String)
  return keys(OPFpbs[scenario].ds.bus)
end

function get_links(OPFpbs::OPFProblems, scenario::String)
  return keys(OPFpbs[scenario].ds.link)
end

function get_node_elems(bus::String, OPFpbs::OPFProblems, scenario::String)
  return keys(OPFpbs[scenario].ds.bus[bus])
end

function get_link_elems(link::Link, OPFpbs::OPFProblems, scenario::String)
  return keys(OPFpbs[scenario].ds.link[link])
end

function get_elem_type(elem::String, bus::String, OPFpbs::OPFProblems, scenario::String)
  return typeof(OPFpbs[scenario].ds.bus[bus][elem])
end

function get_elem_type(elem::String, link::Link, OPFpbs::OPFProblems, scenario::String)
  return typeof(OPFpbs[scenario].ds.link[link][elem])
end

function get_links_in(bus::String, gs::GridStructure)
  return haskey(gs.node_linksin, bus) ? gs.node_linksin[bus] : SortedSet{Any}()
end

function get_links_out(bus::String, gs::GridStructure)
  return haskey(gs.node_linksout, bus) ? gs.node_linksout[bus] : SortedSet{Any}()
end

function get_Var(elem_type::String, bus::String, OPFpbs::OPFProblems, scenario::String)
  return OPFpbs[scenario].gs.node_vars[bus][elem_type]
end

function get_Var(elem_type::String, link::Link, OPFpbs::OPFProblems, scenario::String)
  return OPFpbs[scenario].gs.link_vars[link][elem_type]
end

## Modeling Choice accessors
function get_elem_formulation(elem_id::String, bus::String, OPFpbs::OPFProblems, scenario::String)
  return OPFpbs[scenario].mc.node_formulations[bus][elem_id]
end

function get_elem_formulation(elem_id::Type, link::Link, OPFpbs::OPFProblems, scenario::String)
  return OPFpbs[scenario].mc.node_formulations[link][elem_id]
end
