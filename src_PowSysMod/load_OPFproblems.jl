"""
    load_OPFproblems(input_type, instance_path::String)

Read provided files in `instance_path` function of `input_type` and create a structure ds of type DataSource for each scenario. Also create empty structures gs and mp for each scenario
For each scenario, treat data in ds to complete gs : from bus and links in ds, generate node_linksin and node_linksout
For each scenario, complete mp structure
Return a dictonary Scenario => Scenario structure


# Arguments
- `input_type::T` : type of entry format for example MatpowerInput, GOCInput
- `instance_path::String` : path to the instance to treat, provided the instance path must be associated to an instance of format input_type

# Output
- `OPFpbs`


# Examples
```
julia > OPFpbs = load_OPFproblems(MatpowerInput, "instances\\matpower\\WB2.m")
julia > OPFpbs
Dict{String,Scenario} with 1 entry:
  "BaseCase" => Scenario
```
"""

function load_OPFproblems(input_type, instance_path::String)
  # raw = "powersystem.raw"
  # gen = "generator.csv"
  # con = "contingency.csv"
  #
  # rawfile = joinpath(instance_path,raw)
  # genfile = joinpath(instance_path, gen)
  # confile = joinpath(instance_path, con)
  #
  # OPFpbs = read_GOCfiles(rawfile, genfile, confile)

  OPFpbs = read_input(input_type, instance_path)
  # baseMVA, bus, link, node_elems, link_elems, elem_to_typekey = read_input(input_type, filenames)
  # ds = DataSource(baseMVA, bus, link)

  for scenario in keys(OPFpbs)
    node_linksin, node_linksout = Dict{String, Set{Link}}(), Dict{String, Set{Link}}()

    for (linkname, _) in OPFpbs[scenario].ds.link
      if !haskey(node_linksin, linkname.dest)
        node_linksin[linkname.dest] = Set{Link}()
      end
      union!(node_linksin[linkname.dest], [linkname])

      if !haskey(node_linksout, linkname.orig)
        node_linksout[linkname.orig] = Set{Link}()
      end
      union!(node_linksout[linkname.orig], [linkname])
    end

    OPFpbs[scenario].gs.node_linksin = node_linksin
    OPFpbs[scenario].gs.node_linksout = node_linksout

    node_vars=Dict{String, Dict{String, Variable}}(busname => Dict{String, Variable}() for busname in keys(OPFpbs[scenario].ds.bus))
    link_vars=Dict{Link, Dict{String, Variable}}(link => Dict{String, Variable}() for link in keys(OPFpbs[scenario].ds.link))

    node_formulations = Dict{String, Dict{String, Symbol}}()
    link_formulations = Dict{Link, Dict{String, Symbol}}()
    for bus in keys(OPFpbs[scenario].ds.bus)
      node_formulations[bus] = Dict{String, Symbol}(elem => :NbMinVar for elem in keys(OPFpbs[scenario].ds.bus[bus]))
    end
    for link in keys(OPFpbs[scenario].ds.link)
      link_formulations[link] = Dict{String, Symbol}(elem => :NbMinVar for elem in keys(OPFpbs[scenario].ds.link[link]))
    end

    OPFpbs[scenario].mp = MathematicalProgramming(node_formulations, link_formulations, node_vars, link_vars)
  end

  return OPFpbs
end


function load_OPFproblems(rawfile::String, genfile::String, confile::String)
  OPFpbs = read_GOCfiles(rawfile, genfile, confile)

  for scenario in keys(OPFpbs)
    node_linksin, node_linksout = Dict{String, Set{Link}}(), Dict{String, Set{Link}}()

    for (linkname, _) in OPFpbs[scenario].ds.link
      if !haskey(node_linksin, linkname.dest)
        node_linksin[linkname.dest] = Set{Link}()
      end
      union!(node_linksin[linkname.dest], [linkname])

      if !haskey(node_linksout, linkname.orig)
        node_linksout[linkname.orig] = Set{Link}()
      end
      union!(node_linksout[linkname.orig], [linkname])
    end

    OPFpbs[scenario].gs.node_linksin = node_linksin
    OPFpbs[scenario].gs.node_linksout = node_linksout

    node_vars=Dict{String, Dict{String, Variable}}(busname => Dict{String, Variable}() for busname in keys(OPFpbs[scenario].ds.bus))
    link_vars=Dict{Link, Dict{String, Variable}}(link => Dict{String, Variable}() for link in keys(OPFpbs[scenario].ds.link))

    node_formulations = Dict{String, Dict{String, Symbol}}()
    link_formulations = Dict{Link, Dict{String, Symbol}}()
    for bus in keys(OPFpbs[scenario].ds.bus)
      node_formulations[bus] = Dict{String, Symbol}(elem => :NbMinVar for elem in keys(OPFpbs[scenario].ds.bus[bus]))
    end
    for link in keys(OPFpbs[scenario].ds.link)
      link_formulations[link] = Dict{String, Symbol}(elem => :NbMinVar for elem in keys(OPFpbs[scenario].ds.link[link]))
    end

    OPFpbs[scenario].mp = MathematicalProgramming(node_formulations, link_formulations, node_vars, link_vars)
  end

  return OPFpbs
end
