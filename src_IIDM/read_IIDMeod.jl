type IIDMeodInput <: AbstractInput end

################################################################################
################################################################################
## Untested

## Read inputs
function read_input(intput_type::T, filenames::Array{String}) where T<:Type{IIDMeodInput}
  data_branches = readdlm(filenames[1])
  data_buses = readdlm(filenames[2])
  data_generators = readdlm(filenames[3])
  data_limits = readdlm(filenames[4])
  data_loads = readdlm(filenames[5])
  data_substations = readdlm(filenames[6])

  bus = Dict{String, Dict{String, Any}}()

  for i=1:size(data_buses, 1)
      busid = Int(data_buses[i, 1])
      substationid = Int(data_buses[i, 2])
      V = data_buses[i, 4] # pu
      θ = data_buses[i, 5] # radians
      P, Q = data_buses[i, 6:7] # MW, MVar
      basekV, minV, maxV = data_substations[substationid, 4:6] # kV, pu, pu

      busname = "BUS_$busid"
      data_buses[i, 3] == 0 || warn("Buses - $busname : parameter cc!=0 ($(data_buses[i, 3])), not handled")
      data_buses[i, 8] == 0 || warn("Buses - $busname : parameter fault!=0 ($(data_buses[i, 8])), not handled")
      data_buses[i, 9] == 0 || warn("Buses - $busname : parameter curative!=0 ($(data_buses[i, 9])), not handled")

      bus[busname] = Dict{String, Any}()
      bus_elems = bus[busname]
      bus_elems["Volt"] = IIDMVolt(busname, busid, string(substationid), basekV, minV, maxV)
      bus_elems["VoltMeasure"] = IIDMVoltMeasure(V*exp(im*θ), P+im*Q)
  end

  for i=1:size(data_loads, 2)
      busid = data_loads[i, 2]
      P, Q = data_loads[i, 4:5] #MW, MVar
      busname = "BUS_$busid"

      data_loads[i, 6] == 0 || warn("Loads - $busname : parameter fault!=0 ($(data_loads[i, 6])), not handled")
      data_loads[i, 7] == 0 || warn("Loads - $busname : parameter curative!=0 ($(data_loads[i, 7])), not handled")

      bus_elems = bus["$busname"]
      nb = count([ismatch(r"Load_", x) for x in keys(bus_elems)])
      bus_elems["Load_$(nb+1)"] = IIDMLoad(string(nb+1), P+im*Q)
  end

  for i=1:size(data_generators, 1)
      busid = data_generators[i, 2]
      minP, maxP = data_generators[i, 5:6] # MW
      volt_regulator_on, targetV = data_generators[i, 11:12] # bool, pu
      minQmaxP, minQminP, maxQmaxP, maxQminP = data_generators[i, 7:10] # MVar
      S_bounds = [(minP+im*minQminP, minP+im*maxQminP), (maxP+im*minQmaxP, maxP+im*maxQmaxP)]
      busname = "BUS_$busid"

      data_generators[i, 2] == data_generators[i, 3] || warn("Generators - $busname : 'bus' is not 'con. bus', not handled")
      data_generators[i, 15] == 0 || warn("Loads - $busname : parameter fault!=0 ($(data_generators[i, 15])), not handled")
      data_generators[i, 16] == 0 || warn("Loads - $busname : parameter curative!=0 ($(data_generators[i, 16])), not handled")

      bus_elems = bus["$busname"]
      nb = count([ismatch(r"Gen_", x) for x in keys(bus_elems)])
      bus_elems["Gen_$(nb+1)"] = IIDMGenerator(string(nb+1), minP, maxP, targetV, volt_regulator_on, S_bounds)
  end


  link = Dict{Link, Dict{String, Any}}()
  for i=1:size(data_branches, 1)
      busid1, busid2 = data_branches[i, 2:3]
      r, x, g1, g2, b1, b2 = data_branches[i, 7:12] # pu
      p1, p2, q1, q2 = data_branches[i, 16:19] # MW, MW, MVar, MVar


      linkname = Link("BUS_$busid1", "BUS_$busid2")
      if !haskey(link, linkname)
          link[linkname] = Dict{String, Any}()
      end
      link_elems = link[linkname]

      basekV_dest = basekV_orig = lim = 0

      nb = count([ismatch(r"Branch_", x) for x in keys(link_elems)])
      link_elems["IIDMLine_π_$(nb+1)"] = IIDMLine_π(r, x, g1+im*b1, g2+im*b2, basekV_orig, basekV_dest, lim)
      link_elems["IIDMLine_πMeasure_$(nb+1)"] = IIDMLine_πMeasure(p1+im*q1, p2+im*q2)
  end

  ds = DataSource(bus, link)

  ## Building GridStructure
  node_linksin = node_linksout = Dict{String, Set{Link}}()
  node_vars = Dict{String, Dict{String, Variable}}()
  link_vars = Dict{Link, Dict{Type, Variable}}()
  gs = GridStructure("BaseCase", node_linksin, node_linksout)

  node_formulations = Dict{String, Dict{Tuple{Type, String}, Symbol}}()
  link_formulations = Dict{Link, Dict{Tuple{Type, String}, Symbol}}()
  mp = MathematicalProgramming(node_formulations, link_formulations, node_vars, link_vars)
  return OPFProblems("BaseCase"=>Scenario(ds, gs, mp))
end
