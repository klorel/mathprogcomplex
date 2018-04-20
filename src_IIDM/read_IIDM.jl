type IIDMInput <: AbstractInput end

## Read inputs
function read_input(intput_type::T, instance_path::String) where T<:Type{IIDMInput}
  data = parse_file(instance_path)

  ## Cllecting bus information
  bus = SortedDict{String, SortedDict{String, Any}}()
  ind_bus = 1
  for substation in root(data)["substation"]
    substationname = attribute(substation, "id")
    busind = ind_bus
    ind_bus += 1

    voltagelvl = substation["voltageLevel"][1]
    length(substation["voltageLevel"])==1 || warn("Substation $substationname has several voltageLevel, accounting for $(attribute(voltagelvl, "id")) only")

    busname = attribute(voltagelvl, "id")
    length(voltagelvl["busBreakerTopology"][1]["bus"]) == 1 || warn("Substation has several buses, accounting for $busname only..")
    basekV = parse(attribute(voltagelvl, "nominalV"))
    minV = parse(attribute(voltagelvl, "lowVoltageLimit"))
    maxV = parse(attribute(voltagelvl, "highVoltageLimit"))

    busBreakerTopology = voltagelvl["busBreakerTopology"][1]["bus"][1]
    V, θ = parse(attribute(busBreakerTopology, "v")), parse(attribute(busBreakerTopology, "angle"))*π/180

    bus[busname] = SortedDict{String, Any}()
    bus[busname]["Volt"] = IIDMVolt(busname, busind, substationname, basekV, minV, maxV)
    bus[busname]["MeasureVolt"] = IIDMVoltMeasure(V*exp(im*θ), NaN)

    for cur_load in voltagelvl["load"]
      load_id = attribute(cur_load, "id")
      load0 = parse(attribute(cur_load, "p0")) + im*parse(attribute(cur_load, "q0"))
      load = parse(attribute(cur_load, "p")) + im*parse(attribute(cur_load, "q"))

      abs(load-load0) < 1e-5 || warn("$busname $load_id : load - load0 = $(load-load0)")

      bus[busname]["LOAD_$load_id"] = IIDMLoad(load_id, load)
    end

    for cur_gen in voltagelvl["generator"]
      gen_id = attribute(cur_gen, "id")
      minP = parse(attribute(cur_gen, "minP"))
      maxP = parse(attribute(cur_gen, "maxP"))
      targetV = parse(attribute(cur_gen, "targetV"))
      volt_regulator_on = parse(attribute(cur_gen, "voltageRegulatorOn"))

      attribute(cur_gen["property"][1], "value") == "FUEL" || warn("generator $gen_id of $busname not of expected type FUEL.")
      capabilityCurve = cur_gen["reactiveCapabilityCurve"][1]
      minP = parse(attribute(pt, "p"))
      minQminP, maxQminP = parse(attribute(pt, "minQ")), parse(attribute(pt, "maxQ"))
      maxP = parse(attribute(pt, "p"))
      minQmaxP, maxQmaxP = parse(attribute(pt, "minQ")), parse(attribute(pt, "maxQ"))
      S_bounds = [(minP+im*minQminP, minP+im*maxQminP), (maxP+im*minQmaxP, maxP+im*maxQmaxP)]

      bus[busname]["GEN_$gen_id"] = IIDMGenerator(gen_id, minP, maxP, targetV, volt_regulator_on, S_bounds)
    end
  end


  ## Collecting link information
  link = SortedDict{Link, SortedDict{String, Any}}()
  for line in root(data)["line"]
    orig = attribute(line, "voltageLevelId1")
    dest = attribute(line, "voltageLevelId2")
    linkname = Link(orig, dest)
    lineid = attribute(line, "id")
    r = parse(Float64, attribute(line, "r"))
    x = parse(Float64, attribute(line, "x"))
    G1 = parse(Float64, attribute(line, "g1")) + im*parse(Float64, attribute(line, "b1"))
    G2 = parse(Float64, attribute(line, "g2")) + im*parse(Float64, attribute(line, "b2"))
    S1 = parse(Float64, attribute(line, "p1")) + im*parse(Float64, attribute(line, "q1"))
    S2 = parse(Float64, attribute(line, "p2")) + im*parse(Float64, attribute(line, "q2"))
    basekV_orig = bus["$orig"]["Volt"].nomV
    basekV_dest = bus["$dest"]["Volt"].nomV
    lim = parse(Float64, attribute(line["currentLimits1"][1], "permanentLimit"))

    if !haskey(link, linkname)
      link[linkname] = SortedDict{String, Any}()
    end
    link[linkname][lineid] = IIDMLine_π(r, x, G1, G2, basekV_orig, basekV_dest, lim)
    link[linkname]["Measure_$lineid"] = IIDMLine_πMeasure(S1, S2)
  end


  ds = DataSource(bus, link)

  ## Building GridStructure
  node_linksin = node_linksout = SortedDict{String, SortedSet{Link}}()
  node_vars = SortedDict{String, SortedDict{String, Variable}}()
  link_vars = SortedDict{Link, SortedDict{String, Variable}}()
  gs = GridStructure("BaseCase", node_linksin, node_linksout)

  node_formulations = SortedDict{String, SortedDict{Tuple{Type, String}, Symbol}}()
  link_formulations = SortedDict{Link, SortedDict{Tuple{Type, String}, Symbol}}()
  mp = MathematicalProgramming(node_formulations, link_formulations, node_vars, link_vars)
  return OPFProblems("BaseCase"=>Scenario(ds, gs, mp))
end
