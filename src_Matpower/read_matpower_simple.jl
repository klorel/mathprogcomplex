type MatpowerSimpleInput <: AbstractInput end


"""
    read_input(input_type::T, instance_path::String) where T<:Type{MatpowerSimpleInput}

Read instance in `instance_path` depending on `input_type`.\n
Return a structure OPFProblems.
"""

function read_input(intput_type::T, instance_path::String) where T<:Type{MatpowerSimpleInput}
  data = load_matpower(instance_path)

  # DataStructure and Gridstructure data:
  bus = SortedDict{String, SortedDict{String, Any}}()
  link = SortedDict{Link, SortedDict{String, Any}}()

  bus_id_line=SortedDict{Int, Int}()
  bus_id_name=SortedDict{Int, String}()

  checkfor(data, 2, "mpc.baseMVA")
  baseMVA = data[2,3]

  ## Building bus load and shunt information
  i_debut, i_fin = find_numarray(1, data)
  bustype = data[i_debut:i_fin, 2]                # whether bus is 1:"PQ" (generators are not accounted for) or 2:"PV"
  checkfor(data, i_debut-1, "mpc.bus")
  for i=i_debut:i_fin
    id = i-i_debut+1
    busname = bus_name(id)
    ## Adding MatpowerVolt structure (for each bus)
    if !haskey(bus, busname)
        bus[busname] = SortedDict{String, Any}()
    end
    bus[busname][volt_name()] = MatpowerVolt(busname, id, data[i,13], data[i,12])

    ## Matpower Load
    load = data[i,3] + im*data[i,4]
    if load != 0
      bus[busname]["Load"] = MatpowerLoad(load)
    end

    ## Matpower Shunt
    shunt = data[i,5] + im*data[i,6]
    if shunt != 0
      bus[busname]["Shunt"] = MatpowerShunt(shunt)
    end

    bus_id_line[data[i,1]] = id
    bus_id_name[data[i,1]] = busname
  end

  ## Adding bus generator information
  gen2bus = SortedDict{Int, Int}()
  line2busgen = SortedDict{Int, Tuple{Int, Int}}()

  i_debut, i_fin = find_numarray(i_fin+1, data)
  checkfor(data, i_debut-1, "mpc.gen")

  genind, prevgen = 0, 0
  for i=i_debut:i_fin
    gen2bus[i-i_debut+1] = bus_id_line[data[i,1]]
    busname = bus_id_name[data[i,1]]

    S_min = S_max = 0
    if data[i, 1] == prevgen
      genind += 1
      S_min = bus[busname]["Gen"].power_min
      S_max = bus[busname]["Gen"].power_max
    else
      prevgen = data[i, 1]
      genind = 1
    end
    line2busgen[i-i_debut+1] = (data[i, 1], genind)

    if data[i, 8] > 0 #generator is on
      S_min += data[i,10] + im*data[i,5] #Smin = Pmin + i Qmin
      S_max += data[i,9] + im*data[i,4] #Smax = Pmax + i Qmax
    end
    bus[busname]["Gen"] = MatpowerGenerator(S_min, S_max, -1, [-1],true)
  end

  ## building link information
  i_debut, i_fin = find_numarray(i_fin+1, data)
  checkfor(data, i_debut-1, "mpc.branch")
  for i=i_debut:i_fin
    linkname = Link(bus_id_name[data[i,1]], bus_id_name[data[i,2]])
    if data[i, 11] != 0 #Breaker is functional
      rs, xs, bc = data[i,3:5]
      τ, θ = data[i, 9:10]
      if !haskey(link, linkname)
        link[linkname] = SortedDict{String, Any}()
      end

      nb_elem = length(link[linkname])
      link[linkname]["LinkMP_$(nb_elem+1)"] = MatpowerLine_π(baseMVA, rs, xs, bc, τ, θ)
    end
  end

  i_debut, i_fin = find_numarray(i_fin+1, data)
  if data[i_debut-1, 1] == "mpc.areas"
    warn("Not loading mpc.areas data.")
    i_debut = i_fin+1
    i_debut, i_fin = find_numarray(i_fin+1, data)
  end

  bus_line_id = SortedDict([val=>key for (key,val) in bus_id_line])

  ## Adding generator cost information
  checkfor(data, i_debut-1, "mpc.gencost")
  genind, cur_genind = 0, data[i_debut, 1]
  for i=i_debut:i_fin
    buslineid, _ = line2busgen[i-i_debut+1]
    busid = bus_id_line[buslineid]

    busname = "BUS_$busid"

    cost_degree = data[i, 4]
    cost_coeffs = data[i, 5:(5+cost_degree-1)]
    if haskey(bus[busname], "Gen")
      bus[busname]["Gen"].cost_degree = cost_degree
      bus[busname]["Gen"].cost_coeffs = cost_coeffs
    else
      warn("Gen not found at $busname")
    end
  end

  ## Removing generators from PQ buses
  for (busname, bus_elems) in bus
    i = bus_elems["Volt"].busid
    if bustype[i] == 1
      for (elemid, elem) in bus_elems
        if typeof(elem) == MatpowerGenerator
          delete!(bus_elems, elemid)
        end
      end
    end
  end

  ## Adding null generator to reference bus
  refind = find(bustype .== 3)
  length(refind) == 1 || warn("/! refind = $refind")
  bus["BUS_$(refind[1])"]["Gen_reference"] = MatpowerGenerator(0,0,3,[0 0 0],true)

  ds = DataSource(bus, link)

  ## Building GridStructure
  node_linksin = node_linksout = SortedDict{String, SortedSet{Link}}()
  node_vars = SortedDict{String, SortedDict{String, Variable}}()
  link_vars = SortedDict{Link, SortedDict{String, Variable}}()
  gs = GridStructure(basecase_scenario_name(), node_linksin, node_linksout)

  node_formulations = SortedDict{String, SortedDict{String, Symbol}}()
  link_formulations = SortedDict{Link, SortedDict{String, Symbol}}()
  mp = MathematicalProgramming(node_formulations, link_formulations, node_vars, link_vars)
  return OPFProblems(basecase_scenario_name()=>Scenario(ds, gs, mp))
end
