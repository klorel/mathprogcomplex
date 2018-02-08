type MatpowerInput <: AbstractInput end



"""
    read_input(input_type::T, instance_path::String) where T<:Type{MatpowerInput}

Read instance in `instance_path` depending on `input_type`.\n
Return a structure OPFProblems.    
"""
function read_input(input_type::T, instance_path::String) where T<:Type{MatpowerInput}
  data = load_matpower(instance_path)

  # DataStructure and Gridstructure data:
  bus = Dict{String, Dict{String, Any}}()
  link = Dict{Link, Dict{String, MatpowerLine_π}}()

  bus_id_line=Dict{Int, Int}()
  bus_id_name=Dict{Int, String}()

  checkfor(data, 2, "mpc.baseMVA")
  baseMVA = data[2,3]

  ## Building bus load and shunt information
  i_debut, i_fin = find_numarray(1, data)
  bustype = data[i_debut:i_fin, 2]                # weather bus is 1:"PQ" (generators are not accounted for) or 2:"PV"
  checkfor(data, i_debut-1, "mpc.bus")
  for i=i_debut:i_fin
    id = i-i_debut+1
    busname = bus_name(id)
    ## Adding MatpowerVolt structure (for each bus)
    bus[busname] = Dict(volt_name() => MatpowerVolt(busname, id, data[i,13], data[i,12]))

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
  gen2bus = Dict{Int, Int}()
  line2busgen = Dict{Int, Tuple{Int, Int}}()

  i_debut, i_fin = find_numarray(i_fin+1, data)
  checkfor(data, i_debut-1, "mpc.gen")

  genind, prevgen = 0, data[i_debut, 1]
  for i=i_debut:i_fin
    gen2bus[i-i_debut+1] = bus_id_line[data[i,1]]
    busname = bus_id_name[data[i,1]]

    if data[i, 1] == prevgen
      genind += 1
    else
      prevgen = data[i, 1]
      genind = 1
    end
    line2busgen[i-i_debut+1] = (data[i, 1], genind)

    gen_status = true
    if data[i, 8] ≤ 0 #generator is off
      gen_status = false
    end

    S_min = data[i,10] + im*data[i,5] #Smin = Pmin + i Qmin
    S_max = data[i,9] + im*data[i,4] #Smax = Pmax + i Qmax

    bus[busname]["Gen_$genind"] = MatpowerGenerator(S_min, S_max, -1, [-1], gen_status)
  end

  ## building link information
  i_debut, i_fin = find_numarray(i_fin+1, data)
  checkfor(data, i_debut-1, "mpc.branch")
  for i=i_debut:i_fin
    linkname = Link(bus_id_name[data[i,1]], bus_id_name[data[i,2]])
    if data[i, 11] == 0
      warn("$link $(linkname.orig)⟶$(linkname.dest) breaker out of service !")
    else
      rs, xs, bc = data[i,3:5]
      τ, θ = data[i, 9:10]
      if !haskey(link, linkname)
        link[linkname] = Dict{String, MatpowerLine_π}()
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

  bus_line_id = Dict([val=>key for (key,val) in bus_id_line])

  ## Adding generator cost information
  checkfor(data, i_debut-1, "mpc.gencost")
  genind, cur_genind = 0, data[i_debut, 1]
  for i=i_debut:i_fin
    buslineid, genid = line2busgen[i-i_debut+1]
    busid = bus_id_line[buslineid]

    busname = "BUS_$busid"

    cost_degree = data[i, 4]
    cost_coeffs = data[i, 5:(5+cost_degree-1)]
    if haskey(bus[busname], "Gen_$genid")
      bus[busname]["Gen_$genid"].cost_degree == -1 || warn("read_intput(): $busname, gen $genid already has a defined generator cost.")
      bus[busname]["Gen_$genid"].cost_degree = cost_degree
      bus[busname]["Gen_$genid"].cost_coeffs = cost_coeffs
    else
      println("Gen_$genid not found at $busname")
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

  ds = DataSource(bus, link)

  ## Building GridStructure
  node_linksin = node_linksout = Dict{String, Set{Link}}()
  node_vars = Dict{String, Dict{String, Variable}}()
  link_vars = Dict{Link, Dict{String, Variable}}()
  gs = GridStructure(basecase_scenario_name(), node_linksin, node_linksout)

  node_formulations = Dict{String, Dict{Tuple{Type, String}, Symbol}}()
  link_formulations = Dict{Link, Dict{Tuple{Type, String}, Symbol}}()
  mp = MathematicalProgramming(node_formulations, link_formulations, node_vars, link_vars)
  return OPFProblems(basecase_scenario_name()=>Scenario(ds, gs, mp))
end


## Utils functions
function find_numarray(i_start, data)
  i_debut = i_start
  while !isa(data[i_debut, 1], Int)
    i_debut+=1
  end
  i_fin=i_debut
  while !isa(data[i_fin,1], SubString)
    i_fin += 1
  end
  i_debut, i_fin-1
end

checkfor(data, line_ind, name) = (data[line_ind, 1] == name) || error("Expected ", name, " at line ", line_ind, ", got ", data[line_ind,1], " instead.")

function load_matpower(filename)
  instance_name = split(filename, '.')[1]

  touch(instance_name*".temp")
  f = open(filename)
  out = open(instance_name*".temp", "w")

  # removing all ';' at end of lines
  while !eof(f)
    line = readline(f)
    if length(line) > 0 && line[1] != '%' && line[1] != 'f'
      s = split(line, ";")
      println(out, s[1])
    end
  end
  close(f)
  close(out)

  data = readdlm(instance_name*".temp")
  rm(instance_name*".temp")
  return data
end
