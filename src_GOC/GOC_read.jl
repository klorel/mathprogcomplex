type GOCInput <: AbstractInput end


##TODO
function get_baseMVA(bus::String)
    return 100
end


"""
    read_input(input_type::T, instance_path::String) where T<:Type{GOCInput}

Read instance in `instance_path` depending on `input_type`.\n
Return a structure OPFProblems.\n
For `GOCInput`, read files pscopf_data.mat, contingency.csv and generator.csv, store them in arrays and use information to fulfill a DataSource structure for the basecase scenario. Finally take into account contingencies.

"""
function read_input(input_type::T, instance_path::String) where T<:Type{GOCInput}
    #read
    power_data,generator_data,contingency_data = getdataGOC(instance_path)
    ##info in bus
    bus,index = read_data_bus(power_data)
    ##info in generator
    gendata_dict = generator_data_to_dict(generator_data,index)
    ##info in generator
    bus = add_generator_data!(power_data,gendata_dict,bus,index)
    ##info in branch
    link = read_branch_data(power_data,index)
    #info basecase
    ds = DataSource(bus,link)
    node_linksin, node_linksout = SortedDict{String, SortedSet{Link}}(), SortedDict{String, SortedSet{Link}}()
    node_vars = SortedDict{String, SortedDict{String, Variable}}()
    link_vars = SortedDict{Link, SortedDict{String, Variable}}()
    gs = GridStructure("BaseCase", node_linksin, node_linksout)
    node_formulations = SortedDict{String, SortedDict{String, Symbol}}()
    link_formulations = SortedDict{Link, SortedDict{String, Symbol}}()
    mp = MathematicalProgramming(node_formulations, link_formulations, node_vars,link_vars)
    ##read scenarios
    OPFproblems = scenarios_data(ds, gs, mp, contingency_data,index)
    return OPFproblems
end


#########################
## Utils functions
#########################
"""
    getdataGOC(instance_path)

Read 3 files (pscopf_data.mat,generator.csv,contingency.csv) in `instance_path` and returns 3 arrays of data (power_data,generator_data,contingency_data)
Return an error if one the three files does not exist in path_folder.
"""
function getdataGOC(instance_path) ##String instance+scenario
    power_data = getpowerdata(instance_path)
    generator_data = getgeneratordata(instance_path)
    contingency_data = getcontingencydata(instance_path)
    return power_data,generator_data,contingency_data
end

function getpowerdata(instance_path)
    path_folder = instance_path
    if !("pscopf_data.mat" ∈ readdir(path_folder))
        error("missing file : pscopf_data.mat must be in $path_folder")
    end
    try
        return read(matopen(joinpath(path_folder,"pscopf_data.mat")))
    catch
        warn("Not possible to read matopen(pscopf_data.mat)")
        return
    end
end

function getgeneratordata(instance_path)
    path_folder = instance_path
    if !("generator.csv" ∈ readdir(path_folder))
        error("missing file : generator.csv must be in $path_folder")
    end
    try
        return readdlm(joinpath(path_folder,"generator.csv"),',',skipstart=1)
    catch
        warn("Not possible to readdlm generator.csv")
        return
    end

end

function getcontingencydata(instance_path)
    path_folder = instance_path
    if !("contingency.csv" ∈ readdir(path_folder))
        error("missing file : contingency.csv must be in $path_folder")
    end
    try
        return readdlm(joinpath(path_folder,"contingency.csv"),',',skipstart=1)
    catch
        warn("Not possible to readdlm contingency.csv")
    end
end



#########################
function get_bus_data(bus_data, listofkeys, id)
    listofvalues = []
    for key in listofkeys
        if !haskey(bus_data, key)
            error("bus_data dictionary does not have a key $key")
        else
            append!(listofvalues, bus_data[key][id])
        end
    end
    return listofvalues
end

function get_bus_index(power_data)
    nb_bus = Int(power_data["totalbus"])
    all_busdata_key = filter(x->ismatch(r"all_busdata", x), collect(keys(power_data)))[1]
    raw_bus_data = power_data[all_busdata_key]
    nb_lines = size(raw_bus_data,1)
    if nb_lines!=nb_bus
        error("number of lines of all_busdata not equal to number of bus in power_data at key totalbus")
    end
    index = SortedDict{Int64,Int64}()
    for line in 1:nb_lines
        id_bus = Int(raw_bus_data[line,1])
        index[id_bus] = line
    end
    return index
end


"""
    read_data_bus(power_data)

Reads data in `power_data` to construct the DataSource field bus without generator data (only VOLTM, load and shunt data)
"""
function read_data_bus(power_data)
    nb_bus = Int(power_data["totalbus"])
    index = get_bus_index(power_data)
    bus_data = power_data["bus"]
    # all_busName_key = filter(x->ismatch(r"all_busName", x), collect(keys(power_data)))[1]
    # busname_data = power_data[all_busName_key]
    all_busdata_key = filter(x->ismatch(r"all_busdata", x), collect(keys(power_data)))[1]
    raw_bus_data = power_data[all_busdata_key]
    nb_lines = size(raw_bus_data,1)
    if nb_lines!=length(bus_data["BaseKV"])
        error("sizes of all_busName and bus not equal")
    end
    bus = SortedDict{String, SortedDict{String,Any}}()
    # bus = SortedDict(bus_name(line) => SortedDict{String,Any}() for line in 1:nb_bus)
    for line in 1:nb_lines
        id_bus = Int(raw_bus_data[line,1])
        # index[id_bus] = line
        busname = bus_name(index[id_bus])
        baseKV = bus_data["BaseKV"][line]
        baseMVA = 1.0
        voltage_magnitude_min , voltage_magnitude_max  = get_bus_data(bus_data, ["VoltageMagnitudeMin", "VoltageMagnitudeMax"], line)
        if !haskey(bus, busname)
            bus[busname] = SortedDict{String,Any}()
        end
        bus[busname][volt_name()] = GOCVolt(busname, baseKV, baseMVA, voltage_magnitude_min, voltage_magnitude_max)
    end

    all_loadID_key = filter(x->ismatch(r"all_loadID", x), collect(keys(power_data)))[1]
    load_id = power_data[all_loadID_key]
    all_loaddata_key = filter(x->ismatch(r"all_loaddata", x), collect(keys(power_data)))[1]
    load_data = power_data[all_loaddata_key]

    for line in 1:size(load_data,1)
        id_bus = Int(load_data[line,1])
        busname = bus_name(index[id_bus])
        load = load_data[line,5] + im * load_data[line,6]
        id_load = id_elem(load_id,line)
        loadname = load_name(id_load)
        bus[busname][loadname] = GOCLoad(busname,loadname,load)
    end

    all_fixedshuntID_key = filter(x->ismatch(r"all_fixedshuntID", x), collect(keys(power_data)))[1]
    shunt_id = power_data[all_fixedshuntID_key]
    all_fixedshuntdata_key = filter(x->ismatch(r"all_fixedshuntdata", x), collect(keys(power_data)))[1]
    shunt_data = power_data[all_fixedshuntdata_key]
    for line in 1:size(shunt_data,1)
        id_bus = Int(shunt_data[line,1])
        busname = bus_name(index[id_bus])
        shunt = shunt_data[line,3] + im * shunt_data[line,4]
        id_shunt = id_elem(shunt_id,line)
        shuntname = shunt_name(id_shunt)
        bus[busname][shuntname] = GOCShunt(busname,shuntname,shunt)
    end

    return bus,index
end

function id_elem(array_id,line)
    if typeof(array_id) == String
        temp_string = array_id
    else
        temp_string = array_id[line]
        if typeof(temp_string)==Float64
            temp_string = Int(temp_string)
        end
    end
    return temp_string
    #return matchall(r"\S+" ,temp_string)[1]
end

#########################
function get_generator_data(generator_data, listofkeys, id)
    listofvalues = []
    for key in listofkeys
        if !haskey(generator_data, key)
            error("generator_data dictionary does not have a key $key")
        else
            append!(listofvalues, generator_data[key][id])
        end
    end
    return listofvalues
end


function generator_data_to_dict(generator_data,index) ##conversion array to dict
    generator_data_dict = SortedDict{String, SortedDict{String,SortedDict{Int64,Float64}}}()
    for line in 1:size(generator_data,1)
        bus_id = Int(generator_data[line,1])
        busname = bus_name(index[bus_id])
        gen_id = generator_data[line,2]
        gen_id = remove_simple_quotes_if_present(gen_id)
        if typeof(gen_id)==Float64
            gen_id = Int(gen_id)
        end
        genname = generator_name(gen_id)
        if !haskey(generator_data_dict, busname)
            generator_data_dict[busname] = SortedDict{String,SortedDict{Int64,Float64}}()
        end
        if !haskey(generator_data_dict[busname],genname)
            generator_data_dict[busname][genname] = SortedDict{Int64,Float64}()
        end
        term_id = Int(generator_data[line,3])
        value = generator_data[line,4]
        if term_id!=9
            generator_data_dict[busname][genname][term_id] = value
        end
    end
    return generator_data_dict
end

function remove_simple_quotes_if_present(id)
    if typeof(id)==SubString{String}
        if id[end-1]==' '
            return id[2:(end-2)]
        else
            return id[2:(end-1)]
        end
    else
        if id[end]==' '
            return id[1:(end-1)]
        else
            return id
        end
    end
end

"""
    add_generator_data!(power_data,generator_data_dict,bus,index)

Add generator data contained in `power_data` in DataSource field `bus` depending on `generator_data_dict` and `index`.\n
Return `bus`.

"""
function add_generator_data!(power_data,generator_data_dict,bus,index)
    nb_bus = Int(power_data["totalbus"])
    generator_data = power_data["generator"]
    for gen in 1:length(generator_data["BusNum"])
        bus_id = Int(get_generator_data(generator_data,["BusNum"],gen)[1])
        # all_busName_key = filter(x->ismatch(r"all_busName", x), collect(keys(power_data)))[1]
        # busname_data = power_data[all_busName_key]
        busname = bus_name(index[bus_id])
        gen_id, Pmin, Qmin, Pmax, Qmax, participation_factor =  get_generator_data(generator_data, ["UnitID","RealPowerMin","ReactivePowerMin","RealPowerMax","ReactivePowerMax", "ParticipationFactor"], gen)
        power_min = Pmin + im*Qmin
        power_max = Pmax + im*Qmax
        # cost_degrees = generator_data["RealPowerCostExponent"][gen,:]
        # cost_coeffs = generator_data["RealPowerCostCoefficient"][gen,:]
        gen_id = remove_simple_quotes_if_present(gen_id)
        if typeof(gen_id)==Float64
            gen_id = Int(gen_id)
        end
        generatorname = generator_name(gen_id)
        dict_obj_coeffs = SortedDict{Int64,Float64}()
        for (degree, value) in generator_data_dict[busname][generatorname]
            if (degree ∈ [0,1,2])
                dict_obj_coeffs[degree] = value
            end
        end
        if length(dict_obj_coeffs)<3
            #warn("A 2-order polynom cost with less than 3 values for generator $generatorname at bus $busname, default value 0")
            if !haskey(dict_obj_coeffs,0) dict_obj_coeffs[0] = 0 end
            if !haskey(dict_obj_coeffs,1) dict_obj_coeffs[1] = 0 end
            if !haskey(dict_obj_coeffs,2) dict_obj_coeffs[2] = 0 end
       end
        bus[busname][generatorname] = GOCGenerator(busname,generatorname,power_min,power_max,participation_factor,dict_obj_coeffs)
    end
    return bus
end

#########################
function get_branch_data(branch_data, listofkeys, id)
    listofvalues = []
    for key in listofkeys
        if !haskey(branch_data, key)
            error("branch_data dictionary does not have a key $key")
        else
            append!(listofvalues, branch_data[key][id])
        end
    end
    return listofvalues
end

"""
    read_branch_data(power_data, index)

Read branch data contained in `power_data` to fulfill DataSource field `link` depending on `index`. \n
Return dictionary `link`.
"""
function read_branch_data(power_data, index)
    nb_branch = Int(power_data["totalbranch"])
    nb_bus = Int(power_data["totalbus"])
    # all_busName_key = filter(x->ismatch(r"all_busName", x), collect(keys(power_data)))[1]
    # busname_data = power_data[all_busName_key]
    all_lineCKT_key = filter(x->ismatch(r"all_lineCKT", x), collect(keys(power_data)))[1]
    lines_data = power_data[all_lineCKT_key]
    nb_lines = length(lines_data)
    link = SortedDict{Link, SortedDict{String,Any}}()
    branch_data = power_data["branch"]

    #lines without transformer
    for branch in 1:nb_lines
        origin = Int(get_branch_data(branch_data, ["Origin"], branch)[1])
        destination = Int(get_branch_data(branch_data, ["Destination"], branch)[1])
        linkname = Link(bus_name(index[origin]),bus_name(index[destination]))
        branch_id = remove_simple_quotes_if_present(lines_data[branch])
        if typeof(branch_id)==Float64
            branch_id = Int(branch_id)
        end
        resistance, reactance, susceptance, power_magnitude_max = get_branch_data(branch_data, ["SeriesResistance","SeriesReactance","ChargingSusceptance","PowerMagnitudeMax"], branch)
        if !haskey(link, linkname)
            link[linkname] = SortedDict{String,Any}()
        end
        if resistance == reactance == 0
            println(linkname, " nullimpedance line without transformer")
            name_line = nullimpedance_notransformer_name(branch_id)
            link[linkname][name_line] = GOCNullImpedance_notransformer(linkname,name_line,susceptance, power_magnitude_max)
        else
            name_line = linepi_notransformer_name(branch_id)
            link[linkname][name_line] = GOCLineπ_notransformer(linkname,name_line,resistance,reactance,susceptance, power_magnitude_max)
        end
    end

    all_transformerCKT_key = filter(x->ismatch(r"all_transformerCKT", x), collect(keys(power_data)))[1]
    transformername_data = power_data[all_transformerCKT_key]

    #lines with transformers
    for branch in (nb_lines+1):nb_branch
        origin = Int(get_branch_data(branch_data, ["Origin"], branch)[1])
        destination = Int(get_branch_data(branch_data, ["Destination"], branch)[1])
        linkname = Link(bus_name(index[origin]),bus_name(index[destination]))
        branch_id = remove_simple_quotes_if_present(transformername_data[branch-nb_lines])
        if typeof(branch_id)==Float64
            branch_id = Int(branch_id)
        end
        resistance, reactance, susceptance, power_magnitude_max = get_branch_data(branch_data, ["SeriesResistance","SeriesReactance","ChargingSusceptance","PowerMagnitudeMax"], branch)
        transfo_ratio, transfo_phase = get_branch_data(branch_data, ["TapRatio", "PhaseShift"], branch)
        if !haskey(link, linkname)
            link[linkname] = SortedDict{String,Any}()
        end
        if resistance == reactance == 0
            println(linkname, " nullimpedance line with transformer")
            name_line = nullimpedance_withtransformer_name(branch_id)
            link[linkname][name_line] = GOCNullImpedance_withtransformer(linkname,name_line,susceptance, transfo_ratio,transfo_phase, power_magnitude_max)
        else
            name_line = linepi_withtransformer_name(branch_id)
            link[linkname][name_line] = GOCLineπ_withtransformer(linkname,name_line,resistance,reactance,susceptance, transfo_ratio,transfo_phase, power_magnitude_max)
        end
    end

    return link
end

#########################
"""
    scenarios_data(ds,gs,mp,power_data, contingency_data,index)

Use dictionaries `ds`, `gs` and `mp` to create `Scenario` structures for each contingency contained in `contingency_data`.\n
Return an `OPFproblems` structure : "scenario" => Scenario

"""
function scenarios_data(ds,gs,mp,contingency_data,index)
    output = SortedDict("BaseCase" => Scenario(ds,gs,mp))
    ds = output["BaseCase"].ds
    nb_contingencies = size(contingency_data,1)
    for ct in 1:nb_contingencies
        id_contingency = contingency_data[ct,1]
        scenario_name = scenarioname(id_contingency)
        type_contingency = contingency_data[ct,2]
        #B, T or G (Branch, Transformer, Generator)
        gs_scenario = copy(gs)
        gs_scenario.scenario = scenario_name
        ds_scenario_bus = SortedDict{String, SortedDict{String,Any}}()
        for busname in keys(ds.bus)
            ds_scenario_bus[busname] = SortedDict{String,Any}()
        end
        ds_scenario_link = SortedDict{Link, SortedDict{String,Any}}()
        for linkname in keys(ds.link)
            ds_scenario_link[linkname] = SortedDict{String,Any}()
        end
        ds_scenario = DataSource(ds_scenario_bus,ds_scenario_link)
        if type_contingency == "G"
            println(type_contingency)
            bus_id = contingency_data[ct,3]
            busname = bus_name(index[bus_id])
            gen_num = contingency_data[ct,4]
            CID = contingency_data[ct,5]
            if typeof(CID)==SubString{String}
                gen_ID = CID[2:(end-1)]
            else
                gen_ID = string(CID)
            end
            for (bus,dict) in ds.bus
                if bus==busname
                    for (buslabel,data) in dict
                        if buslabel != generator_name(gen_ID)
                            ds_scenario.bus[bus][buslabel] = ds.bus[bus][buslabel]
                        end
                    end
                else
                    for (buslabel,data) in dict
                        ds_scenario.bus[bus][buslabel] = ds.bus[bus][buslabel]
                    end
                end
            end
            for (linkname,dict) in ds_scenario.link
                for (linklabel,data) in ds.link[linkname]
                    ds_scenario.link[linkname][linklabel] = ds.link[linkname][linklabel]
                end
            end
        elseif type_contingency=="B" || type_contingency=="T"
            #println("Type of contingency $ct is $type_contingency")
            origin = Int(contingency_data[ct,3])
            destination = Int(contingency_data[ct,4])
            link_to_remove = Link(bus_name(index[origin]),bus_name(index[destination]))
            CID = contingency_data[ct,5]
            if typeof(CID)==SubString{String}
                id_line = CID[2:(end-1)]
            else
                id_line = string(CID)
            end
            for (busname,dict) in ds.bus
                for (buslabel,data) in dict
                    ds_scenario.bus[busname][buslabel] = ds.bus[busname][buslabel]
                end
            end
            for (linkname,dict) in ds_scenario.link
                if linkname.orig == link_to_remove.orig && linkname.dest == link_to_remove.dest
                    #println("link to remove : ", link_to_remove)
                    for (linklabel,data) in ds.link[linkname]
                        if !contains(linklabel,id_line)
                            ds_scenario.link[linkname][linklabel] = ds.link[linkname][linklabel]
                        end
                    end
                    if isempty(ds_scenario.link[linkname])
                        delete!(ds_scenario.link,linkname)
                    end
                else
                    for (linklabel,data) in ds.link[linkname]
                    ds_scenario.link[linkname][linklabel] = ds.link[linkname][linklabel]
                    end
                end
            end
        end
        output[scenario_name] = Scenario(ds_scenario,gs_scenario,mp)
    end


    return output
end
