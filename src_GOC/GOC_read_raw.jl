#
# include(joinpath(pwd(),"src_PowSysMod", "PowSysMod_body.jl"))
#
# raw = "powersystem.raw"
# gen = "generator.csv"
# con = "contingency.csv"
#
# instance_path = joinpath(pwd(), "data_GOC", "Phase_0_IEEE14_1Scenario","scenario_1")
# instance_path = joinpath(pwd(), "data_GOC", "Phase_0_RTS96","scenario_1")
# instance_path = joinpath(pwd(), "data_GOC", "Phase_0_Feas179","scenario_1")
#
# read(matopen(joinpath(instance_path, "pscopf_data.mat")))
#
# rawfile = joinpath(instance_path,raw)
# genfile = joinpath(instance_path, gen)
# confile = joinpath(instance_path, con)
#
# baseMVA, bus_data, load_data, shunt_data, generator_data, branch_data, transformer_data = read_rawfile(rawfile)
# gen_data_csv = read_genfile(genfile)
# con_data_csv = read_confile(confile)
#
# read_GOCfiles(rawfile, genfile,confile)

function read_genfile(genfile)
    gen_data = readdlm(genfile,',',skipstart=1)
    return gen_data
end

function read_confile(confile)
    con_data = readdlm(confile,',',skipstart=1)
    return con_data
end

function read_rawfile(rawfile)
    raw_data = readdlm(rawfile, ',')
    baseMVA = raw_data[1,2]
    i_begin = 2
    while raw_data[i_begin,4] == ""
        i_begin +=1
    end
    ##BUS DATA
    i_end = get_data_iend(i_begin,raw_data, "0 / END OF BUS DATA")
    bus_data = raw_data[i_begin:i_end,:]
    ##LOAD DATA
    i_begin = i_end + 2
    i_end = get_data_iend(i_begin,raw_data, "0 / END OF LOAD DATA")
    load_data = raw_data[i_begin:i_end,:]
    ## FIXED SHUNT DATA
    i_begin = i_end + 2
    i_end = get_data_iend(i_begin,raw_data, "0 / END OF FIXED SHUNT DATA")
    shunt_data = raw_data[i_begin:i_end,:]
    ## GENERATOR DATA
    i_begin = i_end + 2
    i_end = get_data_iend(i_begin,raw_data, "0 / END OF GENERATOR DATA")
    generator_data = raw_data[i_begin:i_end,:]
    ## BRANCH DATA
    i_begin = i_end + 2
    i_end = get_data_iend(i_begin,raw_data, "0 / END OF BRANCH DATA")
    branch_data = raw_data[i_begin:i_end,:]
    ## TRANSFORMER DATA
    i_begin = i_end + 2
    i_end = get_data_iend(i_begin,raw_data, "0 / END OF TRANSFORMER DATA")
    transformer_data = raw_data[i_begin:i_end,:]
    return baseMVA, bus_data, load_data, shunt_data, generator_data, branch_data, transformer_data
end



function get_data_iend(i_begin, data, end_sentence)
    i = i_begin
    imax = size(data,1)
    while data[i,1] != end_sentence && i<imax
        i+=1
    end
    return i-1
end


###read bus data
function read_data_bus_fromraw(bus_data, load_data, shunt_data, baseMVA)
    nb_bus = size(bus_data,1)
    # index = SortedDict(bus_data[i,1] => i for i in 1:nb_bus)
    index = SortedDict(bus_data[i,1] => bus_data[i,1] for i in 1:nb_bus)
    bus = SortedDict{String, SortedDict{String,Any}}()
    for i in 1:nb_bus
        bus[bus_name(index[bus_data[i,1]])] = SortedDict{String,Any}()
    end
    for i in 1:nb_bus
        id_bus = bus_data[i,1]
        busname = bus_name(index[id_bus])
        baseKV = bus_data[i,3]
        baseMVA = baseMVA
        voltage_magnitude_min = bus_data[i,11]
        voltage_magnitude_max  = bus_data[i,10]
        bus[busname][volt_name()] = GOCVolt(id_bus, baseKV, baseMVA, voltage_magnitude_min, voltage_magnitude_max)
    end
    ### LOAD DATA
    nb_load = size(load_data,1)
    for i in 1:nb_load
        id_bus = load_data[i,1]
        busname = bus_name(index[id_bus])
        load_re = load_data[i,6]
        load_im = load_data[i,7]
        if load_re!=0 || load_im!=0
            load = load_re  + im * load_im
            id_load = load_data[i,2]
            id_load = remove_simple_quotes_and_spaces_if_present(id_load)
            loadname = load_name(id_load)
            bus[busname][loadname] = GOCLoad(id_bus,loadname,load)
        end
    end
    ### SHUNT DATA
    nb_shunt = size(shunt_data,1)
    for i in 1:nb_shunt
        id_bus = shunt_data[i,1]
        busname = bus_name(index[id_bus])
        shunt_re = shunt_data[i,4]
        shunt_im = shunt_data[i,5]
        if shunt_re!=0 || shunt_im!=0
            shunt = shunt_re + im * shunt_im
            id_shunt = shunt_data[i,2]
            id_shunt = remove_simple_quotes_and_spaces_if_present(id_shunt)
            shuntname = shunt_name(id_shunt)
            bus[busname][shuntname] = GOCShunt(id_bus,shuntname,shunt)
        end
    end
    return bus,index
end

function remove_simple_quotes_and_spaces_if_present(id)
    if typeof(id)==SubString{String}
        if id[1]==''' && id[end]=='''
            id_without_spaces = split(id[2:end-1])[1]
        else
            id_without_spaces = split(id)[1]
    end
        return id_without_spaces
    else
        return id
end

end

function convert_gen_data_csv_into_dict(gen_data_csv, index)
    nb_lines = size(gen_data_csv, 1)
    gen_data_csv_dict = SortedDict{String, SortedDict{String,SortedDict{Int64,Float64}}}()
    for i in 1:nb_lines
        id_bus = Int64(gen_data_csv[i,1])
        busname = bus_name(index[id_bus])
        gen_id = gen_data_csv[i,2]
        gen_id = remove_simple_quotes_and_spaces_if_present(gen_id)
        if typeof(gen_id)==Float64
            gen_id = Int64(gen_id)
        end
        generatorname = generator_name(gen_id)
        term_id = Int64(gen_data_csv[i,3])
        coeff = gen_data_csv[i,4]
        if !haskey(gen_data_csv_dict, busname)
            gen_data_csv_dict[busname] = SortedDict{String,SortedDict{Int64,Float64}}()
        end
        if !haskey(gen_data_csv_dict[busname], generatorname)
            gen_data_csv_dict[busname][generatorname] = SortedDict{Int64,Float64}()
        end
        gen_data_csv_dict[busname][generatorname][term_id] = coeff
    end
    return gen_data_csv_dict
end



function add_generator_data_fromraw!(generator_data, gen_data_csv_dict, bus, index)
    nb_gen = size(generator_data,1)
    for i in 1:nb_gen
        id_bus = generator_data[i,1]
        busname = bus_name(index[id_bus])
        gen_id = generator_data[i,2]
        gen_id2 = remove_simple_quotes_and_spaces_if_present(gen_id)
        if typeof(gen_id2)==Float64
            gen_id2 = Int64(gen_id)
        end
        generatorname = generator_name(gen_id2)
        Pmin = generator_data[i,18]
        Qmin = generator_data[i,6]
        Pmax = generator_data[i,17]
        Qmax = generator_data[i,5]
        participation_factor = gen_data_csv_dict[busname][generatorname][9]
        power_min = Pmin + im*Qmin
        power_max = Pmax + im*Qmax
        dict_obj_coeffs = SortedDict{Int64,Float64}()
        for (degree, value) in gen_data_csv_dict[busname][generatorname]
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
        bus[busname][generatorname] = GOCGenerator(id_bus,gen_id,power_min,power_max,participation_factor,dict_obj_coeffs)
    end
    return bus
end

function modify_transformer_data(transformer_data)
    nb_lines, nb_col = size(transformer_data)
    nb_transfo = nb_lines/4
    if !isinteger(nb_transfo)
        exit("Transformer data is not on the format : 4 lines for 1 transformer")
    end
    nb_transfo = Int64(nb_transfo)
    transfo_data = Array{Any}(nb_transfo,nb_col*4)

    for i in 1:nb_transfo
        line = (i-1)*4+1
        transfo_data[i,1:nb_col] = transformer_data[line,:]
        transfo_data[i,(nb_col+1):(2*nb_col)] = transformer_data[(line+1),:]
        transfo_data[i,(2*nb_col+1):(3*nb_col)] = transformer_data[(line+2),:]
        transfo_data[i,(3*nb_col+1):(4*nb_col)] = transformer_data[(line+3),:]
    end
    return transfo_data
end



function read_data_branch_fromraw(branch_data, transfo_data, index)
    nb_branch = size(branch_data,1)
    link = SortedDict{Link, SortedDict{String,Any}}()
    for i in 1:nb_branch
        orig = Int64(branch_data[i,1])
        dest = Int64(branch_data[i,2])
        linkname = Link(bus_name(index[orig]),bus_name(index[dest]))
        br_id = branch_data[i,3]
        branch_id = remove_simple_quotes_and_spaces_if_present(branch_data[i,3])
        if typeof(branch_id)==Float64
            branch_id = Int64(branch_id)
        end
        resistance = branch_data[i,4]
        reactance = branch_data[i,5]
        susceptance = branch_data[i,6]
        power_magnitude_max = branch_data[i,7]

        if !haskey(link, linkname)
            link[linkname] = SortedDict{String,Any}()
        end
        if resistance == reactance == 0
            println(linkname, " nullimpedance line without transformer")
            name_line = nullimpedance_notransformer_name(branch_id)
            link[linkname][name_line] = GOCNullImpedance_notransformer(orig,dest,br_id,susceptance, power_magnitude_max)
        else
            name_line = linepi_notransformer_name(branch_id)
            link[linkname][name_line] = GOCLineπ_notransformer(orig,dest,br_id,resistance,reactance,susceptance, power_magnitude_max)
        end
    end

    nb_transfo = size(transfo_data,1)
    nb_col = size(branch_data,2)
    for i in 1:nb_transfo
        orig = Int64(transfo_data[i,1])
        dest = Int64(transfo_data[i,2])
        linkname = Link(bus_name(index[orig]),bus_name(index[dest]))
        br_id = transfo_data[i,4]
        branch_id = remove_simple_quotes_and_spaces_if_present(transfo_data[i,4])
        if typeof(branch_id)==Float64
            branch_id = Int(branch_id)
        end
        resistance = transfo_data[i,(nb_col+1)]
        reactance = transfo_data[i,(nb_col+2)]
        susceptance = 0.0
        power_magnitude_max = transfo_data[i,(2*nb_col+4)]
        transfo_ratio = transfo_data[i,(2*nb_col+1)]
        transfo_phase = transfo_data[i,(2*nb_col+3)]
        if !haskey(link, linkname)
            link[linkname] = SortedDict{String,Any}()
        end
        if resistance == reactance == 0
            println(linkname, " nullimpedance line with transformer")
            name_line = nullimpedance_withtransformer_name(branch_id)
            link[linkname][name_line] = GOCNullImpedance_withtransformer(orig,dest,br_id,susceptance, transfo_ratio,transfo_phase, power_magnitude_max)
        else
            name_line = linepi_withtransformer_name(branch_id)
            link[linkname][name_line] = GOCLineπ_withtransformer(orig,dest,br_id,resistance,reactance,susceptance, transfo_ratio,transfo_phase, power_magnitude_max)
        end
    end

    return link
end


function read_GOCfiles(rawfile, genfile,confile)
    baseMVA, bus_data, load_data, shunt_data, generator_data, branch_data, transformer_data = read_rawfile(rawfile)
    gen_data_csv = read_genfile(genfile)
    con_data_csv = read_confile(confile)
    bus,index = read_data_bus_fromraw(bus_data, load_data, shunt_data, baseMVA)
    gen_data_csv_dict = convert_gen_data_csv_into_dict(gen_data_csv, index)
    bus = add_generator_data_fromraw!(generator_data, gen_data_csv_dict, bus, index)
    transfo_data = modify_transformer_data(transformer_data)
    link = read_data_branch_fromraw(branch_data, transfo_data, index)

    ds = DataSource(bus,link)
    node_linksin, node_linksout = SortedDict{String, SortedSet{Link}}(), SortedDict{String, SortedSet{Link}}()
    node_vars = SortedDict{String, SortedDict{String, Variable}}()
    link_vars = SortedDict{Link, SortedDict{String, Variable}}()
    gs = GridStructure("BaseCase", node_linksin, node_linksout)
    node_formulations = SortedDict{String, SortedDict{String, Symbol}}()
    link_formulations = SortedDict{Link, SortedDict{String, Symbol}}()
    mp = MathematicalProgramming(node_formulations, link_formulations, node_vars,link_vars)
    ##read scenarios
    OPFproblems = scenarios_data(ds, gs, mp, con_data_csv,index)
    return OPFproblems
end
