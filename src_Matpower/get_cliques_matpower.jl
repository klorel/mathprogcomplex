include(joinpath(pwd(),"src_PowSysMod", "PowSysMod_body.jl"))


function get_cliques_matpower(instance_path)
    instances_path, instance_name = splitdir(instance_path)
    data_Matpower_path, ~ = splitdir(instances_path)
    blocks = readdlm(joinpath(data_Matpower_path, "blocks_AMD_clique", "$(instance_name[1:end-2])_sdp_blocks.txt"))
    cliques = SortedDict{String, SortedSet{Variable}}()
    println(blocks)
    println(size(blocks,1))
    for i in 1:size(blocks,1)
        block = blocks[i,1]
        var = blocks[i,2]
        VOLT, bus, re_or_im = split(var,"_")
        variable = Variable(variable_name("VOLT", bus_name(parse(bus)), "", basecase_scenario_name())*"_"*re_or_im,Real)
        if !haskey(cliques, block)
            cliques[block] = SortedSet{Variable}()
        end
        push!(cliques[block], variable)
    end
    return cliques
end


function get_cliques_matpower_forQCQP(instance_path)
    instances_path, instance_name = splitdir(instance_path)
    data_Matpower_path, ~ = splitdir(instances_path)
    blocks = readdlm(joinpath(data_Matpower_path, "blocks_AMD_clique", "$(instance_name[1:end-2])_sdp_blocks.txt"))
    cliques = SortedDict{String, SortedSet{Variable}}()
    println(blocks)
    println(size(blocks,1))
    for i in 1:size(blocks,1)
        block = blocks[i,1]
        var = blocks[i,2]
        variable = Variable(var,Real)
        if !haskey(cliques, block)
            cliques[block] = SortedSet{Variable}()
        end
        push!(cliques[block], variable)
    end
    return cliques
end

# instance = "case14.m"
# instance_path = joinpath(pwd(),"..","data", "data_Matpower", "matpower", instance)
#
# OPFpbs = load_OPFproblems(MatpowerInput, instance_path)
# ## Building optimization problem
# pb_global = build_globalpb!(OPFpbs)
# pb_global_real = pb_cplx2real(pb_global)
# cliques = get_cliques_matpower(instance_path)
# println(cliques)
