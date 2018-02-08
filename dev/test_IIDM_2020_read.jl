ROOT=pwd()
include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))

path = joinpath("instances", "IIDM", "etat_reseau")
data_branches = readdlm(joinpath(path, "eod2020_network_branches.txt"))
data_buses = readdlm(joinpath(path, "eod2020_network_buses.txt"))
data_generators = readdlm(joinpath(path, "eod2020_network_generators.txt"))
data_limits = readdlm(joinpath(path, "eod2020_network_limits.txt"))
data_loads = readdlm(joinpath(path, "eod2020_network_loads.txt"))
data_substations = readdlm(joinpath(path, "eod2020_network_substations.txt"))

bus = Dict{String, Dict{String, Any}}()

for i=1:size(data_buses, 1)
    busid = Int(data_buses[i, 1])
    substationid = Int(data_buses[i, 2])
    V = data_buses[i, 4] # p.u.
    θ = data_buses[i, 5] # radians
    P, Q = data_buses[i, 6:7]
    basekV, minV, maxV = data_substations[substationid, 4:6]

    busname = "BUS_$busid"
    data_buses[i, 3] == 0 || warn("Buses - $busname : parameter cc!=0 ($(data_buses[i, 3])), not handled")
    data_buses[i, 8] == 0 || warn("Buses - $busname : parameter fault!=0 ($(data_buses[i, 8])), not handled")
    data_buses[i, 9] == 0 || warn("Buses - $busname : parameter curative!=0 ($(data_buses[i, 9])), not handled")

    bus[busname] = Dict{String, Any}()
    bus[busname]["Volt"] = IIDMVolt(busname, busid, string(substationid), basekV, minV, maxV, V, θ, P, Q)
end

for i=1:size(data_loads, 2)
    busid = data_loads[i, 2]
    P, Q = data_loads[i, 4:5]
    busname = "BUS_$busid"

    data_loads[i, 6] == 0 || warn("Loads - $busname : parameter fault!=0 ($(data_loads[i, 6])), not handled")
    data_loads[i, 7] == 0 || warn("Loads - $busname : parameter curative!=0 ($(data_loads[i, 7])), not handled")

    bus_elems = bus["$busname"]
    nb = count([ismatch(r"Load_", x) for x in keys(bus_elems)])
    bus_elems["Load_$(nb+1)"] = IIDMLoad(string(nb+1), P+im*Q)
end

for i=1:size(data_generators, 1)
    busid = data_generators[i, 2]
    minP, maxP = data_generators[i, 5:6]
    volt_regulator_on, targetV = data_generators[i, 11:12]
    minQmaxP, minQminP, maxQmaxP, maxQminP = data_generators[i, 7:10]
    S_bounds = Set([(minP+im*minQminP, minP+im*maxQminP), (maxP+im*minQmaxP, maxP+im*maxQmaxP)])
    busname = "BUS_$busid"

    data_generators[i, 2] == data_generators[i, 3] || warn("Generators - $busname : 'bus' is not 'con. bus', not handled")
    data_generators[i, 15] == 0 || warn("Loads - $busname : parameter fault!=0 ($(data_generators[i, 15])), not handled")
    data_generators[i, 16] == 0 || warn("Loads - $busname : parameter curative!=0 ($(data_generators[i, 16])), not handled")

    bus_elems = bus["$busname"]
    nb = count([ismatch(r"Gen_", x) for x in keys(bus_elems)])
    bus_elems["Gen_$(nb+1)"] = IIDMGeneratorFuel(string(nb+1), minP, maxP, targetV, volt_regulator_on, S_bounds)
end


link = Dict{Link, Dict{String, IIDMLine_π}}()
for i=1:size(data_branches, 1)
    busid1, busid2 = data_branches[i, 2:3]
    r, x, g1, g2, b1, b2 = data_branches[i, 7:12]
    p1, p2, q1, q2 = data_branches[i, 16:19]


    linkname = Link("BUS_$busid1", "BUS_$busid2")
    if !haskey(link, linkname)
        link[linkname] = Dict{String, IIDMLine_π}()
    end
    link_elems = link[link_elems]

    nb = count([ismatch(r"Branch_", x) for x in keys(bus_elems)])
    bus_elems["Branch_$(nb+1)"] = IIDMLine_π(r, x, g1+im*b1, g2+im*b2, p1+im*q1, p2+im*q2, basekV_orig, basekV_dest, lim)
end

find(data_branches[:, 13] .!= 1)
find(data_branches[:, 14] .!= -1)
find(data_branches[:, 15] .!= -1)


data = readdlm(filename)

find(data[:, 4] .== -1)


path = joinpath("instances", "IIDM", "etat_reseau")

data = readdlm(filename)

find(data[:, 1] .!= 0)
