

function parse_log(knitro_log_path)
    f = open(knitro_log_path, "r")
    lines = readlines(f)
    close(f)
    nb_it = Dict{Int64, Any}()
    feaserror = Dict{Int64, Any}()
    opterror = Dict{Int64, Any}()
    cpu = Dict{Int64, Any}()
    i = 0
    j = 0
    k = 0
    l = 0
    for line in lines
        if ismatch(r"Final feasibility error", line)
            j +=1
            feaserror[j] = split(line)[end]
        end
        if ismatch(r"Final optimality error", line)
            k +=1
            opterror[k] = split(line)[end]
        end
        if ismatch(r"# of iterations", line)
            i +=1
            nb_it[i] = split(line)[end]
        end
        if ismatch(r"CPU time", line)
            l +=1
            cpu[l] = matchall(r"\d+.\d+",line)[end]
        end

    end
    return nb_it, feaserror, opterror, cpu
end

# knitro_log_path = joinpath(pwd(), "..", "knitro_runs", "IEEE14_scenario_1", "Knitro_18_Apr_27_15_46_54.log")
#
# nb_it, feaserror, opterror = parse_log(knitro_log_path)
# println(nb_it)
# println(feaserror)
# println(opterror)
