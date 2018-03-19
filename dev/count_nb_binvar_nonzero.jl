###count the number of binary variables at 1

function treat_sol(csv_file)
    tab = readdlm(csv_file, ';')
    nb_lines = size(tab,1)
    vars = Dict{String, Float64}()
    for i in 1:nb_lines
        varname = tab[i,1]
        if ismatch(r"BinVolt", varname)
            vars[varname] = tab[i,2]
        end
    end
    return vars
end

function return_nonzerovars(vars,epsilon)
    nonzerosvars = Dict{String, Float64}()
    for (varname, value) in vars
        if abs(value) > epsilon
            nonzerosvars[varname] = value
        end
    end
    return nonzerosvars
end

function count_nonzerobinvars(csv_file, epsilon)
    vars = treat_sol(csv_file)
    nonzerosvars = return_nonzerovars(vars,epsilon)
    return length(nonzerosvars)
end

pwd()
instances_folder = joinpath(pwd(),"..","data_GOC")
folders = ["Phase_0_IEEE14","Phase_0_Modified_IEEE14","Phase_0_OriginalDataset_IEEE14", "Phase_0_RTS96","Phase_0_Modified_RTS96","Phase_0_OriginalDataset_RTS96"]
nb_scenarios = [100,100,75,100,100,30]
epsilon = 1e-6
dict = Dict{Tuple{String,String},Int64}()

for i in 1:6
    folder = folders[i]
    nb_scen = nb_scenarios[i]
    for j in 1:nb_scen
        scenario = "scenario_$j"
        instance_path = joinpath(pwd(), "..", instances_folder,folder, scenario)
        instance = "$(folder[9:end])_$(scenario)"
        csv_file = joinpath(pwd(), "..","knitro_runs",instance, "solution_GOC_point.csv")
        nb = count_nonzerobinvars(csv_file, epsilon)
        if nb!=0
            dict[(folder, scenario)] = nb
        end
    end
end


folder = "Phase_0_IEEE14"
scenario = "scenario_1"
instance_path = joinpath(pwd(), "..", instances_folder,folder, scenario)
instance = "$(folder[9:end])_$(scenario)"
csv_file = joinpath(pwd(), "..","knitro_runs",instance, "solution_GOC_point.csv")
nb = count_nonzerobinvars(csv_file, epsilon)
