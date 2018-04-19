### test feasibility of a solution
function read_arguments(ARGS)
    if length(ARGS)==0 || length(ARGS) != 4
        error("Arguments must be :
         julia test_feasiblity_GOC_solution.jl <instances folder> <instance id> <scenario id> <solution path>\n
         instance id = 14_1 for Phase_0_IEEE14_1Scenario
                     = 14 for Phase_0_IEEE14
                     = 14m for Phase_0_Modified_IEEE14
                     = 96 for Phase_0_RTS96
                     = 96m for Phase_0_Modified_RTS96
                     = 179f for Phase_0_Feas179
                     = 179nf for Phase_0_Infeas179
        ")
    else
        instances_folder = ARGS[1]
        id_instance = ARGS[2]
        if id_instance == "14_1"
            instance_name = "Phase_0_IEEE14_1Scenario"
        elseif id_instance == "14"
            instance_name = "Phase_0_IEEE14"
        elseif id_instance == "14m"
            instance_name = "Phase_0_Modified_IEEE14"
        elseif id_instance == "96"
            instance_name = "Phase_0_RTS96"
        elseif id_instance == "96m"
            instance_name = "Phase_0_Modified_RTS96"
        elseif id_instance == "179f"
            instance_name = "Phase_0_Feas179"
        elseif id_instance == "179nf"
            instance_name = "Phase_0_Infeas179"
        else
            error("id_instance must be one of the following ones :
            instance id = 14_1 for Phase_0_IEEE14_1Scenario
                        = 14 for Phase_0_IEEE14
                        = 14m for Phase_0_Modified_IEEE14
                        = 96 for Phase_0_RTS96
                        = 96m for Phase_0_Modified_RTS96
                        = 179f for Phase_0_Feas179
                        = 179nf for Phase_0_Infeas179   ")
        end
        scenario_id = ARGS[3]
        instance_path = joinpath(pwd(), instances_folder, instance_name, "scenario_$scenario_id")
        solution_path = ARGS[4]
    end

    return instance_path, solution_path
end


instance_path, solution_path = read_arguments(ARGS)

using MAT
include(joinpath(pwd(),"..","src_PowSysMod", "PowSysMod_body.jl"))


function test_feasiblity_GOC_solution(instance_path, solution_path)
    raw = "powersystem.raw"
    gen = "generator.csv"
    con = "contingency.csv"
    rawFile = joinpath(instance_path,raw)
    genFile = joinpath(instance_path, gen)
    contFile = joinpath(instance_path, con)

    OPFpbs = load_OPFproblems(rawFile, genFile, contFile)
    introduce_Sgenvariables!(OPFpbs)
    ## Building optimization problem
    pb_global = build_globalpb!(OPFpbs)

    ## convert into real problem
    pb_global_real = pb_cplx2real(pb_global)

    ## read solution
    sol_txt = read_solution_point_GOC(instance_path, solution_path)
    pt_txt = cplx2real(sol_txt)

    return get_minslack(pb_global_real, pt_txt)
end

println("instance_path: ", instance_path)
println("solution_path: ", solution_path)
println("get_minslack point from txt files:", test_feasiblity_GOC_solution(instance_path, solution_path))
