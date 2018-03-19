using ArgParse
"""Test
Test to launch in complex-modeler/src_test folder for example D:\repo\complex-modeler\src_test>julia global_test.jl

# Arguments
- type of instance (GOCInput, MatpowerInput, IIDMInput)
if GOCInput : julia global_test.jl GOCInput <instances folder> <instance name> <scenario name> <epsilon>
if MatpowerInput or IIDMInput : julia global_test.jl MatpowerInput <instances folder> <instance name> <epsilon>


"""
#TODO:: option sur les print
function read_arguments(ARGS)
    if length(ARGS)==0
        error("Arguments must be :
        Argument 1 is the type of instance (GOCInput, MatpowerInput or IIDMInput)\n
         julia global_test.jl GOCInput <instances folder> <instance name> <scenario name> <epsilon>
         julia global_test.jl MatpowerInput <instances folder> <instance name> <epsilon>
         example : julia global_test.jl MatpowerInput ..\data_Matpower\matpower case9.m 1e-20
         julia global_test.jl IIDMInput <instances folder> <instance name> <epsilon>

        ")
    end
    typeofinput = ARGS[1]
    if typeofinput=="GOCInput"
        typeofinput = GOCInput
        if length(ARGS)!=5
            error("Arguments must be : julia global_test.jl GOCInput <instances folder> <instance name> <scenario name> <epsilon>")
        else
            instance_path = joinpath(pwd(), ARGS[2], ARGS[3],ARGS[4])
            eps = ARGS[5]
        end

    elseif typeofinput=="MatpowerInput"
        typeofinput = MatpowerInput
        if length(ARGS)!=4
            error("Arguments must be : julia global_test.jl MatpowerInput <instances folder> <instance name> <epsilon> ")
        else
            instance_path = joinpath(pwd(), ARGS[2], ARGS[3])
            eps = ARGS[4]
        end

    elseif typeofinput=="IIDMInput"
        typeofinput = IIDMInput
        if length(ARGS)!=4
            error("Arguments must be : julia global_test.jl IIDMInput <instances folder> <instance name> <epsilon> ")
        else
            instance_path = joinpath(pwd(), ARGS[2], ARGS[3])
            eps = ARGS[4]
        end

    else
        error()
    end
    return typeofinput, instance_path, eps
end

ROOT = pwd()
include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))
# include("para_fct.jl")
typeofinput, instance_path, eps = read_arguments(ARGS)


function build_and_solve_instance(typeofinput, instance_path)
    nb_variables = nb_constraints = 0
    obj = feas = 0.0

    tic()
    OPFpbs = load_OPFproblems(typeofinput, instance_path)

    ## Introducing coupling constraints on generator output
    (typeofinput != GOCInput) || introduce_Sgenvariables!(OPFpbs)

    ## Bulding optimization problem
    pb_global = build_globalpb!(OPFpbs)
    pb_global_real = pb_cplx2real(pb_global)

    ## Reading GOC initial point
    init_point = Point()
    (typeofinput != GOCInput) || (init_point = solution_point(instance_path))
    init_point_real = cplx2real(init_point)

    ## Exporting real problem
    if typeofinput != GOCInput
        instance_name  = String(split(instance_path,'\\')[end])
        amplexportpath = joinpath("..","knitro_runs", "$(instance_name[1:end-2])")
    else
        folder,scenario = split(instance_path, '\\')[end-1:end]
        folder = String(folder)
        scenario = String(scenario)
        amplexportpath = joinpath("..","knitro_runs", "$(folder[9:end])_$(scenario)")
    end
    export_to_dat(pb_global_real, amplexportpath, init_point_real)
    t_buildexport = toq()


    _, t_knitro, _ = @timed run_knitro(amplexportpath, joinpath(pwd(),"..","src_ampl"))
    pt_knitro, pt_GOC = read_Knitro_output(amplexportpath, pb_global_real)
    diff_GOC = norm(pt_GOC - pt_knitro, 2)
    feas,ctr = get_minslack(pb_global_real, pt_knitro)
    obj = get_objective(pb_global_real, pt_knitro)
    obj2 = get_objective(pb_global_real, pt_GOC)

    nb_variables = length(pb_global_real.variables)
    nb_constraints = length(pb_global_real.constraints)

    return String(split(instance_path,'\\')[end]) => (nb_variables, nb_constraints, obj, feas, t_buildexport, t_knitro, diff_GOC, obj2)
end

function main(ARGS)
    typeofinput, instance_path, eps = read_arguments(ARGS)

    (instance,(nb_variables, nb_constraints, obj, feas, t_buildexport, t_knitro,diff_GOC,obj2)) = build_and_solve_instance(typeofinput, instance_path)
    println("NB variables : ",nb_variables)
    println("NB_constraints : ",nb_constraints)
    println("Objective : ",obj)
    println("Feasibility : ",feas)
    println("t_buildexport : ", t_buildexport)
    println("t_knitro : ", t_knitro)
    println("Diff GOC point : ", diff_GOC)
    println("Objective GOC : ", obj2)

end

main(ARGS)
