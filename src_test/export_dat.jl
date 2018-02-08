ROOT = pwd()
include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))

"""Test
Test to launch in complex-modeler/src_test folder for example D:\repo\complex-modeler\src_test>julia global_test.jl

# Arguments
- type of instance (GOCInput, MatpowerInput, IIDMInput)
if GOCInput : julia global_test.jl GOCInput <instances folder> <instance name> <scenario name>
if MatpowerInput or IIDMInput : julia global_test.jl MatpowerInput <instances folder> <instance name>


"""
function read_arguments(ARGS)
    if length(ARGS)==0
        error("Arguments must be :
        Argument 1 is the type of instance (GOCInput, MatpowerInput or IIDMInput)\n
         julia global_test.jl GOCInput <instances folder> <instance name> <scenario name>
         julia global_test.jl MatpowerInput <instances folder> <instance name>
         example : julia global_test.jl MatpowerInput ..\data_Matpower\matpower case9.m
         julia global_test.jl IIDMInput <instances folder> <instance name>

        ")
    end
    inputtype = ARGS[1]
    if inputtype=="GOCInput"
        inputtype = GOCInput
        if length(ARGS)!=5
            error("Arguments must be : julia global_test.jl GOCInput <instances folder> <instance name> <scenario name> <otuput folder>")
        else
            instance_path = joinpath(pwd(), ARGS[2], ARGS[3],ARGS[4])
            output_path = ARGS[3]
        end

    elseif inputtype=="MatpowerInput"
        inputtype = MatpowerInput
        if length(ARGS)!=4
            error("Arguments must be : julia global_test.jl MatpowerInput <instances folder> <instance name> <otuput folder>")
        else
            instance_path = joinpath(pwd(), ARGS[2], ARGS[3])
            output_path = ARGS[4]
        end

    elseif inputtype=="IIDMInput"
        inputtype = IIDMInput
        if length(ARGS)!=4
            error("Arguments must be : julia global_test.jl IIDMInput <instances folder> <instance name> <otuput folder>")
        else
            instance_path = joinpath(pwd(), ARGS[2], ARGS[3])
            output_path = ARGS[4]
        end
    else
        error("Unexpected input type.")
    end
    return inputtype, instance_path, output_path
end

function main(ARGS)
    inputtype, instance_path, output_path = read_arguments(ARGS)

    ## Building network description
    OPFpbs = load_OPFproblems(inputtype, instance_path)

    ## Setting Generators labels for coupling constraints if GOC instance
    inputtype == GOCInput && introduce_Sgenvariables!(OPFpbs)

    ## Building polynomila optim problem
    pb_global = build_globalpb!(OPFpbs)

    ## Exporting
    export_to_dat(pb_global, output_path)
end

main(ARGS)
