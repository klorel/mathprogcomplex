"""Test
Test to launch in complex-modeler/src_test folder for example D:\repo\complex-modeler\src_test>julia test_feasibility_GOC.jl
# Arguments
- instances path
- instance folder
- scenario folder
- epsilon
Test to check feasibility of solution point given by GOC of VOLTM, Genbounds, Smax and coupling constraints
Print the names of the constraints violated along with a message

"""
function read_arguments(ARGS)
    if length(ARGS)!=4
        error("# Arguments
        - instances path
        - instance folder
        - scenario folder
        - epsilon")
    else
        instance_path = joinpath(pwd(),ARGS[1],ARGS[2],ARGS[3])
        epsilon = float(ARGS[4])
        return instance_path, epsilon
    end

end

instance_path, epsilon = read_arguments(ARGS)


ROOT = pwd()
include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))

violations = test_feasibility_GOC(instance_path,epsilon)
for (scenario, dict) in violations
    for (ctrname, msg) in dict
        @printf("%30s %50s\n", ctrname, msg)
    end
end
