include(joinpath(ROOT, "src_PowSysMod", "PowSysMod_body.jl"))

function buildPOP_1v1c()
    z = Variable("z", Complex)
    problemraw = Problem()
    add_variable!(problemraw, z)
    set_objective!(problemraw, imag(z))
    add_constraint!(problemraw, "ineq", abs2(z) << 4)
    add_constraint!(problemraw, "ineq_rot", real(z*exp(-im*π/4)) >> 0)
    return problemraw
end

"""
    problemraw = buildPOP_2v3c

    Elliptic example problemp from Josz, Molzahn 2018 paper.
"""
function buildPOP_2v3c()
    z = Variable("z", Complex)
    problemraw = Problem()
    add_variable!(problemraw, z)
    set_objective!(problemraw, imag(z))
    add_constraint!(problemraw, "ineq", abs2(z) << 4)
    add_constraint!(problemraw, "ineq_rot", real(z*exp(-im*π/4)) >> 0)
    return problemraw
end

function build_WB2()
    OPFpbs = load_OPFproblems(MatpowerInput, joinpath("data_Matpower", "matpower", "WB2.m"))
    return build_globalpb!(OPFpbs)
end