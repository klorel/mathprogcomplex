include(joinpath(ROOT, "src_PowSysMod", "PowSysMod_body.jl"))

function buildPOP_1v1c()
    z = Variable("z", Complex)
    problemraw = Problem()
    add_variable!(problemraw, z)
    set_objective!(problemraw, -real(z))
    add_constraint!(problemraw, "ineq", abs2(z) << 4)
    return problemraw
end

function buildPOP_1v2c()
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
    z1 = Variable("z1", Complex)
    z2 = Variable("z2", Complex)
    problemraw = Problem()
    add_variable!(problemraw, z1); add_variable!(problemraw, z2);
    set_objective!(problemraw, 3-abs2(z1)-0.5*im*z1*conj(z2)^2+0.5im*z2^2*conj(z1))
    add_constraint!(problemraw, "eq1", (abs2(z1)-0.25*z1^2-0.25*conj(z1)^2) == 1)
    add_constraint!(problemraw, "eq2", (abs2(z1)+abs2(z2)) == 3)
    add_constraint!(problemraw, "eq3", (im*z2-im*conj(z2)) == 0)
    add_constraint!(problemraw, "ineq", (z2+conj(z2)) >> 0)
    return problemraw
end

function buildPOP_WB2()
    OPFpbs = load_OPFproblems(MatpowerInput, joinpath("data_Matpower", "matpower", "WB2.m"))
    return build_globalpb!(OPFpbs)
end

function buildPOP_WB2()
    z1 = Variable("z1", Complex)
    z2 = Variable("z2", Complex)
    problemraw = Problem()
    add_variable!(problemraw, z1); add_variable!(problemraw, z2);
    set_objective!(problemraw, 8*abs2(z2-z1))
    add_constraint!(problemraw, "VOLTM1", 0.9025 << abs2(z1) << 1.1025)
    add_constraint!(problemraw, "VOLTM2", 0.9025 << abs2(z2) << 1.1025)

    add_constraint!(problemraw, "BAL1", ((2+10im)*z1*conj(z2) + (2-10im)*z2*conj(z1) - 4*abs2(z2)) == 350)
    add_constraint!(problemraw, "BAL2", ((-10+2im)*z1*conj(z2) + (-10-2im)*z2*conj(z1) + 20*abs2(z2)) == -350)
    return problemraw
end