include(joinpath(ROOT, "src_PowSysMod", "PowSysMod_body.jl"))

function buildPOP_1v1c()
    z = Variable("z", Complex)
    problem = Problem()
    add_variable!(problem, z)
    set_objective!(problem, -real(z))
    add_constraint!(problem, "ineq", abs2(z) << 4)
    return problem
end

function buildPOP_1v2c()
    z = Variable("z", Complex)
    problem = Problem()
    add_variable!(problem, z)
    set_objective!(problem, imag(z))
    add_constraint!(problem, "ineq", abs2(z) << 4)
    add_constraint!(problem, "ineq_rot", real(z*exp(-im*π/4)) >> 0)
    return problem
end

"""
    problem = buildPOP_2v3c

    Elliptic example problemp from Josz, Molzahn 2018 paper.
"""
function buildPOP_2v3c()
    z1 = Variable("z1", Complex)
    z2 = Variable("z2", Complex)
    problem = Problem()
    add_variable!(problem, z1); add_variable!(problem, z2);
    set_objective!(problem, 3-abs2(z1)-0.5*im*z1*conj(z2)^2+0.5im*z2^2*conj(z1))
    add_constraint!(problem, "eq1", (abs2(z1)-0.25*z1^2-0.25*conj(z1)^2) == 1)
    add_constraint!(problem, "eq2", (abs2(z1)+abs2(z2)) == 3)
    add_constraint!(problem, "eq3", (im*z2-im*conj(z2)) == 0)
    add_constraint!(problem, "ineq", (z2+conj(z2)) >> 0)
    return problem
end

function buildPOP_WB2()
    OPFpbs = load_OPFproblems(MatpowerInput, joinpath("data_Matpower", "matpower", "WB2.m"))
    return build_globalpb!(OPFpbs)
end

function buildPOP_WB2_expl()
    z1 = Variable("z1", Complex)
    z2 = Variable("z2", Complex)
    problem = Problem()
    add_variable!(problem, z1); add_variable!(problem, z2);
    set_objective!(problem, 8*abs2(z2-z1))
    add_constraint!(problem, "VOLTM1", 0.9025 << abs2(z1) << 1.1025)
    add_constraint!(problem, "VOLTM2", 0.9025 << abs2(z2) << 1.1025)

    add_constraint!(problem, "BAL1", ((2+10im)*z1*conj(z2) + (2-10im)*z2*conj(z1) - 4*abs2(z2)) == 350)
    add_constraint!(problem, "BAL2", ((-10+2im)*z1*conj(z2) + (-10-2im)*z2*conj(z1) + 20*abs2(z2)) == -350)
    return problem
end


function buildPOPR_2v2cbis()
    x1 = Variable("x1", Real)
    x2 = Variable("x2", Real)
    problem = Problem()
    add_variable!(problem, x1); add_variable!(problem, x2)
    set_objective!(problem, -x2)
    add_constraint!(problem, "eq", (x1-0.5*x2) == 0)
    add_constraint!(problem, "ineq1", (x1 + 1) >> 0)
    add_constraint!(problem, "ineq2", (x1 + x2) << 0)
    add_constraint!(problem, "ineq_bnd", (x1^2 + x2^2) << 1)
    return problem
end

function buildPOPR_2v2c()
    x1 = Variable("x1", Real)
    x2 = Variable("x2", Real)
    problem = Problem()
    add_variable!(problem, x1); add_variable!(problem, x2)
    set_objective!(problem, -1*x1 - x2)
    add_constraint!(problem, "ineq1", -1 << x1 << 1)
    add_constraint!(problem, "ineq2", -1 << x2 << 1)
    add_constraint!(problem, "eq_up", (x1-2*x2) << 0)
    add_constraint!(problem, "eq_lo", (x1-2*x2) >> 0)
    # add_constraint!(problem, "eq", (x1-2*x2) == 0)
    return problem
end