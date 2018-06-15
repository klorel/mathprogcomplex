using KNITRO
ROOT = pwd()
include(joinpath("..", "src_PolynomialOptim", "PolynomialOptim.jl"))

# problem = Problem()
# x1 = Variable("x1", Real); add_variable!(problem, x1)
# x2 = Variable("x2", Real); add_variable!(problem, x2)
# x3 = Variable("x3", Real); add_variable!(problem, x3)
#
# set_objective!(problem, -2*x1+x2-x3)
# add_constraint!(problem, "ctr1", (x1*(4*x1-4*x2+4*x3-20)+x2*(2*x2-2*x3+9)+x3*(2*x3-13)+24) >> 0)
# add_constraint!(problem, "ctr2", (x1+x2+x3) << 4)
# add_constraint!(problem, "ctr3", (3*x2+x3) << 6)
# add_constraint!(problem, "def_x1", 0 << x1 << 2)
# add_constraint!(problem, "def_x2", 0 << x2)
# add_constraint!(problem, "def_x3", 0 << x3 << 3)
#
# print(problem)



# ###polynomial problem with complex variables
V1 = Variable("VOLT_1",Complex)
p_obj = 0.5*(V1+conj(V1))-0.5*im*(V1-conj(V1))
p_ctr1 = abs2(V1)
problem_poly=Problem()
add_variable!(problem_poly,V1)
add_constraint!(problem_poly, "ctr1", 1 << p_ctr1 << 1 )
set_objective!(problem_poly, p_obj)
print(problem_poly)
###########################################
pb_poly_real = pb_cplx2real(problem_poly)
println(pb_poly_real)

mysolver = KnitroSolver(KTR_PARAM_OUTLEV=4)
m, variables_jump, ~, ~ = get_JuMP_cartesian_model(pb_poly_real, mysolver)
print(m)
solve(m)
