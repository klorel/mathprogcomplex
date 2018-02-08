using Poly
using Modeler

## Build problem objects
V1 = Variable("VOLT_1", Complex)
V2 = Variable("VOLT_2", Complex)

obj= 192 * V1*conj(V1) - (96.2+481im) * conj(V2)*V1 + (-96.2+481im) * conj(V1)*V2


Load2 = (-96.2 -481im) * conj(V1)*V2 + ( 96.2 +481im) * conj(V2)*V2
Load2_lb = -350 + 350im; Load2_ub = -350 + 350im
Unit1 = ( 96.2 +481im) * conj(V1)*V1 + (-96.2 -481im) * conj(V2)*V1
Unit1_lb = -400im; Unit1_ub = 600 + 400im
Voltm1 = abs2(V1)
Voltm1_lb = 0.90; Voltm1_ub = 1.10
Voltm2 = abs2(V2)
Voltm2_lb = 0.90; Voltm2_ub = 1.05

## Assemble problem
pb = Problem()
add_variable!(pb, V1)
add_variable!(pb, V2)

set_objective!(pb, obj)

add_constraint!(pb, "UNIT_1", Unit1_lb << Unit1 << Unit1_ub)
add_constraint!(pb, "LOAD_2", Load2_lb << Load2 << Load2_ub)
add_constraint!(pb, "VOLTM_1", Voltm1_lb << Voltm1 << Voltm1_ub)
add_constraint!(pb, "VOLTM_2", Voltm2_lb << Voltm2 << Voltm2_ub)

println(pb)

## Applying at a local solution
pt = Point([V1, V2], [-0.635304 + 0.706321im, 0.380961 + 0.898656im])

println("⟶ Constraints at Knitro local OPF solution:")
for (cstrName, cstr) in get_constraints(pb)
  println(cstrName, ": ", cstr.lb, " ≤ ", evaluate(cstr.p, pt), " ≤ ", cstr.ub)
end

println("⟶ Objective value: ", evaluate(pb.objective, pt))
println("⟶ Corresponding slacks: ")
print(get_slacks(pb, pt))
println("... and min slack: ", get_minslack(pb, pt))



## In real variables:
pb_R = pb_cplx2real(pb)
pt_R = Point([Variable("VOLT_1_r", Real),  Variable("VOLT_1_i", Real),  Variable("VOLT_2_r", Real),  Variable("VOLT_2_i", Real)], [-0.635304, 0.706321, 0.380961, 0.898656])

println("\n⟶ Constraints at Knitro local OPF solution, real problem:")
for (cstrName, cstr) in get_constraints(pb_R)
  println(cstrName, ": ", cstr.lb, " ≤ ", evaluate(cstr.p, pt_R), " ≤ ", cstr.ub)
end

println("⟶ Objective value, real problem: ", evaluate(pb_R.objective, pt_R))
