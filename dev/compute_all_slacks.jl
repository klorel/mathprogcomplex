"""
    Check feasibility of Matpower solution
"""

include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))
include("read_QCQP.jl")

instances = [name[1:end-4] for name in readdir("instances")]

for instance_name in instances
  # Building complex and real problems
  filepath = abspath(joinpath("instances", instance_name*".dat"))
  pb = read_QCQP(filepath)
  pb_R = pb_cplx2real(pb)

  # Building complex and real solutions
  filepath = abspath(joinpath("solutions", instance_name*".sol"))
  pt_sol_C, pt_sol_R = build_QCQP_sol(filepath)

  # Compute and build minSlack and objective for both cplx and real problems
  minSlack, minCstrName = get_minslack(pb, pt_sol_C)
  val = get_objective(pb, pt_sol_C)
  @printf("%19s min slack: %.3e at %10s\t objective value : %.10e + im %.10e\n", instance_name, minSlack, minCstrName[1], real(val), imag(val))

  minSlack_R, minCstrName_R = get_minslack(pb_R, pt_sol_R)
  val_R = get_objective(pb_R, pt_sol_R)
  @printf("%19s min slack: %.3e at %10s\t objective value : %.10e + im %.10e\n", "C2R "*instance_name, minSlack_R, minCstrName_R[1], real(val_R), imag(val_R))
end
