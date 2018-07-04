"""
  pb = pb_cplx2real(pb_C::Problem)

Convert a complex Polynomial Optimization Problem into a real polynomial
problem in real variables.
"""
function pb_cplx2real(pb_C::Problem)
  pb = Problem()
  for (varName, varType) in get_variables(pb_C)
    if varType <: Complex
      varName_real, varName_imag = varname_cplx2real(varName)
      add_variable!(pb, Variable(varName_real, Real))
      add_variable!(pb, Variable(varName_imag, Real))
    else
      add_variable!(pb, Variable(varName, varType))
    end
  end
  (realPart, imagPart) = cplx2real(pb_C.objective)
  set_objective!(pb, realPart)
  for (cstrName, cstr) in get_constraints(pb_C)
    realPart, imagPart = cplx2real(cstr.p)
    cstrName_real, cstrName_imag = varname_cplx2real(cstrName)


    if length(realPart) != 0
      # Set bounds to proper infty, easier to detect... TODO : Missing attribute ctr_kind ?
      lb = real(cstr.lb)==-Inf ? -Inf-im*Inf : real(cstr.lb)
      ub = real(cstr.ub)== Inf ? +Inf+im*Inf : real(cstr.ub)
      cstrreal = lb << realPart << ub
      cstr.precond != :none && (cstrreal.precond = cstr.precond)
      add_constraint!(pb, cstrName_real, cstrreal)
    end
    if length(imagPart) != 0
      lb = imag(cstr.lb)==-Inf ? -Inf-im*Inf : imag(cstr.lb)
      ub = imag(cstr.ub)== Inf ? +Inf+im*Inf : imag(cstr.ub)
      cstrimag = lb << imagPart << ub
      cstr.precond != :none && (cstrimag.precond = cstr.precond)
      add_constraint!(pb, cstrName_imag, cstrimag)
    end
  end
  return pb
end


function pb_cplx2real_add(pb_C::Problem)
  # println("---- pb_cplx2real_add ----")
  pb = Problem()
  # println("variables")
  # tic()
  for (varName, varType) in get_variables(pb_C)
    if varType <: Complex
      varName_real, varName_imag = varname_cplx2real(varName)
      add_variable!(pb, Variable(varName_real, Real))
      add_variable!(pb, Variable(varName_imag, Real))
    else
      add_variable!(pb, Variable(varName, varType))
    end
  end
  # toc()
  # println("objective")
  # tic()
  (realPart, imagPart) = cplx2real_add(pb_C.objective)
  set_objective!(pb, realPart)
  # toc()
  # println("constraints")
  # println("nb = ", length(get_constraints(pb_C)))
  # tic()
  constraints = get_constraints(pb_C)
  for (cstrName, cstr) in constraints
    realPart, imagPart = cplx2real_add(cstr.p)
    cstrName_real, cstrName_imag = varname_cplx2real(cstrName)


    if length(realPart) != 0
      # Set bounds to proper infty, easier to detect... TODO : Missing attribute ctr_kind ?
      lb = real(cstr.lb)==-Inf ? -Inf-im*Inf : real(cstr.lb)
      ub = real(cstr.ub)== Inf ? +Inf+im*Inf : real(cstr.ub)
      cstrreal = lb << realPart << ub
      cstr.precond != :none && (cstrreal.precond = cstr.precond)
      add_constraint!(pb, cstrName_real, cstrreal)
    end
    if length(imagPart) != 0
      lb = imag(cstr.lb)==-Inf ? -Inf-im*Inf : imag(cstr.lb)
      ub = imag(cstr.ub)== Inf ? +Inf+im*Inf : imag(cstr.ub)
      cstrimag = lb << imagPart << ub
      cstr.precond != :none && (cstrimag.precond = cstr.precond)
      add_constraint!(pb, cstrName_imag, cstrimag)
    end
  end
  # toc()
  # println("--------------------")
  return pb
end


function infer_problem!(pb::Problem, pt::Point)
  pb.objective = evaluate(pb.objective, pt, partial=true)

  println(pt)

  for var in keys(pt)
    if haskey(pb.variables, var.name)
      delete!(pb.variables, var.name)
    end
  end

  for (ctrname, ctr) in pb.constraints
    pb.constraints[ctrname].p = evaluate(ctr.p, pt, partial=true)
    if isa(ctr.p, Number)
      ctr.p ≤ ctr.ub || warn("infer_problem!(): Constraint $ctrname has body $(ctr.p) with bounds ($(ctr.lb), $(ctr.ub).\nProblem may be infeasable.")
      ctr.lb ≤ ctr.p || warn("infer_problem!(): Constraint $ctrname has body $(ctr.p) with bounds ($(ctr.lb), $(ctr.ub).\nProblem may be infeasable.")
    end
  end
end
