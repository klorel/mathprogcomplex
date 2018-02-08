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

  realPart, imagPart = cplx2real(pb_C.objective)
  set_objective!(pb, realPart)

  for (cstrName, cstr) in get_constraints(pb_C)
    realPart, imagPart = cplx2real(cstr.p)
    cstrName_real, cstrName_imag = varname_cplx2real(cstrName)
    if length(realPart) != 0
      cstrreal = real(cstr.lb) << realPart << real(cstr.ub)
      cstr.precond != :none && (cstrreal.precond = cstr.precond)
      add_constraint!(pb, cstrName_real, cstrreal)
    end
    if length(imagPart) != 0
      cstrimag = imag(cstr.lb) << imagPart << imag(cstr.ub)
      cstr.precond != :none && (cstrimag.precond = cstr.precond)
      add_constraint!(pb, cstrName_imag, cstrimag)
    end
  end
  return pb
end
