# Constraint
function <<(p::T, bound::Number) where T<: AbstractPolynomial
  return Constraint(p, -Inf-im*Inf, bound)
end

function <<(bound::Number, p::T) where T<: AbstractPolynomial
  return Constraint(p, bound, Inf+im*Inf)
end

function >>(p::T, bound::Number) where T<:AbstractPolynomial
  return bound << p
end

function >>(bound::Number, p::T) where T<:AbstractPolynomial
  return p << bound
end

function <<(bound::Number, cstr::Constraint)
  if real(bound) > real(cstr.ub) || imag(bound) > imag(cstr.ub)
    warn("<<(bnd, cstr): Creating a constraint with lower bound ", bound, " and upper bound ", cstr.ub)
  end
  return Constraint(cstr.p, bound, cstr.ub)
end
>>(bound::Number, cstr::Constraint) = cstr << bound

function <<(cstr::Constraint, bound::Number)
  if real(cstr.lb) > real(bound) || imag(cstr.lb) > imag(bound)
    warn("<<(cstr, bnd): Creating a constraint with lower bound ", cstr.lb, " and upper bound ", bound)
  end
  return Constraint(cstr.p, cstr.lb, bound)
end
>>(cstr::Constraint, bound::Number) = bound << cstr


function ==(p::T, x::Number) where T<:AbstractPolynomial
  return Constraint(p, x, x)
end

function Base.print(io::IO, cstr::Constraint)
  if cstr.lb == cstr.ub
    print(io, cstr.p, " = ", cstr.ub)
  else
    if cstr.lb != (-Inf-im*Inf)
      print(io, cstr.lb, " < ")
    end
    print(io, cstr.p)
    if cstr.ub !=(Inf+im*Inf)
      print(io, " < ", cstr.ub)
    end
  end
end

# Problem
Problem() = Problem(Polynomial(), Dict{String, Constraint}(), Dict{String, Type}())

function Base.print(io::IO, pb::Problem)
  print(io, "▶ variables: ")
  for (varName, typ) in sort(collect(pb.variables), by=x->x[1])
    print(io, Variable(varName, typ), " ")
  end
  println(io, "\n▶ objective: ", pb.objective)
  println(io, "▶ constraints: ")
  for (cstrName, cstr) in sort(collect(pb.constraints), by=x->x[1])
    @printf(io, "%10s: ", cstrName)
    println(io, cstr)
  end
end
