## objective
get_objective(pb::Problem) = pb.objective

"""
get_objective(pb::Problem, pt::Point)

Return the polynomial objective evaluated at `pt`.
"""
get_objective(pb::Problem, pt::Point) = evaluate(pb.objective, pt)

set_objective!(pb::Problem, p::Polynomial) = pb.objective = p


## variables
get_variables(pb::Problem) = pb.variables

function get_variabletype(pb::Problem, varName::String)
  if !haskey(pb.variables, varName)
    error("get_variable(): Current Problem has no variable named ", varName)
  end
  pb.variables[varName]
end

has_variable(pb::Problem, var::Variable) = haskey(pb.variables, var.name) && pb.variables[var.name] == var.kind

function add_variable!(pb::Problem, x::Variable)
  if haskey(pb.variables, x.name) && !(x.kind <: pb.variables[x.name])
    error("add_variable!(): Attempting to add variable ", x.name, " (", x.kind, ") when ", x, " (", get_variabletype(pb, x.name), ") already exists.")
  end
  pb.variables[x.name] = x.kind
end

function add_variable!(pb::Problem, x::Pair{String, T}) where T
  return add_variable!(pb, Variable(x[1], x[2]))
end

## constraints
has_constraint(pb::Problem, cstrName::String) = haskey(pb.constraints, cstrName)

get_constraints(pb::Problem) = pb.constraints

"""
  cstr = get_constraint(pb::Problem, cstrname::String)

Return the `cstrname` constraint from `pb`.
"""
get_constraint(pb::Problem, cstrname::String) = pb.constraints[cstrname]

"""
  cstr = get_constraint(pb::Problem, pt::Point)

Return the point (dict) of constraints names to their body's value when
evaluated at `pt`.
"""
function get_constraints(pb::Problem, pt::Point)
  img = Point()
  for (cstrName, cstr) in pb.constraints
    img[Variable(cstrName, Complex)] = evaluate(cstr.p, pt)
  end
  return img
end

"""
  add_constraint!(pb::Problem, cstrName::String, cstr::Constraint)

Add the constraint `cstr` under the `cstrname` name in the `pb` problem.
"""
function add_constraint!(pb::Problem, cstrName::String, cstr::Constraint)
  if haskey(pb.constraints, cstrName)
    warn("add_constraint!(): A constraint with that name already exists ($cstrName)")
  end
  if real(cstr.ub) < real(cstr.lb) || imag(cstr.ub) < imag(cstr.lb)
    warn("add_constraint!(): ", cstrName, " Lower bound is higher than upper bound ($(cstr.lb) - $(cstr.ub))")
  end
  pb.constraints[cstrName] = cstr
end

rm_constraint!(pb::Problem, cstrName::String) = pop!(pb.constraints, cstrName)

"""
  pt = get_slacks(pb::Problem, pt::Point)

Return a point associating a constraint name to its slack, i.e. the minimum
algebraic distance between the real and imaginary parts of the body's value at
`pt` and its bounds.
"""
function get_slacks(pb::Problem, pt::Point)
  var_arr = Variable[]
  val_arr = Complex[]
  for (cstrName, cstr) in pb.constraints
    val = evaluate(cstr.p, pt)
    isa(val, Number) || error("get_slacks(): constraint $cstrName not fully evaluated at provided point.\nEvaluated value is $val.")
    push!(var_arr, Variable(cstrName, Complex))
    push!(val_arr, min(real(val-cstr.lb), real(cstr.ub-val)) + min(imag(val-cstr.lb), imag(cstr.ub-val))*im)
  end
  return Point(var_arr, val_arr)
end

function get_minslack(pb::Problem, pt::Point)
  minSlack = +Inf
  minCstrName = ""
  slacks = get_slacks(pb, pt)
  for (cstrName, slack) in slacks
    if minSlack > min(real(slack), imag(slack))
      minSlack = min(real(slack), imag(slack))
      minCstrName = cstrName
    end
  end
  return minSlack, minCstrName
end
