function get_varsconj(exp::Exponent)
  varconj = Vector{Variable}()
  for (var, expo) in exp
    for i=1:expo.conjvar
      push!(varconj, var)
    end
  end
  varconj
end
function get_varslin(exp::Exponent)
  varlin = Vector{Variable}()
  for (var, expo) in exp
    for i=1:expo.explvar
      push!(varlin, var)
    end
  end
  varlin
end

function print_string(io, strng, len; alignright = true)
  if alignright
    print(io, " "^(len-length(strng)), strng, " ")
  else
    print(io, strng, " "^(len-length(strng)), " ")
  end
end

function print_dat_line(io, linetype, cstrname, var1, var2, val1, val2, maxvarlen, maxcstrlen)
  @printf(io, "%8s ", linetype)
  print_string(io, cstrname, maxcstrlen)
  print_string(io, var1, maxvarlen)
  print_string(io, var2, maxvarlen)
  @printf(io, "% 23.16e % 23.16e\n", val1, val2)
end

function print_quad_expo(io, expo::Exponent, cat::String, coeff, maxvarlen, maxcstrlen)
  vars_conj = get_varsconj(expo)
  vars_lin = get_varslin(expo)
  if length(vars_conj) == 1 && length(vars_lin) == 1
    print_dat_line(io, "QUAD", cat, vars_conj[1].name, vars_lin[1].name, real(coeff), imag(coeff), maxvarlen, maxcstrlen)
  elseif length(vars_conj) == 0 && length(vars_lin) == 2 && vars_lin[1].kind <: Real && vars_lin[2].kind <: Real
    print_dat_line(io, "QUAD", cat, vars_lin[1].name, vars_lin[2].name, real(coeff), imag(coeff), maxvarlen, maxcstrlen)
  elseif length(vars_conj) == 1 && length(vars_lin) == 0
    print_dat_line(io, "LIN", cat, vars_conj[1].name, "NONE", real(coeff), imag(coeff), maxvarlen, maxcstrlen)
  elseif length(vars_conj) == 0 && length(vars_lin) == 1
    print_dat_line(io, "LIN", cat, "NONE", vars_lin[1].name, real(coeff), imag(coeff), maxvarlen, maxcstrlen)
  elseif length(vars_conj) == 0 && length(vars_lin) == 0
    print_dat_line(io, "CONST", cat, "NONE", "NONE", real(coeff), imag(coeff), maxvarlen, maxcstrlen)
  else
    warn("export_to_dat(): Exponent $expo not supported.")
  end
end

function print_constraint(io::IO, cstrname::String, cstr::Constraint, maxvarlen::Int, maxcstrlen::Int, expos::Dict{Exponent, String})
  print_poly!(io, cstr.p, cstrname, maxvarlen, maxcstrlen, expos)

  if ismatch(r"_Re", cstrname)
    if real(cstr.lb) != -Inf
      print_dat_line(io, "LB", cstrname, "NONE", "NONE", real(cstr.lb), imag(cstr.lb), maxvarlen, maxcstrlen)
    end
    if real(cstr.ub) != Inf
      print_dat_line(io, "UB", cstrname, "NONE", "NONE", real(cstr.ub), imag(cstr.ub), maxvarlen, maxcstrlen)
    end
  elseif ismatch(r"_Im", cstrname)
    if imag(cstr.lb) != -Inf
      print_dat_line(io, "LB", cstrname, "NONE", "NONE", real(cstr.lb), imag(cstr.lb), maxvarlen, maxcstrlen)
    end
    if imag(cstr.ub) != Inf
      print_dat_line(io, "UB", cstrname, "NONE", "NONE", real(cstr.ub), imag(cstr.ub), maxvarlen, maxcstrlen)
    end
  else
    if imag(cstr.lb) != -Inf-im*Inf
      print_dat_line(io, "LB", cstrname, "NONE", "NONE", real(cstr.lb), imag(cstr.lb), maxvarlen, maxcstrlen)
    end
    if imag(cstr.ub) != +Inf+im*Inf
      print_dat_line(io, "UB", cstrname, "NONE", "NONE", real(cstr.ub), imag(cstr.ub), maxvarlen, maxcstrlen)
    end
  end
end
nb_from_str(string::String) = parse(matchall(r"\d+", string)[1])
get_scenario(string::String) = String(split(string, "_")[1])

function print_doc(io, filename)
  println(io, "# $filename - exported $(now())")
  println(io, "# Description of a QCQP optimization problem. 3 sections:")
  println(io, "#  - List of all variables in the problem, along with type, name, and possibly complex value.")
  println(io, "#    Grouped by \"VAR_TYPE\" tag (col 1). Type is \"C\" Complex, \"R\" Real, \"BOOL\" Boolean (col 2). ")
  println(io, "#    Name (col 3). Real and imag part of value (col 5 and 6 resp.).")
  println(io, "#  - Description of the objective function (one quadratic polynomial).")
  println(io, "#    Sum of monomials, grouped by \"OBJ\" (col 2), either:")
  println(io, "#                order 0 (\"CONST\" col 1) : real and imag part of coef. resp. at col 5 and 6.")
  println(io, "#                order 1 (\"LIN\" col 1) : real and imag part of coef. resp. at col 5 and 6, either:")
  println(io, "#                         variable name in col 3 hence conjugate variable in polynomial,")
  println(io, "#                         variable name in col 4 hence variable in polynomial.")
  println(io, "#                order 2 (\"QUAD\" col 1) : real and imag part of coef. resp. at col 5 and 6, ")
  println(io, "#                         monomial is product of conjugate variable at col 3 and variable at col 4.")
  println(io, "#  - Description of the constraints : one polynomial and two complex numeric bounds.")
  println(io, "#    Constraints are grouped by name (col 2). Quadratic body is described as the objective.")
  println(io, "#    Lower bound (\"LB\" col 1) and upper bound (\"UB\" col 1) have their complex value in col 5 and 6.")
  println(io, "#")
end

"""
  print_variables(io::IO, variables, pt::Point)

  Write a .dat description of `variables` variables to `io`.
"""
function print_variables(io::IO, variables, pt::Point, maxvarlen, maxcstrlen)
  for varname in sort(collect(keys(variables)))
    var = Variable(varname, variables[varname])
    if iscomplex(var)
      var_type = "CPLX"
    elseif isbool(var)
      var_type = "BOOL"
    elseif isreal(var)
      var_type = "REAL"
    else
      error("Export_to_dat(): unsuported variable type $(var.kind)")
    end
    val = 0

    if haskey(pt, var)
      val = pt[var]
    end

    print_dat_line(io, "VAR_TYPE", var_type, var.name, "NONE", real(val), imag(val), maxvarlen, maxcstrlen)
  end
end

"""
  print_poly!(io::IO, p::AbstractPolynomial, cat::String, maxvarlen, maxcstrlen, expos::Dict{Exponent, String})

Print the `p` polynomial corresponding to the constraint or objective
`cat` (category) to `io`. Each of the polynomial's exponent is printed in a line,
either explicitly if its degree is 2 or less, or implicitly by defining an
exponent name in the `expos` dict if it required, and print the exponent name
with the coefficient (which essentially is a linear term).
"""
function print_poly!(io::IO, p::AbstractPolynomial, cat::String, maxvarlen, maxcstrlen, expos::Dict{Exponent, String})
  constval = 0

  for expo in sort(collect(keys(p)))
    coeff = p[expo]

    exp_globaldeg = expo.degree.explvar + expo.degree.conjvar
    if exp_globaldeg > 2
      if !haskey(expos, expo)
        expos[expo] = "MONO_$(length(expos))"
      end
      print_dat_line(io, "MONO", cat, expos[expo], "NONE", real(coeff), imag(coeff), maxvarlen, maxcstrlen)
    elseif exp_globaldeg == 0
      constval = coeff
    else
      print_quad_expo(io, expo, cat, coeff, maxvarlen, maxcstrlen)
    end
  end
  if constval != 0
    print_dat_line(io, "CONST", cat, "NONE", "NONE", real(constval), imag(constval), maxvarlen, maxcstrlen)
  end
end

"""
  export_to_dat(pb_optim::Problem, outpath::String, pt::Point = Point())

Write the `pb_optim` problem to the `outpath` folder, with the dispensory initial
point `pt`. Output files are `real_minlp_instance.dat` and
`real_minlp_precond_cstrs.dat`.
"""
function export_to_dat(pb_optim::Problem, outpath::String, pt::Point = Point())
  ## Get max length varname
  maxvarlen = -1
  for var in pb_optim.variables
    if length(var[1]) > maxvarlen
      maxvarlen = length(var[1])
    end
  end

  ## Get max length dat constraint name
  maxcstrlen = -1
  for cstrname in keys(pb_optim.constraints)
    if length(cstrname) > maxcstrlen
      maxcstrlen = length(cstrname)
    end
  end

  # Container for monomials definition, to be written lastly
  expos = Dict{Exponent, String}()
  precond_cstrs = Set{String}()

  isdir(outpath) || mkpath(outpath)
  filename = joinpath(outpath, "real_minlp_instance.dat")
  touch(filename)
  outfile = open(filename, "w")

  ## Comment section
  print_doc(outfile, filename)

  ## Print variables
  variables = get_variables(pb_optim)
  print_variables(outfile, variables, pt, maxvarlen, maxcstrlen)

  ## Print objective
  print_poly!(outfile, pb_optim.objective, "OBJ", maxvarlen, maxcstrlen, expos)

  ## Print constraints
  for cstrname in sort(collect(keys(pb_optim.constraints)))
    cstr = pb_optim.constraints[cstrname]
    print_constraint(outfile, cstrname, cstr, maxvarlen, maxcstrlen, expos)

    if cstr.precond != :none
      push!(precond_cstrs, cstrname)
    end
  end

  ## Print collected monomials definition
  for (expo, exponame) in expos
    for (var, deg) in expo
      print_dat_line(outfile, "MONO_DEF", exponame, var.name, "NONE", deg.explvar, deg.conjvar, maxvarlen, maxcstrlen)
    end
  end
  close(outfile)

  ## Print constraints with preconditionning
  filename = joinpath(outpath, "real_minlp_precond_cstrs.dat")
  touch(filename)
  outfile = open(filename, "w")

  print_string(outfile, "#cstrname", maxcstrlen)
  @printf(outfile, "%10s\n", "Precondtype")
  for cstrname in sort(collect(precond_cstrs))
    if get_constraint(pb_optim, cstrname).precond == :sqrt
      print_string(outfile, cstrname, maxcstrlen)
      @printf(outfile, "%10s\n", "SQRT")
    else
      warn("export_dat(): Unknown preconditionning for cstr $cstrname.")
    end
  end
  close(outfile)
end




###
###
function export_matpower_to_dat(QCQP::Problem, filename::String, pt::Point = Point())
  ## Get max length varname
  maxvarlen = -1
  for var in QCQP.variables
    if length(var.name) > maxvarlen
      maxvarlen = length(var.name)
    end
  end

  ## Sort constraints by type (voltm, unit and rest), and build dat constraint name
  cstrs_keys = sort(collect(keys(QCQP.constraints)))

  cstr_keys = Set(keys(QCQP.constraints))
  voltm_keys = filter(x->ismatch(r"VOLTM", x), cstr_keys)
  unit_keys = filter(x->ismatch(r"UNIT", x), cstr_keys)
  load_keys = setdiff(cstr_keys, union(unit_keys, voltm_keys))
  id_to_loadkey = Dict(nb_from_str(str)=> (str, "LOAD_$(string(nb_from_str(str)))") for str in load_keys)
  id_to_voltmkey = Dict(nb_from_str(str)=> (str, String(split(str, "_")[2])) for str in voltm_keys)
  id_to_unitkey = Dict(nb_from_str(str)=> (str, "UNIT_$(string(nb_from_str(str)))") for str in unit_keys)

  ## Get max length dat constraint name
  maxcstrlen = -1
  for (_, val) in union(id_to_loadkey, id_to_voltmkey, id_to_unitkey)
    if length(val[2]) > maxcstrlen
      maxcstrlen = length(val[2])
    end
  end

  touch(filename)
  outfile = open(filename, "w")

  ## Comment section
  print_doc(outfile, filename)

  ## Print variables
  variables = get_variables(QCQP)
  print_variables(outfile, variables, pt, maxvarlen, maxcstrlen)

  ## Print objective
  obj_keys = sort(collect(keys(QCQP.objective)))
  const_printed = false
  const_val = 0
  for expo in obj_keys
      coeff = QCQP.objective[expo]

      if length(expo) != 0
        print_quad_expo(outfile, expo, "OBJ", coeff, maxvarlen, maxcstrlen)
      else
        const_val = coeff
      end
  end
  print_dat_line(outfile, "CONST", "OBJ", "NONE", "NONE", real(const_val), imag(const_val), maxvarlen, maxcstrlen)

  for i in sort(collect(keys(id_to_loadkey)))
    cstr = QCQP.constraints[id_to_loadkey[i][1]]
    print_constraint(outfile, id_to_loadkey[i][2], cstr, maxvarlen, maxcstrlen, expos)
  end
  for i in sort(collect(keys(id_to_voltmkey)))
    cstr = QCQP.constraints[id_to_voltmkey[i][1]]
    print_constraint(outfile, id_to_voltmkey[i][2], cstr, maxvarlen, maxcstrlen, expos)
  end
  for i in sort(collect(keys(id_to_unitkey)))
    cstr = QCQP.constraints[id_to_unitkey[i][1]]
    print_constraint(outfile, id_to_unitkey[i][2], cstr, maxvarlen, maxcstrlen, expos)
  end
  close(outfile)
end
