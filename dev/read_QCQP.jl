include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))

## This should be useless.
# TODO: remove
# function add_to_dict!(dict, key, val)
#   if !haskey(dict, key)
#     dict[key] = 0
#   end
#   dict[key] += val
# end


"""
pb = read_QCQP(filepath::String)

Read a QCQP file (e.g. WB2.dat) and builds a corresponding Problem type object.
"""
function read_QCQP(filepath::String)
  pb = Problem()

  data = readdlm(filepath)
  line_ind = 1
  while data[line_ind, 1] == "VAR_TYPE"
    add_variable!(pb, Variable(data[line_ind, 3], Complex))
    line_ind += 1
  end

  obj = Polynomial()
  while data[line_ind, 2] == "OBJ"
    obj = obj + build_monomial_from_line(data[line_ind, :], pb)
    line_ind += 1
  end
  set_objective!(pb, obj)

  ## Add constraints
  while line_ind < size(data, 1)
    line = data[line_ind, :]
    current_cstrName = String(line[2])
    current_p = Polynomial()
    current_lb = 0
    current_ub = 0
    while string(line[2]) == current_cstrName && line_ind <= size(data, 1)
      if line[1] == "UB"
        current_ub = line[5] + im*line[6]
      elseif line[1] == "LB"
        current_lb = line[5] + im*line[6]
      else
        current_p = current_p + build_monomial_from_line(line, pb)
      end
      line_ind += 1
      if line_ind <= size(data, 1)
        line = data[line_ind, :]
      end
    end
    add_constraint!(pb, current_cstrName, current_lb << current_p << current_ub)
  end
  pb
end

function build_monomial_from_line(line::Vector{Any}, pb::Problem)
  var1 = read_var(String(line[3]), pb)
  var2 = read_var(String(line[4]), pb)
  (line[5] + im * line[6]) * conj(var1) * var2
end


"""
var = read_var(varName, pb)

Return the Problem variable corresponding to varName if it exists in pb, or else
warn, add it to pb and return the variable.
"""
function read_var(varName::String, pb::Problem)
  var = Variable()
  if varName == "NONE"
    var = Variable("1", Complex)
  elseif !has_variable(pb, varName)
    var = Variable(varName, Complex)
    warn("Variable ", varName, " is not known in the current Problem. Adding ", var, ".")
    add_variable!(pb, var)
  else
    var = Variable(varName, pb.variables[varName])
  end
  var
end


extractWord(l, x) = String(matchall(r"\S+", l)[x])
