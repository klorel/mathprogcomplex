# Conversion of all complex variables to real ones in the Poly structures

"""
  realPart, imagPart = cplx2real(expo::Exponent)

  Convert a complex Exponent in complex variables into `realPart` and
  `imagPart` polynomials of twice as many variables, real and imag parts of
  `expo` variables. Done recursively with the `cplx2real_rec` function.
"""
function cplx2real(expo::Exponent)
  vars, inds = collect(keys(expo.expo)), collect(values(expo.expo))
  return cplx2real_rec(vars, inds, Polynomial()+1, Polynomial()+0, length(expo)+1, Degree(0,0))
end

"""
realPart, imagPart = cplx2real_rec(vars::Array{Variable}, degs::Array{Degree}, realPart::Polynomial, imagPart::Polynomial, cur_ind::Int, cur_deg::Degree)

Transform recursively the complex exponent represented by the `vars` and `degs`
arrays into its real and imag parts, functions of its imag and real part
variables.
`cur_ind` decreases to 0, `cur_deg` decreases to Degree(0,0) for each step of
cur_ind. Terminaison case is reached at 0, Degree(0,0).
Initial arrays `vars` and `degs` are read only.

### Arguments
- vars::Array{Variable}
- degs::Array{Degree}
- realPart::Polynomial
- imagPart::Polynomial
- cur_ind::Int
- cur_deg::Degree
"""
function cplx2real_rec(vars::Array{Variable}, degs::Array{Degree}, realPart::Polynomial, imagPart::Polynomial, cur_ind::Int, cur_deg::Degree)
  ## Final case:
  if cur_ind == 1 && cur_deg == Degree(0,0)
    return (realPart, imagPart)
  ## One less variable to deal with:
  elseif cur_deg == Degree(0,0)
    return cplx2real_rec(vars, degs, realPart, imagPart, cur_ind-1, degs[cur_ind-1])
  ## Recursion rule, decrease the current variable exponent until it reaches Degree(0,0):
  else
    var = vars[cur_ind]
    if iscomplex(var)
      var_real, var_imag = varname_cplx2real(var.name)
      var_R, var_I = Variable(var_real, Real), Variable(var_imag, Real)
      if cur_deg.explvar > 0
        return cplx2real_rec(vars, degs, var_R * realPart - var_I * imagPart, var_R * imagPart + var_I * realPart, cur_ind, Degree(cur_deg.explvar-1, cur_deg.conjvar))
      elseif cur_deg.conjvar > 0
        cur_deg.explvar == 0 || warn("cur_deg.explvar should be 0 (and not $(cur_deg.explvar)), set to this value")
        return cplx2real_rec(vars, degs, var_R * realPart + var_I * imagPart, var_R * imagPart - var_I * realPart, cur_ind, Degree(0, cur_deg.conjvar-1))
      end
    elseif isbool(var)
      return cplx2real_rec(vars, degs, var*realPart, var*imagPart, cur_ind, Degree(0,0))
    else
      return cplx2real_rec(vars, degs, Exponent(SortedDict(var=>cur_deg))*realPart, Exponent(SortedDict(var=>cur_deg))*imagPart, cur_ind, Degree(0,0))
    end
  end
end

"""
  realPart, imagPart = cplx2real(pol::Polynomial)

  Convert a complex polynomial in complex variables into `realPart` and
  `imagPart` polynomials of twice as many variables, real and imag parts of
  `pol` variables.
"""
function cplx2real(pol::Polynomial)
  realPart = Polynomial()
  imagPart = Polynomial()

  for (expo, λ) in pol
    realexpo, imagexpo = cplx2real(expo)

    realPart += realexpo*real(λ) - imagexpo*imag(λ)
    imagPart += imagexpo*real(λ) + realexpo*imag(λ)
  end
  return (realPart, imagPart)
end


function cplx2real(pt_C::Point)
  pt = Point()
  for (var, val) in pt_C
    if var.kind <: Complex
      var_real, var_imag = varname_cplx2real(var.name)
      if real(val) != 0
        pt[Variable(var_real, Real)] = real(val)
      end
      if imag(val) != 0
        pt[Variable(var_imag, Real)] = imag(val)
      end
    else
      if val != 0
        pt[var] = real(val)
      end
    end
  end
  return pt
end


function real2cplx(pt::Point)
    ptC = Point()
    for (var, val) in pt
        if ismatch(r"(_Re|_Im)", var.name)
            if ismatch(r"_Re", var.name)
                varname = var.name[1:(end-3)]
                ptC[Variable(varname, Complex)] = val + im*pt[Variable(varname*"_Im", Real)]
            end
        else
            ptC[var] = val
        end
    end
    return ptC
end
