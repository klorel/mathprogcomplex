varname_cplx2real(varname::String) = (varname*"_Re",varname*"_Im")

function Base.show(io::IO, x::Variable)
    print(io, x.name)
end

function Base.show(io::IO, d::Degree)
    print(io, "(", d.explvar, ",", d.conjvar, ")")
end

function Base.show(io::IO, exp::Exponent)
  if exp == Exponent()
    print(io, "1")
  else
    expo = exp.expo
    i=length(expo)
    for (var, deg) in expo
      if var.kind <: Complex && deg.conjvar>0
        print(io, "conj(", var, ")")
        if deg.conjvar > 1
          print(io, "^", deg.conjvar)
        end
        if deg.explvar > 0
          print(io, " * ")
        end
      end
      if deg.explvar == 1
        print(io, var)
      elseif deg.explvar > 1
        print(io, var, "^", deg.explvar)
      end
      if i > 1
        print(io, " * ")
      end
      i -= 1
    end
  end
end

function Base.print(io::IO, exp::Exponent)
  if exp == Exponent()
    print(io, "1")
  else
    expo = exp.expo
    i=length(expo)
    for (var, deg) in expo
      if var.kind <: Complex && deg.conjvar>0
        print(io, "conj(", var, ")")
        if deg.conjvar > 1
          print(io, "^", deg.conjvar)
        end
        if deg.explvar > 0
          print(io, " * ")
        end
      end
      if deg.explvar == 1
        print(io, var)
      elseif deg.explvar > 1
        print(io, var, "^", deg.explvar)
      end
      if i > 1
        print(io, " * ")
      end
      i -= 1
    end
  end
end

function Base.print(io::IO, P::Polynomial)
  poly = P.poly
  i = length(poly)
  un = Exponent()
  for (expo, λ) in poly

    if expo == un
      print(io, "$λ")
    else
      if λ != 1
        print(io, "(", λ, ")*")
      end
      print(io, expo)
    end
    if i > 1
      print(io, " + ")
    end
    i -= 1
  end
end


function Base.print(io::IO, pt::Point)
  for (var, val) in pt
    println(io, var, " ", val)
  end
end
