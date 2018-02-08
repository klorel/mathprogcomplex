varname_cplx2real(varname::String) = (varname*"_Re",varname*"_Im")

function Base.show(io::IO, x::Variable)
    print(io, x.name)
end

function Base.show(io::IO, d::Degree)
    print(io, "(", d.explvar, ",", d.conjvar, ")")
end

function Base.print(io::IO, exp::Exponent)
  if exp == Exponent()
    print(io, "1")
  else
    expo = exp.expo
    sortedCollec = sort(collect(expo), by=(x)->x[1].name)
    i=length(sortedCollec)
    for (var, deg) in sortedCollec
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
  sorted_keys = sort(collect(keys(poly)))
  for expo in sorted_keys
    λ = poly[expo]
    if λ != 0
      print(io, "(", λ, ")")
    end
    if expo.degree != Degree(0,0)
      print(io, "*")
      print(io, expo)
      if i > 1
        print(io, " + ")
      end
      i -= 1
    end
  end
end


function Base.print(io::IO, pt::Point)
  for (var, val) in sort(collect(pt), by=x->x[1].name)
    println(io, var, " ", val)
  end
end
