"""
  add_coord!(pt::Point, var::Variable, val::Number)

  Sparsely add `val` to the `var` coordinate of `pt`.
"""
function add_coord!(pt::Point, var::Variable, val::Number)
  return add_to_dict!(pt.coords, var, val)
end

## Overloads for transparent iteration on `coords` attribute
start(pt::Point) = start(pt.coords)
next(pt::Point, state) = next(pt.coords, state)
done(pt::Point, state) = done(pt.coords, state)
length(pt::Point) = length(pt.coords)
haskey(pt::Point, key) = haskey(pt.coords, key)
keys(pt::Point) = keys(pt.coords)
values(pt::Point) = values(pt.coords)
getindex(pt::Point, var::Variable) = pt.coords[var]
function setindex!(pt::Point, var::Variable, val::Number)
  if val != 0
    setindex!(pt.coords, var, val)
  end
  return
end


function add!(pt1::Point, pt2::Point)
  for (var, val) in pt2.coords
    add_coord!(pt1, var, val)
  end
end

function add(pt1::Point, pt2::Point)
  pt_coord = copy(pt1.coords)
  pt = Point(pt_coord)
  add!(pt, pt2)
  return pt
end

+(pt1::Point, pt2::Point) = add(pt1, pt2)
-(pt1::Point, pt2::Point) = pt1 + (-1)*pt2

function merge(x::Point ...)
  pt_res = Point()
  for pt in x
    add!(pt_res, pt)
  end
  return pt_res
end

function *(pt1::Point, λ::Number)
  pt = Point()
  if λ == 0
    return pt
  end
  for (var, val) in pt1
    pt[var] = λ*val
  end
  return pt
end
*(λ::Number, pt1::Point) = pt1*λ


function norm(pt::Point, p::Real = 2)
    p > 0 || error("norm(): p should be positive ($p ≤ 0).")
    if p == Inf
        maxi = -Inf
        for (var, val) in pt
            # iscomplex(var) && warn("norm(): var $(var) is not real, splitting real and imag parts.")
            if !iscomplex(var) && (val > maxi)
                maxi = val
            elseif iscomplex(var) && max(real(val), imag(val)) > maxi
                maxi = max(real(val), imag(val))
            end
        end
        return maxi
    elseif p == 1
      s = 0
      for (var, val) in pt
          if !iscomplex(var)
            s += abs(val)
          else
            # warn("norm(): var $(var) is not real, splitting real and imag parts.")
            s += abs(real(val)) + abs(imag(val))
          end
      end
      return s
    else
      s = 0
      for (var, val) in pt
          if !iscomplex(var)
            s += val^p
          else
            # warn("norm(): var $(var) is not real, splitting real and imag parts.")
            s += real(val)^p + imag(val)^p
          end
      end
      return s^(1/p)
  end
end
