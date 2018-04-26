##################################
## Internal operations
##################################

"""
    add_to_dict!(dict::SortedDict{Any, V}, key, val::V) where V<:Number

    *Sparsely* add `val` to the `key` entry of `dict` dictionnary. That is creates
    the entry if needed, deletes it if the resulting value is null.
"""
function add_to_dict!(dict::SortedDict{U, Number}, key::U, val::T) where T<:Number where U
    if !haskey(dict, key)
        dict[key] = 0
    end
    dict[key] += val
    if dict[key] == 0
        delete!(dict, key)
    end
end

function setindex!(pt::Point, val::Number, var::Variable)
    pt.coords[var] = val
end

## Empty constructors
Point() = Point(SortedDict{Variable, Number}())
Exponent() = Exponent(SortedDict{Variable, Degree}())
Polynomial() = Polynomial(SortedDict{Exponent, Number}())

Exponent(x::Variable) = Exponent(SortedDict{Variable, Degree}(x=>Degree(1,0)))
function Point(vars::Array{Variable}, vals::Array{<:Number})
  if length(vars) != length(vals)
    error("Point(): input arrays must have same size.")
  end

  pt = Point()
  for i=1:length(vars)
    var, val = vars[i], vals[i]
    if isreal(var) val = real(val) end
    if isbool(var) && val != 0
      val = Int((val/abs(val)+1)/2)
    end
    add_coord!(pt, vars[i], val)
  end
  return pt
end


## copy methods
function copy(p::Polynomial)
    pdict = deepcopy(p.poly)
    return Polynomial(pdict)
end

## convert functions
function convert(::Type{Polynomial}, expo::Exponent)
    return Polynomial(SortedDict{Exponent, Number}(expo=>1.0))
end

function convert(::Type{Polynomial}, x::Variable)
    return Polynomial(SortedDict{Exponent, Number}(Exponent(SortedDict(x=>Degree(1,0)))=>1.0))
end

function convert(::Type{Polynomial}, λ::Number)
    return Polynomial(SortedDict{Exponent, Number}(Exponent()=>λ))
end

function convert(::Type{Point}, pt::SortedDict{Variable, T}) where T
    return Point(convert(SortedDict{Variable, Number}, pt))
end

function convert(::Type{AbstractPolynomial}, λ::Number)
    return Polynomial(SortedDict{Exponent, Number}(Exponent()=>λ))
end

##
isreal(x::Variable) = x.kind <: Real
isbool(x::Variable) = x.kind <: Bool
iscomplex(x::Variable) = x.kind <: Complex


function conj(d::Degree)
    return Degree(d.conjvar, d.explvar)
end

function conj(x::Variable)
    if x.kind<:Complex
        return Exponent(SortedDict(x=>Degree(0,1)))
    else
        return Exponent(SortedDict(x=>Degree(1,0)))
    end
end

function conj(expo::Exponent)
    expodict = SortedDict{Variable, Degree}()
    for (var, deg) in expo.expo
        if iscomplex(var)
            expodict[var] = conj(deg)
        else
            expodict[var] = deg
        end
    end
    return Exponent(expodict)
end

function conj(p::Polynomial)
    pdict = SortedDict{Exponent, Number}()
    for (expo, λ) in p
        pdict[conj(expo)] = conj(λ)
    end
    return Polynomial(pdict)
end


## real, imag
real(p::Polynomial) = (p + conj(p))/2
imag(p::Polynomial) = (p - conj(p))/(2im)

function real(p::T) where T<:AbstractPolynomial
    return (p + conj(p))/2
end

function imag(p::T) where T<:AbstractPolynomial
    return (p + -conj(p))/(2im)
end

function abs2(p::T) where T<:AbstractPolynomial
    return p*conj(p)
end



## Evaluate
# NOTE: Point are sparse vectors, hence non referenced variables are assumed to be null.
function evaluate(p::Polynomial, pt::Point)
    res=0
    expldeg = conjdeg = 0
    for (expo, λ) in p
        res += λ*evaluate(expo, pt)
    end
    return res
end

function evaluate(expo::Exponent, pt::Point)
    res=1
    for (var, deg) in expo.expo
        if haskey(pt, var)
            val = evaluate(var, pt)
            res *= (val^deg.explvar) * (conj(val)^deg.conjvar)
        #TODO: correct evaluation
        else
            val = evaluate(var, pt)
            res *= (val^deg.explvar) * (conj(val)^deg.conjvar)
        end
    end
    return res
end

function evaluate(x::Variable, pt::Point)
    if !haskey(pt, x)
        return 0
    elseif x.kind<:Complex
        return pt[x]
    else
        return real(pt[x])
    end
end
