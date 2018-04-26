# module Poly
using DataStructures

import Base: ==, !=, isless, isconst, isreal, isnull, isequal
import Base: +, -, *, /, ^, conj, conj!, abs2, norm, real, imag
import Base: show, print, convert, copy, hash, merge
import Base: start, next, done, length, setindex!, getindex, haskey, keys, values


abstract type  AbstractPolynomial end


"""
    Degree(explvar::Int, conjvar::Int)

Define a mathematical degree, which is used to define a variable exponent to
`explvar` and its conjugate exponent to `conjvar`.

### Attributes
- `explvar`: Int
- `conjvar`: Int
"""
mutable struct Degree
    explvar::Int
    conjvar::Int
end


"""
    Variable(varname::String, Complex)

Define a mathematical variable `varname` of a certain mathematical kind,
`Complex` here.

### Attributes
- `name`: String
- `kind`: a type among `Complex`, `Real` and `Bool`
"""
mutable struct Variable <: AbstractPolynomial
    name::String
    kind::Type

    function Variable(name, kind)
        if kind ∉ Set([Complex, Real, Bool]) || typeof(name) ∉ Set([String, SubString{String}])
            error("Variable() : attempting to define a variable $name of type $kind, supported types are {Complex, Real, Bool}")
        end
        return new(String(name), kind)
    end
end

# function isless(t1::Type, t2::Type)
#     if t1 == t2
#         return false
#     end
#     if t1 == Real
#         if t2 == Bool
#             return false
#         elseif t2 == Complex
#             return true
#         end
#     elseif t1 == Bool
#         if t2 == Real
#             return true
#         elseif t2 == Complex
#             return true
#         else
#             println(t1, t2)
#         end
#     elseif t1 == Complex
#         if t2 == Bool
#             return false
#         elseif t2 == Real
#             return false
#         else
#             println(t1, t2)
#         end
#     else
#         println(t1, t2)
#     end
# end



"""
    Exponent(expo::SortedDict{Variable, Degree})

Define a mathematical exponent, that is a product of `Variable` and conjugated
`Variable`.

### Attributes
- `expo` : a dictionary associating `Variable` objects to a `Degree` objects,
their exponent and the exponent of their conjugates.
- `degree`: a `Degree` object indicating the global degree of the exponent :
the sum of the explicit variables degrees and of the conjugated variables
degrees.

### Exemple
```julia
julia > a, b = Variable("a", Complex), Variable("b", Complex)
julia > Exponent(SortedDict(a=>Degree(1,0), b=>Degree(1,1))) == a*b*conj(b)
```
"""
struct Exponent <: AbstractPolynomial
    expo::SortedDict{Variable, Degree}
    degree::Degree

    function Exponent(expo::SortedDict{Variable, Degree})
        degexpl, degconj = 0,0
        for (var, degree) in expo
            ((degree.explvar < 0) || (degree.conjvar < 0)) && error("Exponent(): Expected non negative exponent for variable $var (got $degree)")
            (isreal(var) && degree.conjvar != 0) && error("Exponent(): Expected nul conj exponent for real variable $var (got $degree)")
            (isbool(var) && degree.explvar ∉ SortedSet([0,1])) && error("Exponent(): Expected boolean exponent for bool $var (got $degree)")
            if degree != Degree(0,0)
                degexpl = max(degexpl, degree.explvar)
                degconj = max(degconj, degree.conjvar)
            else
                delete!(expo, var)
            end
        end
        return new(expo, Degree(degexpl, degconj))
    end
end


"""
    Polynomial(poly::SortedDict{Exponent, Number}, degree::Degree)

Define a mathematical polynomial, that is a linear combinaison of product of
`Variable` and conjugated `Variable`.

### Attributes
- `poly` : a dictionary associating an `Exponent` to a a number.
- `degree`: a `Degree` object indicating the global degree of the polynomial :
the maximum of the explicit exponent degrees and of the conjugated exponent
degrees.

### Exemple
```julia
julia > a, b = Variable("a", Complex), Variable("b", Complex)
julia > p = a*b*conj(b) + (2+3im) * b^3
julia > p.poly = SortedDict(a*b*conj(b)=>1, b^3=>2+3im)
julia > p.degree = (3,1)
```
"""
struct Polynomial <: AbstractPolynomial
    poly::SortedDict{Exponent, Number}
    degree::Degree

    function Polynomial(poly::SortedDict{Exponent, Number})
        degexpl, degconj = 0, 0
        for (expo, λ) in poly
            if λ!=0
                degexpl = max(degexpl, expo.degree.explvar)
                degconj = max(degconj, expo.degree.conjvar)
            else
                delete!(poly, expo)
            end
        end
        return new(poly, Degree(degexpl, degconj))
    end
end

"""
    Point(coords::SortedDict{Variable, Number})

Define a mathematical point, that is a pairing of variables and numbers.

### Attributes
- `coords` : a dictionary associating a `Variable` to a a number.

### Exemple
```julia
julia > a, b = Variable("a", Complex), Variable("b", Complex)
julia > pt = Point(SortedDict(a=>π, b=>e+7im))
```
"""
struct Point
    coords::SortedDict{Variable, Number}

    function Point(dict::SortedDict)
        dict_pt = SortedDict{Variable, Number}()
        for (var, val) in dict
            if !isa(var, Variable) || !isa(val, Number)
                error("Point(): Expected pair of (Variable, Number), got ($var, $val) of type ($typeof(var), $typeof(val)) instead.")
            end
            if isbool(var)
                booled = 0
                if val != 0
                    booled = Int((val/abs(val) + 1) / 2)
                end
                (booled ≠ val) && warn("Point(): $var is $(var.kind), provided value is $val, $booled stored.")
                add_to_dict!(dict_pt, var, booled)
            elseif isreal(var)
                realed = real(val)
                (realed ≠ val) && warn("Point(): $var is $(var.kind), provided value is $val, $realed stored.")
                add_to_dict!(dict_pt, var, realed)
            else
                add_to_dict!(dict_pt, var, val)
            end
        end
        return new(dict_pt)
    end
end


"""
    Constraint(p::AbstractPolynomial, lb::Number, ub::Number)

Define a mathematical complex constraint, from a polynomial body and two complex
bounds. Equivalent to the two real inequality produced from real and imag part
of `p`, `lb` and `ub`.

### Attributes
- `p::Polynomial` : the polynomial body of the constraint.
- `lb::Number` : complex lower bound.
- `ub::Number` : complex upper bound.
- `precond::Symbol` : a symbol describing the preconditioning method to be
applied for this constraint. Default value is `:none`.

### Exemple
```julia
julia > a, b = Variable("a", Complex), Variable("b", Complex)
julia > cstr = abs2(a^2) + abs2(3*b) << 25
julia > cstr.precond = :sqrt
```
"""
mutable struct Constraint
  p::Polynomial
  lb::Number
  ub::Number
  precond::Symbol

  Constraint(p::AbstractPolynomial, lb::Number, ub::Number) = new(p, lb, ub, :none)
end


"""
    Problem()

Define an empty mathematical polynomial optimization problem, characterized by
a polynomial `objective`, a dictionnary of `constraints` and a dictionnary of
used variables `variables`.

### Attributes
- `objective::Polynomial` : polynomial criterion.
- `constraints::SortedDict{String, Constraint}` : dictionnary of constraint name to
constraint.
- `variables::SortedDict{String, Type}` : dictionnary of variable name to type.

### Exemple
```julia
julia > x, y, z = ...
julia > pb = Problem()
julia > set_objective(pb, abs2(x+y+z))
julia > add_constraint!(pb, abs2(x) << 1)
julia > add_constraint!(pb, y << 0.5+1im)
julia > add_constraint!(pb, 0-1im << z << 1+0im)
```
"""
mutable struct Problem
  objective::Polynomial
  constraints::SortedDict{String, Constraint}
  variables::SortedDict{String, Type}
end


include("Modeler_accessors.jl")
include("Modeler_constructors.jl")
include("Modeler_cplx2real.jl")
include("Poly_Cplx2Real.jl")
include("Poly_operators.jl")
include("Poly_order.jl")
include("Poly_point_algebra.jl")
include("Poly_poly_algebra.jl")
include("Poly_print.jl")
include("utils_ampl.jl")
include("utils_dat_compare.jl")
include("utils_dat_export.jl")
include("utils_dat_import.jl")
include("utils_Poly.jl")
include("utils_jump.jl")

# export Variable, Point, Exponent, Monomial, Polynomial
# export isconst, isone, is_homogeneous
# export evaluate, abs2, conj, add!, add
# export cplx2real, real, imag
# export norm

#end
