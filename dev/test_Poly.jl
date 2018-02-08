"""
    test_Poly
    
"""
ROOT = pwd()
include(joinpath(ROOT, "src_PolynomialOptim", "PolynomialOptim.jl"))

x = Variable("x", Complex)
y = Variable("y", Complex)
z = Variable("z", Real)
b = Variable("b", Bool)

## Polynomial algebra
pt = Point(Dict(x=>2, y=>1+im, z=>0+im, b=>2.1))
print(pt)

p1 = x^2 + 3*y + conj(x) + conj(z) + b
println(p1)
evaluate(p1, pt) == 10 + 3im

p2 = y^6 - y^6 + (-y*x*b) / 4 + π*conj(b)
println(p2)


## Poly operators
@assert evaluate(x, pt) == 2
@assert evaluate(y, pt) == 1+im
@assert evaluate(conj(y), pt) == 1-im
@assert evaluate(z, pt) == 0
@assert evaluate(b, pt) == 1

@assert evaluate(x*y, pt) == 2+2im
@assert evaluate(p2, pt) == (-1-im)/2 + π

@assert evaluate(real(y), pt) == 1
@assert evaluate(imag(y), pt) == 1

@assert evaluate(x*y*conj(y), pt) == 2*(1+im)*(1-im)
@assert evaluate(x*y + conj(y), pt) == 3+im
@assert evaluate(x*y + 1 + 1im, pt) == 3+3im



## Point algebra
pt = Point(Dict(x=>2, y=>1+im, z=>0+im, b=>2.1))
pt1 = Point(Dict(z=>0+im))
pt2 = Point(Dict(x=>3, y=>-1+2im, b=>-1))
@assert merge(pt, pt1, pt2) == Point(Dict(x=>5, y=>3im, z=>0+im, b=>1))

@assert pt + pt1 - 3*pt2 == Point(Dict(x=>-7, y=>4-5im, z=>0+im, b=>1))

@assert norm(pt2, Inf) == 3
@assert norm(pt2, 1) == 6
@assert norm(pt2, 2) == √(14)


## Poly cplx2real
p = x + y
println("p: $p")

p1 = copy(p)
println("p1: $p1")

add!(p1, y*z)
println("p: $p")
println("p1: $p1")

real_p, imag_p = cplx2real(p1)


expo = Exponent(Dict(x=>Degree(1,0), y=>Degree(1,0), z=>Degree(2,0)))
println(expo)
real_p, imag_p = cplx2real(expo)
println(expo)

println(real_p)
println(imag_p)


println("p1: $p1")
println("real part: $real_p")
println("imag part: $imag_p")


### Exponent generation...
x = Variable("x", Complex)
y = Variable("y", Complex)
z = Variable("z", Complex)
μ = Variable("μ", Complex)

input_vars = [x, y, z, μ]
d = 3
