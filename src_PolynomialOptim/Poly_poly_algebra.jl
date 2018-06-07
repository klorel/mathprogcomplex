## Polynomial iterator
start(poly::Polynomial) = start(poly.poly)
next(poly::Polynomial, state) = next(poly.poly, state)
done(poly::Polynomial, state) = done(poly.poly, state)
length(poly::Polynomial) = length(poly.poly)
haskey(poly::Polynomial, key) = haskey(poly.poly, key)
keys(poly::Polynomial) = keys(poly.poly)
values(poly::Polynomial) = values(poly.poly)
getindex(poly::Polynomial, expo::Exponent) = poly.poly[expo]

## Exponent iterator
start(expo::Exponent) = start(expo.expo)
next(expo::Exponent, state) = next(expo.expo, state)
done(expo::Exponent, state) = done(expo.expo, state)
length(expo::Exponent) = length(expo.expo)
haskey(expo::Exponent, key) = haskey(expo.expo, key)
keys(expo::Exponent) = keys(expo.expo)
values(expo::Exponent) = values(expo.expo)
getindex(expo::Exponent, var::Variable) = expo.expo[var]

##################################
## Addition
##################################
function add!(p::Polynomial, p1::T) where T<:AbstractPolynomial
    return add!(p, convert(Polynomial, p1))
end

function add!(p::Polynomial, p1::Polynomial)
    for (cur_expo, λ) in p1
        λ != 0 || continue
        add_to_dict!(p.poly, cur_expo, λ)
    end
    p.degree.explvar = max(p.degree.explvar, p1.degree.explvar)
    p.degree.conjvar = max(p.degree.conjvar, p1.degree.conjvar)
end

function add(p1::Polynomial, p2::Polynomial)
    p = deepcopy(p1)
    add!(p, p2)
    return p
end

function +(p1::T, p2::U) where T<:AbstractPolynomial where U<:AbstractPolynomial
    return add(convert(Polynomial, p1), convert(Polynomial, p2))
end
function +(p1::Number, p2::T) where T<:AbstractPolynomial
    return add(convert(Polynomial, p1), convert(Polynomial, p2))
end
function +(p1::T, p2::Number) where T<:AbstractPolynomial
    return add(convert(Polynomial, p1), convert(Polynomial, p2))
end

##################################
## Product
##################################
function product(p1::Polynomial, p2::Polynomial)
    p = Polynomial()
    for (expo1, λ1) in p1
        λ1 != 0 || continue
        for (expo2, λ2) in p2
            λ2 != 0 || continue
            expoprod = product(expo1, expo2)
            add_to_dict!(p.poly, expoprod, λ1 * λ2)
        end
    end
    p.degree.explvar = p1.degree.explvar+p2.degree.explvar
    p.degree.conjvar = p1.degree.conjvar+p2.degree.conjvar
    return p
end

"""
    product!(expod, expod1)

    Add the `expod1` degree to `expod` inplace (equivalent to the monomial product)
"""
function product!(expod::Exponent, expod1::Exponent)
    for (var, deg) in expod1
        if !haskey(expod, var)
            expod.expo[var] = Degree(0,0)
        end
        expod.expo[var].explvar += deg.explvar
        expod.expo[var].conjvar += deg.conjvar
        if (expod.expo[var].explvar, expod.expo[var].conjvar) == (0,0)
            delete!(expod, var)
        end
    end
    update_degree!(expod)
    return expod
end

function product(exp1::Exponent, exp2::Exponent)
    expod = Exponent()
    product!(expod, exp1)
    product!(expod, exp2)
    return Exponent(expod)
end

function *(p1::T, p2::U) where T<:AbstractPolynomial where U<:AbstractPolynomial
    return product(convert(Polynomial, p1), convert(Polynomial, p2))
end
function *(p1::Number, p2::T) where T<:AbstractPolynomial
    return product(convert(Polynomial, p1), convert(Polynomial, p2))
end
function *(p1::T, p2::Number) where T<:AbstractPolynomial
    return product(convert(Polynomial, p1), convert(Polynomial, p2))
end

### Soustraction
function -(p::T) where T<:AbstractPolynomial
    return -1 * convert(Polynomial, p)
end

function -(p1::T, p2::U) where T<:AbstractPolynomial where U<:AbstractPolynomial
    return p1 + (-p2)
end
function -(p1::Number, p2::T) where T<:AbstractPolynomial
    return p1 + (-p2)
end
function -(p1::T, p2::Number) where T<:AbstractPolynomial
    return p1 + (-p2)
end

## Division
function divide(p1::Polynomial, p2::Polynomial)
    if length(p2) != 1
        error("/(::Polynomial, ::Polynomial): Only allowed for monomial divisor ($(length(p2))-monomial polynomial here).")
    end
    expo, λ = collect(p2)[1]
    if expo.degree != Degree(0,0)
        error("/(::Polynomial, ::Polynomial): Only allowed for constant divisor ($(expo.degree)-degree monomial here).")
    end
    if λ == 0
        error("/(::Polynomial, ::Polynomial): Only allowed for non null constant divisor.")
    end
    return p1 * (1/λ)
end

function /(p1::T, p2::U) where T<:AbstractPolynomial where U<:AbstractPolynomial
    return divide(convert(Polynomial, p1), convert(Polynomial, p2))
end
function /(p1::Number, p2::T) where T<:AbstractPolynomial
    return divide(convert(Polynomial, p1), convert(Polynomial, p2))
end
function /(p1::T, p2::Number) where T<:AbstractPolynomial
    return divide(convert(Polynomial, p1), convert(Polynomial, p2))
end


## Exponent
function powpoly(p::Polynomial, d::Int)
    if length(p) == 1
        expo, λ = collect(p)[1]

        expod = SortedDict{Variable, Degree}()
        for (var, deg) in expo.expo
            expod[var] = Degree(deg.explvar*d, deg.conjvar*d)
        end
        return Polynomial(SortedDict{Exponent, Number}(Exponent(expod)=>λ^d))
    else
        t_poly = Polynomial(1)
        for i=1:d
            t_poly = product(t_poly, p)
        end
        return t_poly
    end
end

function ^(p::T, d::Int) where T<:AbstractPolynomial
    return powpoly(convert(Polynomial, p), d)
end
