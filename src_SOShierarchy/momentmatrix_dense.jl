import Base.prod!

"""
    MomentMatrix(mm, vars, order)

    Store a moment or localizing matrix of size `order`, corresponding to the `vars` variables in the `mm` dictionnary.
    **Note** that the matrix is indexed by a tuple of exponents, *the first of which contains only conjugated variables*, et second only real ones.
"""
mutable struct MomentMatrix
    mm::Dict{Tuple{Exponent, Exponent}, AbstractPolynomial}
    vars::Set{Variable}
    order::Int
end

"""
    mm = MomentMatrix(vars::Set{Variable}, d)

    Build the moment matrix corresponding to the moment of degree up to `d` of the `vars` polynomial algebra.
"""
function MomentMatrix(vars::Set{Variable}, d::Int)
    mm = Dict{Tuple{Exponent, Exponent}, AbstractPolynomial}()
    realexpos = compute_exponents(vars, d)
    conjexpos = compute_exponents(vars, d, compute_conj=true)
    for cexp in conjexpos
        for rexp in realexpos
            mm[(cexp, rexp)] = cexp*rexp
        end
    end
    return MomentMatrix(mm, copy(vars), d)
end

function copy(mm::MomentMatrix)
    return MomentMatrix(copy(mm.mm), mm.vars, mm.order)
end

function print(io::IO, mm::MomentMatrix)
    for (key, val) in mm.mm
        println(io, "($(key[1]), $(key[2])) ⟶  $val")
    end
end


##########################
# Moment matrix algebra
function prod!(p::T, mm::MomentMatrix) where T<:AbstractPolynomial
    for (key, val) in mm.mm
        mm.mm[key] = p * mm.mm[key]
    end
end
function prod!(λ::Number, mm::MomentMatrix)
    for (key, val) in mm.mm
        mm.mm[key] = λ*val
    end
end

function *(p::T, mm::MomentMatrix) where T<:AbstractPolynomial
    mm_copy = copy(mm)
    prod!(p, mm_copy)
    return mm_copy
end
*(mm::MomentMatrix, p::T) where T<:AbstractPolynomial = p*mm

function *(λ::Number, mm::MomentMatrix)
    mm_copy = copy(mm)
    prod!(λ, mm_copy)
    return mm_copy
end
*(mm::MomentMatrix, λ::Number) = λ*mm

# function add!(mm1::MomentMatrix, mm2::MomentMatrix)
#     (mm1.order == mm2.order) || error("add!(): Different orders")
#     (mm1.vars == mm2.vars) || error("add!(): Different var sets")
#     for key in keys(mm2)
#         add!(mm1[key], mm2[key])
#     end
# end

#######################################
# Conversion to B base

"""
    MomentMatrixBasis(basis, expo2int, int2expo, msize)

    Store the matrix coefficients of the moment variable decomposition of the moment matrix.

    Arguments:
    - basis::Dict{Exponent, AbstractMatrix} : matrice correspondant au moment clé,
    - expo2int : carte donnant les coordonnées de l'exposant clé dans la matrice des moments initiale,
    - int2expo : carte donnant l'xposant correspondant aux coordonnées clé dans la matrice des moments initiale,
    - msize: taille de la matrice des moments initiale (ordre d-k).
"""
mutable struct MomentMatrixBasis
    basis::Dict{Exponent, AbstractMatrix} # Les exposants sont de degré inférieur à d
    # NOTE: on ne pourra peu-être pas couper au Tuple{Expo, Expo}...
    expo2int::Dict{Exponent, Int}   # Carte de la matrice des coefficients d-ki
    int2expo::Dict{Int, Exponent}
    msize::Int                      # size of the matrix
end

function MomentMatrixBasis(vars, d, k)
    expo2int = Dict{Exponent, Int}()
    int2expo = Dict{Int, Exponent}()
    i = 1
    for expo in sort(collect(compute_exponents(vars, d-k)))
        expo2int[expo] = i
        int2expo[i] = expo
        i += 1
    end

    return MomentMatrixBasis(Dict{Exponent, AbstractMatrix}(), expo2int, int2expo, length(expo2int))
end


"""
    mmb = convertMMtobase(mm::MomentMatrix, d, k)

    Compute the projection of the `mm` localizing matrix of the `k`-degree constraint on the order `d` moment basis.
    Yield a MomentMatrixBasis object with *dense* matrices.
"""
function convertMMtobase(mm::MomentMatrix, d, k)
    mmb = MomentMatrixBasis(mm.vars, d, k)
    for (key, poly) in mm.mm
        for (expo, λ) in poly
            if (expo.degree.explvar > d) || (expo.degree.conjvar > d)
                warn("convertMMtobase(): Found exponent of degree $(expo.degree) > $d ($expo, at $key of MM matrix)")
            end

            if !haskey(mmb.basis, expo)
                mmb.basis[expo] = Array{Complex64}(mmb.msize, mmb.msize) * 0
            end

            mmb.basis[expo][mmb.expo2int[conj(key[1])], mmb.expo2int[key[2]]] += λ
        end
    end
    return mmb
end
