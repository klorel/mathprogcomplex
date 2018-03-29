
mutable struct RelaxationContext
    ismultiordered
    issparse
    leveragesymmetries      # Check if equations have certain type of symmetry, to set afterwards some moments to 0 for example
    hierarchykind           # :Complex or :Real
    renamevars              # Replace variables with by shorter named ones
    di
    ki
end


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
    # NOTE: on ne pourra peut-être pas couper au Tuple{Expo, Expo}...
    expo2int::Dict{Exponent, Int}   # Carte de la matrice des coefficients d-ki
    int2expo::Dict{Int, Exponent}
    msize::Int                      # size of the matrix
end


const SDPBody = Dict{Tuple{String, String, Exponent, Exponent}, Dict{Tuple{Exponent, Exponent}, Complex128}}
const SDPRhs = Dict{Tuple{Exponent, Exponent}, Complex128}